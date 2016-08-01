# -*- coding: utf-8 -*-
#Steven Guinn
#01/04/2016
#University of Maryland Center for Environmental Science
#Appalachian Laboratory            http://www.al.umces.edu/ 
#301 689 7146
#301 Braddock Rd    Frostburg MD 215325

import os,sys
from osgeo import gdal
from osgeo import gdal_array
from osgeo.gdal_array import *
from osgeo.gdalconst import *
import tarfile as tar
import gdal_io
import numpy as np
import numpy.ma as ma
import math
import multiprocessing
from copy import deepcopy

########################################  Constants ############################################
TORAD = math.pi/180
INT16_NODATA = -32768
FLOAT32_NODATA = -3.402823466e+38
ACCESS_MODE = 0777
SCALEF = .0001
########################################  Classes   ############################################

######################################## Functions   ###########################################
def cleanMakeDir(out_dir,mode=0777):
	if os.path.exists(out_dir):
		print "Output Directory exist!"
		removeall(out_dir)   #this erases the contents of the directory
		print "Cleaned Output Directory!"
	else:
		os.mkdir(out_dir)
		os.chmod(out_dir,mode)
		print "Created Output Directory!"
def rmgeneric(path, __func__):
    try:
        __func__(path)
        print 'Removed ', path
    except OSError, (errno, strerror):
        print "Error" % {'path' : path, 'error': strerror }           
def removeall(path):
    if not os.path.isdir(path):
        return    
    files=os.listdir(path)
    for x in files:
        fullpath=os.path.join(path, x)
        if os.path.isfile(fullpath):
            f=os.remove
            rmgeneric(fullpath, f)
        elif os.path.isdir(fullpath):
            removeall(fullpath)
            f=os.rmdir
            rmgeneric(fullpath, f)   
def listFiles(dir):
    rootdir = dir
    for root, subFolders, files in os.walk(rootdir):
        for file in files:
            yield os.path.join(root,file)
    return
def make_tarfile(outpath,file_list,MODE="w:gz"):
	tarfile = tar.open(outpath,MODE)
	for file in file_list:
		print "Adding to "+os.path.basename(outpath)+" "+os.path.basename(file)
		tarfile.add(file)
	tarfile.close()
#extract_eos_archive()	A function that extracts a landsat surface reflection archive file into a multiband geotiff. Extent and bands can be set, or entire image is processed
def extract_eos_archive(archive,outdir,ullr=None,bandList=None,PRODUCT='LNDSR'):
	temp_dir=os.path.join(os.path.dirname(archive),'.temp')
	cleanMakeDir(temp_dir)
	if tar.is_tarfile(archive):
		if PRODUCT=='LNDSR':
			tfile = tar.open(archive, 'r:gz')
			tfile.extractall(temp_dir)
			temp_dirlist = os.listdir(temp_dir)
			for d in temp_dirlist:
				path=(os.path.join(temp_dir,d))
				if os.path.splitext(path)[1]=='.hdf':
					eos_hdf_file=path
			image_eos=gdal_io.eos_hdf(eos_hdf_file)
			if ullr:
				image_eos.extent.update(ullr)
			if bandList:
				image_eos.set_bandlist(bandList)
			out_name=image_eos.meta1['product']+'.'+image_eos.meta1['scene_name']+'.tif'
			outpath=os.path.join(outdir,out_name)
			print outpath
			image_eos.write_array(outpath)
			#image_eos=None
			return image_eos
	else:
		print "ERROR: Not a TARFILE"
def SMA(params):
	#process the parameters
	inpath=params[0]	#str
	outpath=params[1]	#str
	end_members_textfile=params[2]	#str
	num_endmembers=params[3]	#int
	atmo_scattering=params[4]	#int the endmember column to use for atmospheric correction
	#get the endmembers and process them
	em_array=extract_endmembers(end_members_textfile,num_endmembers)
	em=np.asmatrix(em_array)
	emt=em.T
	emtem=emt*em
	emtem_inv=emtem.I
	#open the input image dataset and get the input array
	image=gdal_io.multi_band(inpath)
	image.get_multiband_array()
	#check for atmospheric scattering
	if atmo_scattering:
		a_corr=em_array[:,atmo_scattering-1]
		bcount=0
		while bcount < len(a_corr):
			image.array[:,:,bcount]=image.array[:,:,bcount]+a_corr[bcount]
			bcount+=1
	#create the output object useing the input as a template
	sma_image=gdal_io.raster_output(outpath,image)
	sma_image.datatype=6		#set to float32	
	sma_image.set_bandlist(range(1,num_endmembers+3))	#set the number of output bands (number of endmembers + 1 for RMSE calculations +1 mask
	sma_image.Create_Dataset()		#create the GDAL Dataset Object
	#create the output array to store the results
	sma_image.get_multiband_array()
	sma_image.array=sma_image.array.astype('d')	
	mask=((image.array[:,:,-1]<=0).astype('b'))		#create the mask
	image.array=image.array[:,:,:-1]	#remove the FMASK band from input image
	#unmix the pixels one at a time
	rmse_avg=0
	row=0
	while row < image.array.shape[0]:
		#print int((float(row)/image.array.shape[0])*100),"% of image processed"
		col=0
		while col < image.array.shape[1]:
			if mask[row][col]:
				pix=image.array[row,col,:]
				pixmat=np.asmatrix(pix)
				B=pixmat.T
				Bemt=emt*B
				result=emtem_inv*Bemt
				#rmse calculation
				r1=em*result
				r2=B-r1
				r3=r2*r2.T
				r4=r3.diagonal().sum()
				rmse=math.sqrt(r4/6)
				rmse_avg=rmse_avg+rmse
				#output results
				rmsem=np.asmatrix(rmse)
				out=np.vstack((result,rmsem,mask[row][col]))
				sma_image.array[row,col,:]=np.resize(np.array(out),(6,))
				col+=1
			else:
				col +=1
		row+=1
	sma_image.write_array()
	sma_image=None
	image=None
	return 1
def extract_endmembers(file,num_endmembers):
	text=open(file,'r')
	textlist=text.readlines()
	em_array=np.zeros((len(textlist),num_endmembers),'d')
	row=0
	for line in textlist:
		words=line.split(',')
		col=0
		for word in words:
			if col==0:
				pass
			else:
				em_array[row][col-1]=float(word)
			col+=1
		row+=1
	return em_array
def getCosTheta(slope,aspect,suna,sunz):
    if slope.shape==aspect.shape:
        #return np.cos(TORAD*sunz)*np.cos(TORAD*slope)+np.sin(TORAD*sunz)*np.sin(TORAD*slope)*np.cos(TORAD*(aspect-suna))
		return np.cos(TORAD*sunz)*np.cos(TORAD*slope)+np.sin(TORAD*sunz)*np.sin(TORAD*slope)*np.cos(((TORAD*aspect)-(TORAD*suna)))
    else:
        #return -9999
		print "Shape Error!!!!"
	
def SCS(params):		#assumes a image stack with the last band being a mask band with 0= unmasked cells
	image_path=params[0]
	out_path=params[1]
	slope_path=params[2]
	aspect_path=params[3]
	azimuth=params[4]
	zenith=params[5]
	image=gdal_io.multi_band(image_path)
	try:
		ullr=params[6]
	except IndexError:
		ullr=None
	if ullr==None:
		ullr=image.extent.getULLR()
	#open datasets 
	image.extent.update(ullr)
	image.get_multiband_array()
	#ullr=image.extent.getULLR()
	#slope preprocessing
	slope=gdal_io.multi_band(slope_path,None)
	slope.extent.update(ullr)
	slope.get_array()
	slope_arr=slope.array	
	# aspect
	aspect=gdal_io.multi_band(aspect_path,None)
	aspect.extent.update(ullr)
	aspect.get_array()
	aspect_arr=aspect.array
	#create the output object useing the input as a template
	scs_image=scs_image=gdal_io.raster_output(out_path,image)
	gt=scs_image.extent.geotransform
	ngt=(ullr['left'],gt[1],gt[2],ullr['upper'],gt[4],gt[5])
	scs_image.extent.geotransform=ngt
	scs_image.extent.update(ullr)
	scs_image.set_bandlist(range(1,image.numbands+1))	#set the number of output bands (input bands +1 for the final mask)
	#scs_image.set_bandlist(range(1,image.numbands+2))   #write out cosine theta........
	scs_image.Create_Dataset()	
	mask=((image.array[:,:,-1] > 0).astype(np.int16))		#create the mask
	image.array=image.array[:,:,:-1]	#remove the FMASK band from input image array
	slope_arr=slope_arr/math.pi		#new slope correction
	cos_theta_array=getCosTheta((slope_arr), aspect_arr, azimuth, zenith)
	#slope_arr=slope_arr/math.pi		#new slope correction
	cos_theta_maskout=(cos_theta_array < numpy.cos(TORAD*85.0)).astype(np.int16)
	#print "images shape",image.array.shape
	mask=((mask+cos_theta_maskout) > 0).astype(np.int16)
	inv_mask=(mask==0).astype(np.int16)
	bcount=0
	band_array_list=list()
	while bcount < image.array.shape[2]:
		addback_arr=(image.array[:,:,bcount]*mask).astype(np.int16)
		band=numpy.ma.masked_array(image.array[:,:,bcount],mask)
		#scs_band=(band * (numpy.cos(TORAD*slope_arr)/cos_theta_array)).astype(np.int16)
		#scs_band=((band * numpy.cos(TORAD*slope_arr))/cos_theta_array).astype(np.int16)
		scs_band=(band * numpy.cos(TORAD*slope_arr)* numpy.cos(TORAD*zenith)/cos_theta_array).astype(np.int16)
		#scs_band=(band * (cos_theta_array/(numpy.cos(TORAD*slope_arr)* numpy.cos(TORAD*zenith)))).astype(np.int16)
		#scs_band= band*cos_theta_array#/math.pi
		band_array_list.append((scs_band.data*inv_mask)+addback_arr)
		bcount +=1
	band_array_list.append(mask)
	#band_array_list.append(cos_theta_array*10000)
	scs_image.array=np.dstack(band_array_list)
	scs_image.write_array()
	#clean it UP!
	image=None
	scs_image=None
	slope=None
	aspect=None
	return 1
def ReadMetaData(path,PRODUCT='LNDSR'):
	text=open(path,'r')
	textlist=text.readlines()
	subkeyline=textlist.pop(0)[:-1]
	subkeylist=subkeyline.split(',')
	subkeylist=subkeylist[1:]
	meta_dict=dict()
	for line in textlist:
		line_dict=dict()
		line_items=line.split(',')
		key=line_items.pop(0)[:16]
		lcount=0
		for subkey in subkeylist:
			line_dict[subkey.strip()]=line_items[lcount].strip()
			lcount+=1
		meta_dict[key]=line_dict
	return meta_dict
def fwdlogit(m,t):
    s1=np.exp((m[0,2]-t)/m[0,3])
    s2=np.exp((m[0,4]-t)/m[0,5])
    return m[0,0]+np.multiply((m[0,1]-(np.multiply(t,m[0,6]))),((1.0/(1.0+s1))-(1.0/(1+s2))))

def l1_solution(gv,mest,t,wd):
    dataresiduals = np.sqrt(np.power((gv - fwdlogit(mest,t)),2))
    ident=np.mat(np.identity(t.shape[1]))
    Rd = ident/dataresiduals
    wg = wd.T * Rd * wd
    return wg

def calcG(m,t):
	d3=1+np.exp((m[0,2]-t)/m[0,3])
	d5=1+np.exp((m[0,4]-t)/m[0,5])   
	#partial derivatives of new equations
	g1 = np.mat(np.ones(t.shape[1])).T
	g2 = ((1/d3)-(1/d5)).T
	g3 = (np.multiply(np.multiply((1/m[0,3]) , np.exp((m[0,2]-t)/m[0,3])) , (np.multiply(m[0,6],t) -m[0,1]))/np.power(d3,2)).T
	g4 = (np.multiply((1/np.power(m[0,3],2)) , np.multiply((m[0,2]-t) , np.multiply(np.exp((m[0,2]-t)/m[0,3]) , ((m[0,1] - np.multiply(m[0,6],t))/ np.power(d3,2)))))).T
	g5 = (np.multiply((1/m[0,5]) , np.multiply(np.exp((m[0,4]-t)/m[0,5]) , (m[0,1] - np.multiply(m[0,6],t))/ np.power(d5,2))) ).T
	g6 = (np.multiply((1/np.power(m[0,5],2)) , np.multiply((m[0,4]-t) , np.multiply(np.exp((m[0,4]-t)/m[0,5]) , (np.multiply(m[0,6],t) - m[0,1])))) / np.power(d5,2)).T
	g7 = ((t/d5)-(t/d3)).T
	G=np.hstack((g1,g2,g3,g4,g5,g6,g7))
	return G

def misfitSmall(misfit,dresid):
    result = math.sqrt(2/np.pi)*dresid.shape[1]
    return misfit <= result
    
def misfitDelta(misfit,oldmisfit):
    result=math.fabs(misfit-oldmisfit)/oldmisfit
    return result < 0.005

def sigmoid(vegmat,doymat):   
	mdefault=np.mat('.00001;.00001;.001;.00001;.001;.00001; 0.000001').T	##create the initial parameters
	mest=np.mat('0.07; 1.3; 125.0; 7.0; 300.0; 7.0; 0.003').T					#first estimate of the parameters
	mref=mest   
	cin=np.mat('0.05; 0.5; 30.0; 2.0; 30.0; 2.0; 0.002').T						#set up the cin according to the num parameters (damping matrix)
	#Parameter weights
	ident=np.mat(np.identity(7))
	Wm=ident/cin
	C=Wm*Wm
	#data weights
	ident=np.identity(doymat.shape[1])    #this is a proxy for MatLabs diagonal matrix function
	rms=np.zeros(doymat.shape[1])
	rms+=.05
	Wd=np.asmatrix(ident/rms)
	Wg=Wd*Wd  
	dresid = (vegmat-fwdlogit(mest,doymat))
	#set up the stopping values
	misfit = np.sqrt(np.power(Wd*dresid.T,2)).sum()
	oldmisfit=2000
	pert=0
	##this is where the start of the loop goes
	try:
		while (pert < 2) or ((pert < 101) and (not misfitSmall(misfit,dresid)) and (not misfitDelta(misfit,oldmisfit))):
			G=calcG(mest,doymat)
			if G == None:
				return mdefault.T
			GtWgG = G.T*Wg*G
			GtWgGdmp = GtWgG+C        
			dresid = (vegmat-fwdlogit(mest,doymat))
			GtWg = G.T*Wg        
			GtWgdresid = (GtWg*dresid.T)+(C*(mref-mest).T)
			mpert = GtWgGdmp.I * GtWgdresid        
			mest = mest+mpert.T     
			Wg=l1_solution(vegmat,mest,doymat,Wd)
			oldmisfit=misfit
			misfit = np.sqrt(np.power(Wd*dresid.T,2)).sum()
			pert +=1
		return np.hstack([mest,np.array([misfit]).reshape((1,1)),np.array([pert]).reshape((1,1)),np.array([doymat.shape[1]]).reshape((1,1))])
	except numpy.linalg.linalg.LinAlgError:
		print "whoops"
		return np.mat('0.0;0.0;0.0;0.0;0.0;0.0;0.0; 0.0; 0.0; 0.0').T
def Phenology(parameters):
	veg_list=parameters[0]
	ullr=parameters[1]
	num_parameters=parameters[2]
	image_list=list()
	for image_path in veg_list:
		image=gdal_io.multi_band(image_path)
		image_list.append(image)
		#print image.meta1['scene_name'],image.meta1['doy']
	image_list.sort(key=lambda x: x.meta1['doy'])
	mask_list=list()
	veg_list=list()
	doy_list=list()
	for image in image_list:
		#print image.meta1['scene_name'],image.meta1['doy']
		image.extent.update(ullr)
		image.get_array(6)
		mask_list.append(image.array)
		image.get_array(3)
		veg_list.append(image.array)
		doy_list.append(float(image.meta1['doy']))
	mask_array=np.dstack(mask_list)
	veg_array=np.dstack(veg_list)
	doy_array=np.array(doy_list)
	out_arr=np.zeros((veg_array.shape[0],veg_array.shape[1],num_parameters),'d')     #create the output array for each process
	row=0
	while row < veg_array.shape[0]:
		print "row",row
		col=0
		while col < veg_array.shape[1]:
			vslc=veg_array[row,col,:]
			mslc=mask_array[row,col,:]
			dmat=np.matrix(np.compress(mslc,doy_array))
			vmat=np.matrix(np.compress(mslc,vslc))
			if mslc.sum()!=0:
				out_arr[row,col,:]=sigmoid(vmat,dmat).reshape((num_parameters,))
			else:
				out_arr[row,col,:]=FLOAT32_NODATA
			col +=1
		row+=1
	veg_array=None
	mask_array=None
	doy_array=None
	#print "chunk number:",chunk_num,"finished"
	return out_arr

def NDVI(rast_1_path,rast_2_path,mask_path,outpath):
	############################################################ I/O setup, read data sets and process into masked arrays
	#open mask band
	try:
		image_mask=gdal_io.multi_band(mask_path,False)
		image_mask.get_multiband_array()
		mask = image_mask.array   !=4096
		image_mask = None		#free up the memory 
		#open red band
		image_red=gdal_io.multi_band(rast_1_path,False)
		image_red.get_multiband_array()
		#get a masked array from the red data object
		red_array=ma.masked_array(deepcopy(image_red.array),mask)
		#create the output object using the red band as a template
		ndvi_image=gdal_io.raster_output(outpath,image_red,False)
		ndvi_image.datatype=3		#set to int16
		ndvi_image.Create_Dataset()
		ndvi_image.get_multiband_array()
		image_red=None		#free up the memory 
		#open nir band
		image_nir=gdal_io.multi_band(rast_2_path,False)
		image_nir.get_multiband_array()
		nir_array=ma.masked_array(deepcopy(image_nir.array),mask)
		image_nir=None		#free up the memory
		########################################################### NDVI Processing
		ndvi_arr=((((nir_array-red_array)/(nir_array+red_array).astype('float16'))/SCALEF).astype('int16'))
		ma.set_fill_value(ndvi_arr,INT16_NODATA)
		ndvi_image.array =ndvi_arr.data
		ndvi_image.write_array()
		########################################################### Clean Up
		ndvi_image=None
		nir_array=None
		red_array=None
		return 0
	except Exception:
		return -1

def Modis_Mask_Extract(rast_1_path,outpath,MASK=4096):
	try:
		image=gdal_io.multi_band(rast_1_path,False)
		image.get_multiband_array()
		mask = image.array   ==  MASK
		mask_image=gdal_io.raster_output(outpath,image,False)
		mask_image.Create_Dataset()
		mask_image.get_multiband_array()
		mask_image.array = mask
		mask_image.write_array()
		mask_image=None
		return 0
	except Exception:
		print "whoops", rast_1_path,"failed"
		return -1
	
	













