# -*- coding: utf-8 -*-
#Steven Guinn
#04/21/2016
#University of Maryland Center for Environmental Science
#Appalachian Laboratory            http://www.al.umces.edu/ 
#301 689 7146
#301 Braddock Rd    Frostburg MD 215325

from modis_tools import *
import gdal_io
import os,sys
import time
#########################   CONSTANTS   ###########################################################
_NUM_PROC=16
_NUM_PARAMS=10
stime=time.clock()
qc_mask_dir ='/ws/Modis/wd/QC_MASK'
ndvi_dir  = '/ws/Modis/wd/NDVI'
mask_path = '/ws/Modis/wd/LULC_forest_mask.tif'
outpath= '/ws/Modis/wd/phenology.tif'

ndvi_list = os.listdir(ndvi_dir)
qc_mask_list = os.listdir(qc_mask_dir)
mask_image =  gdal_io.multi_band(mask_path,False)
## set the extent for a 101* 20 sub_region
#ullr={'upper':817157,'left':1001821,'lower':793982,'right':1006455} #this is the test area
## Prepare the output image
os.chdir(ndvi_dir)
_image=gdal_io.multi_band(ndvi_list[0],False)
phenology_image=gdal_io.raster_output(outpath,_image,IS_LANDSAT=True)
# phenology_image.extent.update(ullr)
# gt=phenology_image.extent.geotransform
# ngt=(ullr['left'],gt[1],gt[2],ullr['upper'],gt[4],gt[5])
# phenology_image.extent.geotransform=ngt
# phenology_image.extent.window =_image.extent.window
phenology_image.extent.describe()
_image =None
## prep the ndvi vegetation images
ndvi_image_list = list()
for image_path in ndvi_list:
	image=gdal_io.multi_band(image_path,True)
	ndvi_image_list.append(image)

ndvi_image_list.sort(key=lambda x: x.meta1['year'])
ndvi_image_list.sort(key=lambda x: x.meta1['doy'])
## 	create the master ndvi shared array
veg_list=list()
doy_list=list()

for image in ndvi_image_list:
	#image.extent.update(ullr)
	image.get_array(1)
	veg_list.append(image.array)
	doy_list.append(int(image.meta1['doy']))
	image.dataset = None

veg_array=np.dstack(veg_list)
veg_list =None
print "created veg array"
doy_array=np.array(doy_list)		##$& may have to reform this
print "created doy array"
ndvi_image_list=None
	
# prep the qc mask  images
os.chdir(qc_mask_dir)
qc_image_list = list()

for image_path in qc_mask_list:
	image=gdal_io.multi_band(image_path, True)
	qc_image_list.append(image)
qc_image_list.sort(key=lambda x: x.meta1['year'])
qc_image_list.sort(key=lambda x: x.meta1['doy'])
qc_list = list()

for image in qc_image_list:
	#image.extent.update(ullr)
	image.get_array(1)
	qc_list.append(image.array)
	image.dataset=None

qc_array=np.dstack(qc_list)
qc_list=None
print "created qc array"
qc_image_list=None


# get the mask array
#mask_image.extent.update(ullr)
mask_image.get_array(1)
mask_array = mask_image.array
mask_image.extent.describe()

sh_time=time.clock()
########################### Create shared arrays ###########################
#sh_mask_arr
sh_mask_arr_base = mp.Array(ctypes.c_short,mask_array.ravel(order = 'A'))
sh_mask_arr  = np.ctypeslib.as_array(sh_mask_arr_base.get_obj()) 
sh_mask_arr.shape=mask_array.shape
print "shared mask"
esh_time=time.clock()
print 'shared array time ',esh_time-sh_time, ' seconds'
sh_time=time.clock()
#sh_veg_arr
sh_veg_arr_base = mp.Array(ctypes.c_short,veg_array.ravel(order = 'A'))
sh_veg_arr  = np.ctypeslib.as_array(sh_veg_arr_base.get_obj()) 
sh_veg_arr.shape = veg_array.shape
print "shared veg"
esh_time=time.clock()
print 'shared array time ',esh_time-sh_time, ' seconds'
sh_time=time.clock()
#sh_qc_arr
sh_qc_arr_base = mp.Array(ctypes.c_short,qc_array.ravel(order = 'A'))
sh_qc_arr  = np.ctypeslib.as_array(sh_qc_arr_base.get_obj()) 
sh_qc_arr.shape= qc_array.shape
print "shared qc"
esh_time=time.clock()
print 'shared array time ',esh_time-sh_time, ' seconds'
sh_time=time.clock()
#sh_doy_arr
sh_doy_arr_base = mp.Array(ctypes.c_short,doy_array.ravel(order = 'A'))
sh_doy_arr  = np.ctypeslib.as_array(sh_doy_arr_base.get_obj()) 
sh_doy_arr.shape= doy_array.shape
print "shared doy"
esh_time=time.clock()
print 'shared array time ',esh_time-sh_time, ' seconds'


#calculate indes to slice up the data into chunks
num_rows = mask_image.extent.window['ywin']
num_chunks=100
mod_exp=divmod(num_rows,num_chunks)

#slice and dice
slice_list=list()
for i in range(num_chunks):
	start=(i*mod_exp[0])
	end=((i+1)*mod_exp[0])
	if (i+1) == num_chunks:
		end=num_rows
	slice_list.append((start,end))

test_list=list()			##################!!!!!!!!! Test CODE !!!!!!!!!!!!!!!!!!!!
	
############################################################################################################################
############################################  Multiprocessing Function Definitions  ########################################
############################################################################################################################
def shared_array_test(_slc,def_param = sh_veg_arr, ):		#_slc is a tuple with the start/stoprow in slice notation, _shape is the _arrary 2D shape tuple
	_ndarray = def_param[_slc[0]:_slc[1],:]
	return _ndarray[:,-1,:].sum()

def Phenology(_slc,_veg_arr=sh_veg_arr,_qc_arr=sh_qc_arr,_doy_arr=sh_doy_arr,_mask_arr=sh_mask_arr):
	#Extract the array chunck for this child process
	veg_chunk =( _veg_arr[_slc[0]:_slc[1],:,:]) * SCALEF		#convert the veg array chunk to float
	qc_chunk =_qc_arr[_slc[0]:_slc[1],:,:]
	mask_chunk =_mask_arr[_slc[0]:_slc[1],:]
	#Create output array
	_out_arr=np.zeros((veg_chunk.shape[0],_veg_arr.shape[1],_NUM_PARAMS),'d')
	row=0
	while row <  veg_chunk.shape[0]:
		print "row",row
		col=0
		while col < veg_chunk.shape[1]:
			vslc=veg_chunk[row,col,:]
			mslc=qc_chunk[row,col,:]
			dmat=np.matrix(np.compress(mslc,_doy_arr))
			vmat=np.matrix(np.compress(mslc,vslc))
			if mask_chunk[row,col]==1:
				_out_arr[row,col,:]=sigmoid(vmat,dmat).reshape((_NUM_PARAMS,))
			else:
				_out_arr[row,col,:]=FLOAT32_NODATA
			col +=1
		row+=1
	veg_chunk=None
	qc_chunk=None
	mask_chunk=None
	#print "chunk number:",chunk_num,"finished"
	return _out_arr
	
############################################################################################################################
############################################  Main Program Execution  ######################################################
############################################################################################################################
print "lets get busy !"
pool = mp.Pool(_NUM_PROC)
#result = pool.map(shared_array_test, slice_list)
result = pool.map(Phenology, slice_list)
pool.close()
pool.join()
# # # print np.array(result).sum()
# # # print veg_array[:,-1,:].sum()
out_arr = np.vstack(result)
# # row = 2
# # col = 18
# # vslc=veg_array[row,col,:]
# # mslc=qc_array[row,col,:]
# # dmat=np.matrix(np.compress(mslc,doy_array))
# # vmat=np.matrix(np.compress(mslc,vslc))
# # out_arr=np.zeros((veg_array.shape[0],veg_array.shape[1],_NUM_PARAMS),'d')
# # out_arr[row,col,:]=sigmoid(vmat,dmat).reshape((_NUM_PARAMS,))
phenology_image.datatype=6		#set to float32	
phenology_image.set_bandlist(range(1,_NUM_PARAMS+1))	#set the number of output bands (number of parameters + 1 for range
phenology_image.Create_Dataset()		#create the GDAL Dataset Object
phenology_image.get_multiband_array()
phenology_image.array=phenology_image.array.astype('d')
phenology_image.array=out_arr
phenology_image.write_array()
phenology_image=None
etime=time.clock()
print 'Total processing time on ',_NUM_PROC,' of processors is',(etime-stime)/60, ' minutes'