# -*- coding: utf-8 -*-
#Steven Guinn
#02/26/2014
#University of Maryland Center for Environmental Science
#Appalachian Laboratory            http://www.al.umces.edu/ 
#301 689 7146
#301 Braddock Rd    Frostburg MD 215325

import os
import sys
from osgeo import gdal
from osgeo import gdal_array
from osgeo.gdalconst import *
import numpy as np
from copy import deepcopy
# Need to define object structure such as extent and window
# window={'startx':None,'starty':None,'xwin':None,'ywin':None}
# extent={'upper':None,'left':None,'lower':None,'right':None}
# meta1={'sensor':None,'path':None,'row':None,'year':None,'doy':None,'product':None}
########################################  Constants ############################################
NP2GDAL_CONVERSION = { "uint8": 0, "int8": 1, "uint16": 2, "int16": 3, "uint32": 4, "int32": 5, "float32": 6, "float64": 7, "complex64": 10, "complex128": 11, }
DEFAULT_GT = (0.0,1.0,0.0,0.0,0.0,1.0)
INT16_NODATA = -32768
FLOAT32_NODATA = -3.402823466e+38
ACCESS_MODE=0777
########################################  Classes   ############################################
class gdal_file(object):
	def __init__(self,path):
		self.path=path
		self.dataset=gdal.Open(self.path)
		self.driver=self.dataset.GetDriver()
		self.projection=self.dataset.GetProjectionRef()
		self.meta1={'sensor':None,'path':None,'row':None,'year':None,'doy':None,'product':None,'scene_name':None}
class multi_band(gdal_file):	#A class object that wraps a multiband gdal raster object. i.e ENVI,IMG,GeoTIFF
	def __init__(self,path,IS_LANDSAT=True):
		if os.path.isfile(path):
			super(multi_band,self).__init__(path)
			self.extent=Extent(self.dataset)
			self.numbands=self.dataset.RasterCount
			self.set_bandlist(BandList=None)	#this bandlist is a direct acces to the GDAL bandlist convention, [1:numbands]
			self.datatype=self.dataset.GetRasterBand(1).DataType
			self.array=None
			if IS_LANDSAT:
				self.parse_filename_lndsat()
			else:
				try:
					self.parse_filename_modis()
				except Exception:
					pass
		else:
			self.path=path
			self.dataset=None
			self.driver=None
			self.projection=None
			self.meta1={'sensor':None,'path':None,'row':None,'year':None,'doy':None,'product':None,'scene_name':None}
			self.extent=None
			self.numbands=1
			self.set_bandlist(BandList=None)
			self.datatype=None
			self.datatype=None
			self.array=None
			if IS_LANDSAT:
				self.parse_filename_lndsat()
			else:
				try:
					self.parse_filename_modis()
				except Exception:
					pass
	def parse_filename_lndsat(self):
		fname=os.path.basename(os.path.basename(self.path))
		words=os.path.splitext(fname)[0].split('.')
		if len(words) > 1:
			self.meta1['scene_name']=words[1][0:16]
			self.meta1['product']=words[0]
			self.meta1['sensor']=words[1][0:3]
			self.meta1['path']=int(words[1][3:6])
			self.meta1['row']=int(words[1][6:9])
			self.meta1['year']=int(words[1][9:13])
			self.meta1['doy']=int(words[1][13:16])
		else:
			self.meta1['scene_name']=words[0][0:16]
			self.meta1['product']='unknown'
			self.meta1['sensor']=words[0][0:3]
			self.meta1['path']=int(words[0][3:6])
			self.meta1['row']=int(words[0][6:9])
			self.meta1['year']=int(words[0][9:13])
			self.meta1['doy']=int(words[0][13:16])
	def parse_filename_modis(self):
		fname=os.path.basename(os.path.basename(self.path))
		words=os.path.splitext(fname)[0].split('_')
		self.meta1['product']=words[0]
		self.meta1['sensor']='MODIS'
		self.meta1['year']=int(words[1][0:4])
		self.meta1['doy']=int(words[1][4:])
	def get_array(self,Bandnumber=None):		# Bandnumber is base 1 (gdal format)
		if Bandnumber:
			self.array =self.dataset.GetRasterBand(Bandnumber).ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin'],self.extent.window['ywin'])
		else:
			self.array =self.dataset.GetRasterBand(1).ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin'],self.extent.window['ywin'])
	def get_multiband_array(self,BandList=None):
		if BandList:
			blist=BandList
		else:
			blist=self.bandlist
		array_list=list()
		for b in blist:
			array_list.append(self.dataset.GetRasterBand(b).ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin'],self.extent.window['ywin']))
		self.array = np.dstack(array_list)
	def get_pixel_array(self,xy_tuple,BandList=None):
		ullr={'upper':xy_tuple[1],'left':xy_tuple[0],'lower':xy_tuple[1],'right':xy_tuple[0]}
		self.extent.update(ullr)
		if BandList:
			blist=BandList
		else:
			blist=self.bandlist
		array_list=list()
		for b in blist:
			array_list.append(self.dataset.GetRasterBand(b).ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin']+1,self.extent.window['ywin']+1))
		self.extent.reset()
		return np.dstack(array_list).flatten()
	def set_bandlist(self,BandList=None):
		if BandList:
			self.bandlist=BandList
		else:
			self.bandlist=range(1,self.numbands+1)
class eos_hdf(gdal_file): # A class object that wraps a EOS HDF file format containing a Landsat Surface reflectance product derived from the LEDAPS processing framework
	def __init__(self,path):
		super(eos_hdf,self).__init__(path)
		self.subdatasets=dict()		
		self.get_sub_datasets()
		self.extent=Extent(self.subdatasets['b1']['DS'])
		self.parse_filename_lndsat()
		self.set_bandlist()		#this bandlist is a list of dictionary keys
		self.datatype=self.subdatasets['b1']['DS'].GetRasterBand(1).DataType
		self.update_projection()
		self.array=None
	def get_sub_datasets(self):		#retrieves and opens subdataset from the EOS HDF file as GDAL dataset objects
		subdatsets_list=self.dataset.GetSubDatasets()
		self.numbands=len(subdatsets_list)
		bcount=0
		while bcount < self.numbands:
			bandkey='b'+str(bcount+1)
			self.subdatasets[bandkey]={'name':subdatsets_list[bcount][0],'description':subdatsets_list[bcount][1],'DS':gdal.Open(subdatsets_list[bcount][0])}
			bcount+=1
	def parse_filename_lndsat(self):
		fname=os.path.basename(self.path)
		words=fname.split('.')
		self.meta1['scene_name']=words[1][0:16]
		self.meta1['product']=words[0]
		self.meta1['sensor']=words[1][0:3]
		self.meta1['path']=int(words[1][3:6])
		self.meta1['row']=int(words[1][6:9])
		self.meta1['year']=int(words[1][9:13])
		self.meta1['doy']=int(words[1][13:16])
	def get_multiband_array(self,BandList=None):
		if BandList:
			blist=BandList
		else:
			blist=self.bandlist
		array_list=list()
		for band in blist:
			array_list.append(self.subdatasets[band]['DS'].ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin'],self.extent.window['ywin']))
		self.array = np.dstack(array_list)
	def set_bandlist(self,BandList=None):
		if BandList:
			self.bandlist=BandList
			self.numbands=len(BandList)
		else:
			bl=self.subdatasets.keys()
			bl.sort(key=lambda x: int(x[1:]))
			self.bandlist=bl
			self.numbands=len(bl)
	def write_array(self,outpath,format='GTiff'):
		#check to see if image object array exist
		if self.array==None:
			self.get_multiband_array(self.bandlist)
		#create the new GDAL Driver object based on image object extent and number of bands
		image_driver = gdal.GetDriverByName(format)
		projection = self.projection
		#update the geotransform based on new extent (if necessary)
		geotransform = Create_GeoTransform(self.extent)
		#retrieve the Dataset object from the driver based on image object extent and number of bands
		image_out_DS=image_driver.Create(outpath,self.extent.xsize,self.extent.ysize,self.numbands,self.datatype)
		image_out_DS.SetProjection(projection)
		image_out_DS.SetGeoTransform(geotransform)
		count=1
		while count <= image_out_DS.RasterCount:
			image_band_out=image_out_DS.GetRasterBand(count)
			image_band_out.SetNoDataValue(INT16_NODATA)
			image_band_out=gdal_array.BandWriteArray(image_band_out,self.array[:,:,count-1])
			count=count+1
		image_out_DS=None
		os.chmod(outpath,ACCESS_MODE)	
	def update_projection(self):
		self.projection=self.subdatasets['b1']['DS'].GetProjectionRef()
class multi_band_virt(gdal_file):
	def __init__(self,gdal_dataset):
		self.parent_path=gdal_dataset.path
		self.path=''
		self.driver = gdal.GetDriverByName( 'MEM' )
		self.numbands=gdal_dataset.numbands
		self.rasterXSize = gdal_dataset.rasterXSize 
		self.rasterYSize = gdal_dataset.rasterYSize
		self.get_sub_datasets(gdal_dataset)
		self.geotransform = gdal_dataset.geotransform
			
	def get_sub_datasets(self,gdal_dataset):
		self.subdatasets=dict(gdal_dataset.subdatasets)
		bcount=0
		while bcount < self.numbands:
			bandkey='b'+str(bcount+1)
			self.subdatasets[bandkey]['DS']=self.driver.CreateCopy(self.path,self.subdatasets[bandkey]['DS'])
			bcount+=1	

class Extent(object):		
	def __init__(self,gdal_object):
		self.geotransform=Get_GeoTransform(gdal_object)
		self.extent=Get_Extent(gdal_object)
		self.original_extent = dict(self.extent)
		self.get_window()
		self.xsize=self.window['xwin']
		self.ysize=self.window['ywin']
	def get_window(self):
		self.window={'startx':None,'starty':None,'xwin':None,'ywin':None}
		self.window['startx'] = GetRasterCoord((self.extent['left'],self.extent['upper']),self.geotransform)[0]
		self.window['starty'] = GetRasterCoord((self.extent['left'],self.extent['upper']),self.geotransform)[1] 
		self.window['xwin'] = (GetRasterCoord((self.extent['right'],self.extent['upper']),self.geotransform)[0] - self.window['startx'])
		self.window['ywin'] = (GetRasterCoord((self.extent['right'],self.extent['lower']),self.geotransform)[1] - self.window['starty'])
	def update(self,ullr):
		self.extent['upper']=ullr['upper']
		self.extent['left']=ullr['left']
		self.extent['lower']=ullr['lower']
		self.extent['right']=ullr['right']
		self.get_window()
		self.xsize=self.window['xwin']
		self.ysize=self.window['ywin']
	def reset(self):
		self.update(self.original_extent)
	def describe(self):
		for k in self.extent.keys():
			print k,self.extent[k]
		for k in self.window.keys():
			print k,self.window[k]
	def getULLR(self):		#returns a urrl formated dictionary object
		return {'upper':float(self.extent['upper']),'left':float(self.extent['left']),'lower':float(self.extent['lower']),'right':float(self.extent['right'])}
	def resetWindow(self):
		self.window['startx'] = 0
		self.window['starty'] = 0
		self.window['xwin'] = (GetRasterCoord((self.extent['right'],self.extent['upper']),self.geotransform)[0] - self.window['startx'])
		self.window['ywin'] = (GetRasterCoord((self.extent['right'],self.extent['lower']),self.geotransform)[1] - self.window['starty'])
class raster_output(object):
	def __init__(self,path,gdal_object,IS_LANDSAT=True):
		if isinstance(gdal_object,gdal.Dataset):
			self.path=path
			self.driver = gdal_object.GetDriver()
			self.projection = gdal_object.GetProjectionRef()
			self.numbands = gdal_object.RasterCount
			self.datatype=NP2GDAL_CONVERSION[gdal_object.GetRasterBand(1).DataType]
			self.extent=Extent(gdal_object)
			self.set_bandlist(BandList=None)
			self.array=None
			self.meta1={'sensor':None,'path':None,'row':None,'year':None,'doy':None,'product':None,'scene_name':None}
			if IS_LANDSAT:
				self.parse_filename_lndsat()
			#self.Create_Dataset()
		if isinstance(gdal_object,gdal_file):
			self.path=path
			self.driver = gdal_object.driver
			self.projection = gdal_object.projection
			self.extent = deepcopy(gdal_object.extent)
			self.numbands = gdal_object.numbands
			self.datatype = gdal_object.datatype
			self.set_bandlist(BandList=None)
			self.array=None
			self.meta1={'sensor':None,'path':None,'row':None,'year':None,'doy':None,'product':None,'scene_name':None}
			if IS_LANDSAT:
				self.parse_filename_lndsat()
			#self.Create_Dataset()
	def parse_filename_lndsat(self):
		fname=os.path.basename(os.path.basename(self.path))
		words=os.path.splitext(fname)[0].split('.')
		if len(words) > 1:
			self.meta1['scene_name']=words[1][0:16]
			self.meta1['product']=words[0]
			self.meta1['sensor']=words[1][0:3]
			self.meta1['path']=int(words[1][3:6])
			self.meta1['row']=int(words[1][6:9])
			self.meta1['year']=int(words[1][9:13])
			self.meta1['doy']=int(words[1][13:16])
		else:
			self.meta1['scene_name']=words[0][0:16]
			self.meta1['product']='unknown'
			self.meta1['sensor']=words[0][0:3]
			self.meta1['path']=int(words[0][3:6])
			self.meta1['row']=int(words[0][6:9])
			self.meta1['year']=int(words[0][9:13])
			self.meta1['doy']=int(words[0][13:16])	
	def Create_Dataset(self,format='GTiff'):
		image_driver = gdal.GetDriverByName(format)
		projection = self.projection
		geotransform = Create_GeoTransform(self.extent)
		#self.extent.window['startx']=0		####################      HACK               
		#self.extent.window['starty']=0
		self.dataset=image_driver.Create(self.path,self.extent.xsize,self.extent.ysize,self.numbands,self.datatype)
		self.dataset.SetProjection(projection)
		self.dataset.SetGeoTransform(geotransform)
	def set_bandlist(self,BandList=None):
		if BandList:
			self.bandlist=BandList
			self.numbands=len(self.bandlist)
		else:
			self.bandlist=range(1,self.numbands+1)
			self.numbands=len(self.bandlist)
	def get_multiband_array(self,BandList=None):
		if BandList:
			blist=BandList
		else:
			blist=self.bandlist
		array_list=list()
		for b in blist:
			array_list.append(self.dataset.GetRasterBand(b).ReadAsArray(self.extent.window['startx'],self.extent.window['starty'],self.extent.window['xwin'],self.extent.window['ywin']))
		self.array = np.dstack(array_list)
	def write_array(self):
		#check to see if image object array exist
		if self.array==None:
			#self.get_multiband_array(self.bandlist)
			print "No array exists"
		count=1
		while count <= self.dataset.RasterCount:
			image_band_out=self.dataset.GetRasterBand(count)
			if self.datatype==6:
				image_band_out.SetNoDataValue(FLOAT32_NODATA)
			if self.datatype==3:
				image_band_out.SetNoDataValue(INT16_NODATA)
			image_band_out=gdal_array.BandWriteArray(image_band_out,self.array[:,:,count-1])
			count=count+1
		os.chmod(self.path,ACCESS_MODE)
######################################## Functions   ###########################################
def Get_GeoTransform(gdal_object):
	if isinstance(gdal_object,gdal.Dataset):
		geotransform=gdal_object.GetGeoTransform()
		return geotransform
	if isinstance(gdal_object,gdal_file):
		geotransform=gdal_object.geotransform
		return geotransform
	if isinstance(gdal_object,str):
		if os.path.isfile(gdal_object):
			ds=gdal.Open(gdal_object)
			if ds != None:
				geotransform=ds.GetGeoTransform()
				return geotransform
			else:
				print gdal_object," is not a GDAL supported dataset"
				return None
		else:
			print gdal_object," is not a valid path in the file system"
			return None

def Create_GeoTransform(extent):
	geotransform = extent.geotransform
	UL=(extent.window['startx'],extent.window['starty'])
	xy_proj=GetProjCoord(UL,geotransform)
	geotransform=(xy_proj[0],geotransform[1],geotransform[2],xy_proj[1],geotransform[4],geotransform[5])
	return geotransform
	
def Get_Extent(gdal_object):
	if isinstance(gdal_object,gdal.Dataset):
		geotransform=gdal_object.GetGeoTransform()
		rasterxsize=gdal_object.RasterXSize
		rasterysize=gdal_object.RasterYSize
		extent = Get_Corners(geotransform,rasterxsize,rasterysize)
		return extent
	if isinstance(gdal_object,gdal_file):
		geotransform=gdal_object.geotransform
		rasterxsize=gdal_object.rasterXSize
		rasterysize=gdal_object.rasterYSize
		extent = Get_Corners(geotransform,rasterxsize,rasterysize)
		return extent
	if isinstance(gdal_object,str):
		if os.path.isfile(gdal_object):
			ds=gdal.Open(gdal_object)
			if ds != None:
				geotransform=ds.GetGeoTransform()
				rasterxsize=ds.RasterXSize
				rasterysize=ds.RasterYSize
				extent = Get_Corners(geotransform,rasterxsize,rasterysize)
				return extent
			else:
				print gdal_object," is not a GDAL supported dataset"
				return None
		else:
			print gdal_object," is not a valid path in the file system"
			return None	

def Get_Corners(gt,xsize,ysize):
	edges={'upper':None,'left':None,'lower':None,'right':None}
	edges['upper']=gt[3]
	edges['left']=gt[0]
	edges['lower']=(gt[3] + (ysize* gt[5]))
	edges['right']=(gt[0] + (xsize* gt[1]))
	return edges

def GetRasterCoord(xy_tuple,geotransform):		#returns in xy format(col/row)
	gt=geotransform
	xpix=int(round((xy_tuple[0]-gt[0])/abs(gt[1])))
	ypix=int(round((gt[3]-xy_tuple[1])/abs(gt[5])))
	#lets do the back check
	xpix-=1
	ypix-=1
	X=gt[0]+(xpix*gt[1])
	x_range=(X,X+gt[1])
	while not ((xy_tuple[0] >= x_range[0])&(xy_tuple[0] < x_range[1])):
		xpix +=1
		X=gt[0]+(xpix*gt[1])
		x_range=(X,X+gt[1])
		
	Y=gt[3]+(ypix*gt[5])
	y_range=(Y,Y+gt[5])
	while not ((xy_tuple[1] <= y_range[0])&(xy_tuple[1] > y_range[1])):
		ypix+=1
		Y=gt[3]+(ypix*gt[5])
		y_range=(Y,Y+gt[5])
	return (xpix,ypix)
	
def GetProjCoord(xy_tuple,geotransform):	#returns in xy format
	gt=geotransform
	x=(xy_tuple[0]*gt[1])+gt[0]#+(gt[1]/2)
	y=(gt[3]+(xy_tuple[1]*gt[5]))#-(gt[1]/2)
	return (x,y)
