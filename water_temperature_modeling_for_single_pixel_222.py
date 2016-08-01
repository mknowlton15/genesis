'''
Author: Haiyong Ding
School of remote sensing, Nanjing University of Information Science and Technology,
Last edit date: July 13, 2014
'''
import os
from osgeo import gdal
import numpy as np
from numpy import *

import matplotlib.pyplot as plt
import lmfit
from lmfit import minimize, Parameters, report_fit
from pandas.tools.plotting import scatter_matrix
import sys
import osgeo.ogr as ogr
import osgeo.osr as osr
from osgeo import gdal_array
from osgeo.gdalconst import *
import math
from scipy.optimize import curve_fit
from pandas import *
import pandas as pd
import struct
import pylab

def writeFloatRaster(outpath,inArray,template_DS,FLOAT_NODATA=-9999,mode=0777):
	#writes a single band raster to the same extent and projection as the template DS
	image_driver=template_DS.GetDriver()
	image_out_DS=image_driver.Create(outpath,template_DS.RasterXSize,template_DS.RasterYSize,1,GDT_Float32)
	projection=template_DS.GetProjectionRef()
	geotransform=template_DS.GetGeoTransform()
	image_out_DS.SetProjection(projection)
	image_out_DS.SetGeoTransform(geotransform)
	image_band_out=image_out_DS.GetRasterBand(1)
	image_band_out.SetNoDataValue(FLOAT_NODATA)
	image_band_out=gdal_array.BandWriteArray(image_band_out,inArray)
	os.chmod(outpath,mode)
	image_out_DS=None

#os.system("pause")

# get the working directory from cmd line
#wd = sys.argv[1] # get the name of the fold

#################################### using curve_fit model

def Using_WST_curve_fitting(DOY, pixel):
	DOY_normal = (DOY-mean(DOY))/std(DOY)

	WST_return_value  = WST_curve_fitting_model(pixel, DOY, weight=None)
	weighted = np.divide(1, abs(WST_return_value[0]))
	weighted = weighted/weighted.sum()
	WST_return_value  = WST_curve_fitting_model(pixel, DOY, weight=weighted)

	#print "WST fitting results is ", WST_return_value[1:5]
	fitting_para = WST_return_value[1:5]
	print fitting_para

	y_fitted = new_func_2(fitting_para, DOY_normal)

	plot_multi_scale_data(DOY, pixel, y_fitted)

#Using_WST_curve_fitting(DOY, pixel)



#############################################################################
def extrema(x, y):
	max_v = y[0]
	min_v = y[0]
	max_index = 0
	min_index = 0
	for iter in arange(len(y)):
		if y[iter]>max_v:
			max_v = y[iter]
			max_index = iter
		if y[iter]<min_v:
			min_v = y[iter]
			min_index = iter
	return max_v, x[max_index], min_v, x[min_index]

def max_min_water_temperature():
#	 wd = 'F:\Chesapeak_Bay_water_temperature_data'
	file_a = 'F:\\Chesapeak_bay_files_NEW\\A1.tif'
	file_b = 'F:\\Chesapeak_bay_files_NEW\\B1.tif'
	file_c = 'F:\\Chesapeak_bay_files_NEW\\C1.tif'
	#file_d = 'F:\D1_v3.tif'

	g1 = gdal.Open(file_a)
	g2 = gdal.Open(file_b)
	g3 = gdal.Open(file_c)
	#g4 = gdal.Open(file_d)

	row_number = g1.RasterYSize
	col_number = g1.RasterXSize

	max_temperature = np.zeros((row_number, col_number), dtype=float)
	date_max_temp = np.zeros((row_number, col_number), dtype=float)
	min_temperature = np.zeros((row_number, col_number), dtype=float)
	date_min_temp = np.zeros((row_number, col_number), dtype=float)
	row = 0
	col = 0
	while row<row_number:
		col = 0
		print row
		while col < col_number:
			#print row, col
			a = g1.ReadAsArray(col, row, 1, 1)
			b = g2.ReadAsArray(col, row, 1, 1)
			c = g3.ReadAsArray(col, row, 1, 1)
			#d = g4.ReadAsArray(col, row, 1, 1)
			#print len(a)
			if a==0 or b==0 or c==0:
				max_temperature[row, col] = 0
				date_max_temp[row, col] = 0
				min_temperature[row, col] = 0
				date_min_temp[row, col] = 0
				pass
			else:
				doy = linspace(1, 365, 365)
				#doy_normal = (doy-mean(doy))/std(doy)
				y = a*np.sin(2*np.pi*doy/365 + b) + c
				try:
					p = extrema(doy, y[0])
					#print p
				except RuntimeError:
					p = [0, 0, 0, 0]
					continue
				max_temperature[row, col] = p[0]
				date_max_temp[row, col] = p[3]
				min_temperature[row, col] = p[2]
				date_min_temp[row, col] = p[1]
			col +=1
		row +=1

	out_path_A = 'F:\max_temperature.tif'
	writeFloatRaster(out_path_A, max_temperature, g1)
	out_path_B = 'F:\date_max_temp.tif'
	writeFloatRaster(out_path_B, date_max_temp, g1)
	out_path_C = 'F:\min_temperature.tif'
	writeFloatRaster(out_path_C, min_temperature, g1)
	out_path_D = 'F:\date_min_temp.tif'
	writeFloatRaster(out_path_D, date_min_temp, g1)
# run max_min_water_temperature function
#max_min_water_temperature()

def func(x, p):
	#return p[0]*sin(p[1]*x+p[2])+p[3]
	return p[0]*sin(2*np.pi*x/365 + p[1]) + p[2]

def residuals(p, y, x):
	return y-func(x, p)

def new_func(p, x):
	a = p['a'].value
	b = p['b'].value
	c = p['c'].value
	return a*np.sin(2*np.pi*x/365 + b) + c

def new_func_2(par, x):
	#return par[0]*np.sin(par[1]*x + par[2]) + par[3]
	return par[0]*np.sin(2*np.pi*x/365 + par[1]) + par[2]

def new_residual(p, y, x, weight):
    return weight*(y-new_func(p, x))
    
	#if weight == 1:
	#	return y - new_func(p, x)
	#else:
	#	return weight*(y-new_func(p, x))

class water_surface_temperature:
	def __init__(self, file_directory, begin_year=1984,end_year=1990):
		self.directory = file_directory
		self.begin_year = begin_year
		self.end_year = end_year
		file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected\\Land_Water_QA_2011282.tif"
		self.land_water_file = file
		g = gdal.Open(file)
		self.row_number = g.RasterYSize
		self.col_number = g.RasterXSize
	def get_dataset(self):
		wd = self.directory
		# get file list
		dir_list=os.listdir(wd)
		# parse list and extract .tif file only to image_list
		WST = list()		#create a empty list
		cloud = list()

		for file in dir_list:
			p=os.path.splitext(file)
			if p[0][3:5] == 'B6':
				year = int(p[0][6:10])
				if year>= self.begin_year and year <= self.end_year:
					doy = int(p[0][10:])
					WST.append((year, doy, os.path.join(wd,file)))
			if p[0][0:2] == 'CQ':
				year = int(p[0][3:7])
				if year>=self.begin_year and year <= self.end_year:
					#year = int(p[0][3:7])
					doy = int(p[0][7:])
					cloud.append((year, doy, os.path.join(wd,file)))

		#for item in WST:
		#	 print item
		#for item in cloud:
		#	 print item
			# sort the data according to the DOY order
		WST.sort(lambda x,y: cmp((x[1]),(y[1])))
		cloud.sort(lambda x,y: cmp((x[1]),(y[1])))
		
		WST_dataset = list()
		cloud_dataset = list()
		DOY = list()
		
		for item in arange(0, len(WST)):
			p1 = WST[item]
			WST_dataset.append((p1[0],p1[1],gdal.Open(p1[2])))
			DOY.append(p1[1])
		for item in cloud:
			cloud_dataset.append((item[0], item[1], gdal.Open(item[2])))

		self.WST = WST_dataset
		self.cloud = cloud_dataset
		self.data_length = len(WST)
		self.DOY = np.asarray(DOY)
		
	def get_one_pixel(self, row=3635, col=3215):
		#print Data_set.size, Data_set.shape, len(b6_list)
		g = gdal.Open(self.land_water_file)
		land_water_pixel = g.ReadAsArray(col, row, 1, 1)
		print land_water_pixel[0]
		if land_water_pixel[0]==0:
			print 'this pixel is land'
			pass
		else:
			#Data_set = np.zeros((len(self.WST), 4), dtype=int16)
			pixel = []
			doy = []
			for iter in arange(0, len(self.WST)):
				#print iter
				g2 = self.cloud[iter][2]
				cloud_pixel = g2.ReadAsArray(col, row, 1, 1)
				if cloud_pixel[0] == 0: # no cloud
					g1 = self.WST[iter][2]
					#Data_set[iter,3] = cloud_pixel[0]
					water_pixel = g1.ReadAsArray(col, row, 1, 1)
					if water_pixel[0] > -10.0:
						pixel.append(water_pixel[0][0])
						p1 = self.WST[iter]
						doy.append(p1[1])
					else:
						pass
				else:
					pass

		DOY = np.asarray(doy)
		pixel = np.asarray(pixel)
		plt.figure()
		plt.plot(DOY, pixel, 'k*')
		plt.show()
		self.DOY = DOY
		self.pixel = pixel

	def WST_curve_fitting_model(self):
		p1 = lmfit.Parameters()
		p1.add('a', value=10.5)
		#p1.add('b', value=2*np.pi/365, vary=False)
		p1.add('b', value=8.0)
		p1.add('c', value=11.5)

		DOY = self.DOY
		pixel = self.pixel
		
		#print 'pixel is ', pixel, DOY
		mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, 1))
		#print "pixel is "
		weighted = np.divide(1, abs(mi.residual))
		weighted = weighted/weighted.sum()
		mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, weighted))
		
		report_fit(mi.params)
		ci = lmfit.conf_interval(mi)
		lmfit.printfuncs.report_ci(ci)
		
		a = mi.params['a'].value
		b = mi.params['b'].value
		c = mi.params['c'].value
		#d = mi.params['d'].value
			   
		rmse = mi.residual**2
		rmse = np.sqrt(rmse.sum()/len(rmse))
		WST_return_value = (mi.residual, a, b, c, rmse)
		self.fit_result = WST_return_value

	def plot_fitting_curve(self):
		DOY = self.DOY
		fitted_value = func(DOY, self.fit_result[1:4])
		plt.figure()
		plt.plot(self.DOY, self.pixel, 'k*')
		plt.plot(self.DOY, fitted_value)
		plt.show()
		
	def bootstraping_estimate(self):
		residual = self.fit_result[0]

		noise_numb = 1000 # noise number
		bootstrap_para = np.zeros((noise_numb, 3), dtype=float)

		for iter in range(noise_numb):
			residual_1 = np.random.choice(residual, len(residual))
			new_pixel = self.pixel + residual_1

			p1 = lmfit.Parameters()
			p1.add('a', value=10.5, vary=True)
			#p1.add('b', value=2*np.pi/365, vary=False)
			p1.add('b', value=8.0, vary=True)
			p1.add('c', value=11.5, vary=True)

			mi = lmfit.minimize(new_residual, p1, args=(new_pixel, self.DOY, 1))
			#weighted = np.divide(1, abs(mi.residual))
			#weighted = weighted/weighted.sum()
			#mi = lmfit.minimize(new_residual, p1, args=(new_pixel, self.DOY, weighted))

			a = mi.params['a'].value
			b = mi.params['b'].value
			c = mi.params['c'].value
			#d = mi.params['d'].value

			bootstrap_para[iter, 0] = np.abs(a)
			bootstrap_para[iter, 1] = np.abs(b)
			bootstrap_para[iter, 2] = np.abs(c)
			#bootstrap_para[iter, 3] = WST_return_value[4]

		df = DataFrame(bootstrap_para, columns=['a', 'b', 'c'])
		#scatter_matrix(df, alpha=0.2, figsize=(6,6), diagonal='kde')
		scatter_matrix(df, alpha=0.2, figsize=(8, 8), diagonal='hist')
		plt.show()
		print mean(bootstrap_para[:, 0]), std(bootstrap_para[:, 0])
		print mean(bootstrap_para[:, 1]), std(bootstrap_para[:, 1])
		print mean(bootstrap_para[:, 2]), std(bootstrap_para[:, 2])
		
	def get_parameters_for_subimage(self, col_0, row_0, col_size, row_size):
		# the subset is from corner (col=4910, row=4500),have 200*200 size

		g = gdal.Open(self.land_water_file)
		#row_0 = 4500
		#col_0 = 4910
		row_number = row_0+row_size
		col_number = col_0+col_size
		amplitude = np.zeros((row_size, col_size), dtype=float)
		#B = np.zeros((row_number, col_number), dtype=float)
		phase = np.zeros((row_size, col_size), dtype=float)
		mean_t = np.zeros((row_size, col_size), dtype=float)
		rmse_im = np.zeros((row_size, col_size), dtype=float)
		length = self.data_length

		for col in arange(col_0, col_number):
			#print "col is ", col
			for row in arange(row_0, row_number):
				#print "row is ", row
				land_water_pixel = g.ReadAsArray(col, row, 1, 1)
				#print land_water_pixel[0]
				if land_water_pixel[0]==0:
					#print 'This pixel is land'
					amplitude[col-col_0, row-row_0] = 0
					phase[col-col_0, row-row_0] = 0
					mean_t[col-col_0, row-row_0] = 0
					rmse_im[col-col_0, row-row_0] = 0
					pass
				else:
					pixel = []
					DOY = []
									 
					for iter in arange(0, length):
						#print iter
						g2 = self.cloud[iter][2]
						cloud_pixel = g2.ReadAsArray(col, row, 1, 1)
						if cloud_pixel[0][0] == 255:
							pass
						else:
							g1 = self.WST[iter][2]
							wst_pixel = g1.ReadAsArray(col, row, 1, 1)
							if wst_pixel[0][0]<0:
								pass
							else:
								pixel.append(wst_pixel[0][0])
								
								p1 = self.WST[iter]
								DOY.append(p1[1])
				
					pixel = np.asarray(pixel)
					DOY = np.asarray(DOY)
					if len(pixel) < 15:
						a = 0
						b = 0
						c = 0
						rmse = 0
						pass
					else:
						p1 = lmfit.Parameters()
						p1.add('a', value=-10.5, vary=True)
						#p1.add('b', value=2*np.pi/365, vary=False)
						p1.add('b', value=8.0, vary=True)
						p1.add('c', value=11.5, vary=True)

						mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, 1))
						weighted = np.divide(1, abs(mi.residual))
						weighted = weighted/weighted.sum()
						mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, weighted))

						a = mi.params['a'].value
						b = mi.params['b'].value
						c = mi.params['c'].value 

						rmse = mi.residual**2
						rmse = np.sqrt(rmse.sum()/len(rmse))

					amplitude[col-col_0, row-row_0] = np.abs(a)
					phase[col-col_0, row-row_0] = b
					mean_t[col-col_0, row-row_0] = c
					rmse_im[col-col_0, row-row_0] = rmse
		#print "mean of amplitude is ", mean(amplitude)
		#print "mean of phase is ", mean(phase)
		#print "mean of mean_temperature is ", mean(mean_t)
		#print "mean of rmse is ", mean(rmse_im)
		self.sub_max_mean_min = (np.max(amplitude+mean_t), np.mean(mean_t), np.min(mean_t-amplitude))
		
	def get_slope_for_subimage(self, col_0, row_0, col_size, row_size):
		# the subset is from corner (col=4910, row=4500),have 200*200 size
		WST = self.WST
		cloud = self.cloud
		
		#g = gdal.Open(self.land_water_file)
		#row_0 = 4500
		#col_0 = 4910
		row_number = row_0+row_size
		col_number = col_0+col_size
		DOY = self.DOY

		length = self.data_length
		dataset = list()
		for doy in arange(0, length):
			
			#data = list()
			g = cloud[doy][2] # cloud 255, 
			cloud_pixels = g.ReadAsArray(col_0, row_0, col_size, row_size)
			cloud_pixels = cloud_pixels.flatten()
			#print cloud_pixels
			g1 = WST[doy][2]
			water = g1.ReadAsArray(col_0, row_0, col_size, row_size)
			water = water.flatten()
			#print 'water is ', water
			water = water[cloud_pixels==0]
			#print 'water without cloud is ', water
			water = water[water>0]
			#print water
			mean_value = mean(np.asarray(water))
			#print 'mean_value is ', mean_value
			
			if np.isnan(mean_value):
				pass
			else:
				dataset.append((DOY[doy], mean_value))
		#print dataset
		dataset = np.asarray(dataset)
		if dataset.shape[0] < 5:
			self.return_mean_T_for_slope = 0.0
			#print "return mean Temperature is 0"
		else:
			doy_0 = dataset[:, 0]
			pixel = dataset[:, 1]
		
			p1 = lmfit.Parameters()
			p1.add('a', value=-10.5, vary=True)
			#p1.add('b', value=2*np.pi/365, vary=False)
			p1.add('b', value=8.0, vary=True)
			p1.add('c', value=11.5, vary=True)

			mi = lmfit.minimize(new_residual, p1, args=(pixel, doy_0, 1))
			weighted = np.divide(1, abs(mi.residual))
			weighted = weighted/weighted.sum()
			mi = lmfit.minimize(new_residual, p1, args=(pixel, doy_0, weighted))

			# a = mi.params['a'].value
			# b = mi.params['b'].value
			c = mi.params['c'].value

							# rmse = mi.residual**2
							# rmse = np.sqrt(rmse.sum()/len(rmse))

			#print "mean of amplitude is ", mean(amplitude)
			#print "mean of phase is ", mean(phase)
			#print "mean of mean_temperature is ", c
			#print "mean of rmse is ", mean(rmse_im)
			self.return_mean_T_for_slope = c
			
	def get_slope_for_each_pixel(self, col_0, row_0, col_size, row_size):
		# the subset is from corner (col=4910, row=4500),have 200*200 size
		WST = self.WST
		cloud = self.cloud

		DOY = self.DOY

		length = self.data_length
		dataset = list()
		
		pixel = []
		doy = []
		for iter in arange(0, len(self.WST)):
			#print iter
			g2 = cloud[iter][2]
			cloud_pixel = g2.ReadAsArray(col_0, row_0, 1, 1)
			if cloud_pixel[0] == 0: # no cloud
				g1 = WST[iter][2]
				#Data_set[iter,3] = cloud_pixel[0]
				water_pixel = g1.ReadAsArray(col_0, row_0, 1, 1)
				if water_pixel[0] > -10.0:
					pixel.append(water_pixel[0][0])
					p1 = WST[iter]
					doy.append(p1[1])
				else:
					pass
			else:
				pass

		DOY = np.asarray(doy)
		pixel = np.asarray(pixel)
		#plt.figure()
		#plt.plot(DOY, pixel, 'k*')
		#plt.show()

		if len(pixel) < 5:
			self.return_mean_T_for_pixel = 0.0
			#print "return mean Temperature is 0"
		else:
			doy_0 = DOY
	
			p1 = lmfit.Parameters()
			p1.add('a', value=-10.5, vary=True)
			#p1.add('b', value=2*np.pi/365, vary=False)
			p1.add('b', value=8.0, vary=True)
			p1.add('c', value=11.5, vary=True)

			mi = lmfit.minimize(new_residual, p1, args=(pixel, doy_0, 1))
			weighted = np.divide(1, abs(mi.residual))
			weighted = weighted/weighted.sum()
			mi = lmfit.minimize(new_residual, p1, args=(pixel, doy_0, weighted))

			# a = mi.params['a'].value
			# b = mi.params['b'].value
			c = mi.params['c'].value

							# rmse = mi.residual**2
							# rmse = np.sqrt(rmse.sum()/len(rmse))

			#print "mean of amplitude is ", mean(amplitude)
			#print "mean of phase is ", mean(phase)
			#print "mean of mean_temperature is ", c
			#print "mean of rmse is ", mean(rmse_im)
			self.return_mean_T_for_pixel = c	

	def temperature_lagging(self): 
		land_water_file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected\\Land_Water_QA_2011282.tif"
		g = gdal.Open(land_water_file)
		sample_location = [[2502, 4706],[2498, 4709], [2497,4717], [2508, 4736], [2485, 4792]]
		parameter = np.zeros((len(sample_location), 3), dtype=float)

		for it in arange(len(sample_location)):
			row = sample_location[it][0]
			col = sample_location[it][1]
			land_water_pixel = g.ReadAsArray(col, row, 1, 1)
			if land_water_pixel[0] == 0:
				print 'this pixel is land'
				pass
			else:
				pixel = []
				DOY = []
				for iter in arange(0, len(self.WST)):
				
					g2 = self.cloud[iter][2]
					cloud_pixel = g2.ReadAsArray(col, row, 1, 1)
					if cloud_pixel[0] == 0: # no cloud
						g1 = self.WST[iter][2]
						water_pixel = g1.ReadAsArray(col, row, 1, 1)
						#print water_pixel[0][0]
						if water_pixel[0][0] > 0:
							pixel.append(water_pixel[0][0])
							day = self.WST[iter]
							DOY.append(day[1])
						else:
							pass
					else:
						pass
				#print pixel
				#print 'DOY is ', doy
				pixel = np.asarray(pixel)
				DOY = np.asarray(DOY)
				
				p1 = lmfit.Parameters()
				p1.add('a', value=-10.5, vary=True)
				p1.add('b', value=8.0, vary=True)
				p1.add('c', value=11.5, vary=True)

				mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, 1))
				weighted = np.divide(1, abs(mi.residual))
				weighted = weighted/weighted.sum()
				mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, weighted))
			
				a = mi.params['a'].value
				b = mi.params['b'].value
				c = mi.params['c'].value
				   
				rmse = mi.residual**2
				rmse = np.sqrt(rmse.sum()/len(rmse))
				
				parameter[it, 0] = a
				parameter[it, 1] = b
				parameter[it, 2] = c
				
		print 'the fitting parameter is ', parameter
		lable1 = u'a: '+ str(format(np.abs(parameter[0,0]), '.2f'))+', b: '+str(format(parameter[0,1], '.2f'))+\
		', c: '+ str(format(parameter[0,2], '.2f'))
		plt.plot(DOY, func(DOY, (parameter[0,0],parameter[0,1],parameter[0,2])),'k', label=lable1)

		lable2 = u'a: '+ str(format(np.abs(parameter[1,0]), '.2f'))+', b: '+str(format(parameter[1,1], '.2f'))+\
		', c: '+ str(format(parameter[1,2], '.2f'))
		plt.plot(DOY, func(DOY, (parameter[1,0],parameter[1,1],parameter[1,2])),'k--', label=lable2)

		lable3 = u'a: '+ str(format(np.abs(parameter[2,0]), '.2f'))+', b: '+str(format(parameter[2,1], '.2f'))+\
		', c: '+ str(format(parameter[2,2], '.2f'))
		plt.plot(DOY, func(DOY, (parameter[2,0],parameter[2,1],parameter[2,2])),'k-.', label=lable3)

		lable4 = u'a: '+ str(format(np.abs(parameter[3,0]), '.2f'))+', b: '+str(format(parameter[3,1], '.2f'))+\
		', c: '+ str(format(parameter[3,2], '.2f'))
		plt.plot(DOY, func(DOY, (parameter[3,0],parameter[3,1],parameter[3,2])),'k:', label=lable4)


		#lable5 = u'a: '+ str(format(np.abs(parameter[4,0]), '.2f'))+', b: '+str(format(parameter[4,1], '.2f'))+\
		#', c: '+ str(format(parameter[4,2], '.2f'))
		#plt.plot(DOY, func(DOY, (parameter[4,0],parameter[4,1],parameter[4,2])),'k-', label=lable5)

		plt.xlabel('Day of Year')
		plt.ylabel('Water surface temperature')
		#plt.ylim(-5, 30)
		plt.xlim(-5, 370)
		plt.legend(loc='best') 
		plt.show()
	
################## plot for single pixel
#file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected"
#wst_all = water_surface_temperature(file, 1984, 2011)
#wst_all.get_dataset()	
#wst_all.get_one_pixel()
#wst_all.WST_curve_fitting_model()
#wst_all.plot_fitting_curve()
#wst_all.bootstraping_estimate()
#wst_all.temperature_lagging()


####################### regression for each block

import pandas as pd
import scipy
#import statsmodels.api as sm
import statsmodels.api as sm
from xlrd import open_workbook
import xlwt

'''
file = 'F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected'
wst_all = water_surface_temperature(file, 1984, 2011)
wst_all.get_dataset()
#g = gdal.Open(wst_all.land_water_file)


row_number = wst_all.row_number
col_number = wst_all.col_number

length = wst_all.data_length

slope = np.zeros((row_number, col_number), dtype=float)
p_value = np.zeros((row_number, col_number), dtype=float)
xBSize = 40 #7
yBSize = 40 # 7
for row in range(0, row_number, yBSize):
	print 'row = ', row
	if row+yBSize < row_number:
		numRows = yBSize
	else:
		numRows = row_number - row - 1
  
	for col in range(0, col_number, xBSize): 
		print 'col = ', col
		if col+xBSize < col_number:
			numCols = xBSize
		else:
			numCols = col_number - row - 1
		mask_file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected\\Land_Water_QA_2011282.tif"	   
		g = gdal.Open(mask_file)
		land_water_pixel = g.ReadAsArray(col, row, numCols, numRows)
		#print sum(land_water_pixel>0)
		if sum(land_water_pixel>0)<800: # this means in the 7*7 region, there are more than 500 pixels are land pixels
			slope[row:row+numRows, col:col+numCols] = 0.0
			p_value[row:row+numRows, col:col+numCols] = 0.0
			pass
		else:
			year_Tmean = list()
			for year in arange(1984, 2008):
				wst = water_surface_temperature(file, year, year+5)
				wst.get_dataset()
				wst.get_slope_for_subimage(col, row, numCols, numRows)
				#print 'return_mean_T_for_slope is'
				#print wst.return_mean_T_for_slope
				year_Tmean.append((year, wst.return_mean_T_for_slope))
				   # # convert list to dataframe
			df = pd.DataFrame(year_Tmean,  columns= ['Year', 'S_Tmean'])
			   # #print df
			T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
			res = T_water_year.fit()
			slope_0 = res.params[1]
			print 'slope_0 is ', slope_0
			slope[row:row+numRows, col:col+numCols] = slope_0
			p_value[row:row+numRows, col:col+numCols] = res.f_pvalue

out_path_slope = "F:\\slope_9.tif"
writeFloatRaster(out_path_slope, slope, wst_all.WST[0][2])
out_path_p_value = "F:\\p_value_9.tif"
writeFloatRaster(out_path_p_value, p_value, wst_all.WST[0][2])
'''


##########################################################################################
############### rate of temperature on every pixel ###################################
#######################################################################################


file = 'F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected'
wst_all = water_surface_temperature(file, 1984, 2011)
wst_all.get_dataset()
#g = gdal.Open(wst_all.land_water_file)


row_number = wst_all.row_number
col_number = wst_all.col_number

#length = wst_all.data_length

slope = np.zeros((row_number, col_number), dtype=float)
p_value = np.zeros((row_number, col_number), dtype=float)

for row in range(0, row_number):
	print 'row = ', row
  
	for col in range(0, col_number): 
		mask_file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected\\Land_Water_QA_2011282.tif"	   
		g = gdal.Open(mask_file)
		land_water_pixel = g.ReadAsArray(col, row, 1, 1)
		#print land_water_pixel[0][0]
		if land_water_pixel[0][0]==0: 
			slope[row, col] = 0.0
			p_value[row, col] = 0.0
			pass
		else:
			year_Tmean = list()
			for year in arange(1984, 2008):
				wst = water_surface_temperature(file, year, year+5)
				wst.get_dataset()
				wst.get_slope_for_each_pixel(col, row, 1, 1)
				#print 'return_mean_T_for_slope is'
				#print wst.return_mean_T_for_slope
				year_Tmean.append((year, wst.return_mean_T_for_pixel))
				   # # convert list to dataframe
			df = pd.DataFrame(year_Tmean,  columns= ['Year', 'S_Tmean'])
			   # #print df
			T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
			res = T_water_year.fit()
			slope_0 = res.params[1]
			#print 'slope_0 is ', slope_0
			slope[row, col] = slope_0
			p_value[row, col] = res.f_pvalue

out_path_slope = "F:\\slope_for_pixel.tif"
writeFloatRaster(out_path_slope, slope, wst_all.WST[0][2])
out_path_p_value = "F:\\p_value_for_pixel.tif"
writeFloatRaster(out_path_p_value, p_value, wst_all.WST[0][2])


############################################################################################
############################################################################################

################# testing function : get_slope_for_subimage
#file = 'F:\Chesapeak_Bay_water_temperature_data'
#wst_all = water_surface_temperature(file, 1984, 1988)
#wst_all.get_dataset()	 
#wst_all.get_slope_for_subimage(5460, 4320, 500, 100) 
	
'''
####################### get WST for the entire region #################
file = "F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected"
wst_all = water_surface_temperature(file, 1984, 2011)
wst_all.get_dataset()

row_number = wst_all.row_number
col_number = wst_all.col_number
amplitude = np.zeros((row_number, col_number), dtype=float)
phase = np.zeros((row_number, col_number), dtype=float)
mean_t = np.zeros((row_number, col_number), dtype=float)

length = wst_all.data_length

g = gdal.Open(wst_all.land_water_file)
row = 0
col = 0
while row < row_number:
	col = 0
	print row
	while col < col_number:
		#print col
		land_water_pixel = g.ReadAsArray(col, row, 1, 1)
		#print land_water_pixel[0]
		if land_water_pixel[0][0]==0:
			#print "This pixel is land"
			amplitude[row, col] = 0
			phase[row, col] = 0
			mean_t[row, col] = 0

			pass
		else:
			pixel = []
			DOY = []
			for iter in arange(0, length):
				#print iter
			
				g2 = wst_all.cloud[iter][2]
				cloud_pixel = g2.ReadAsArray(col, row, 1, 1)
				if cloud_pixel[0][0] == 255:
					pass
				else:
					g1 = wst_all.WST[iter][2]
					wst_pixel = g1.ReadAsArray(col, row, 1, 1)
					if wst_pixel[0][0]<0:
						pass
					else:
						pixel.append(wst_pixel[0][0])
							
						p1 = wst_all.WST[iter]
						DOY.append(p1[1])
			pixel = np.asarray(pixel)
			DOY = np.asarray(DOY)
			
			if len(pixel) < 15:
				a = 0
				b = 0
				c = 0

				pass
			else:
				p1 = lmfit.Parameters()
				p1.add('a', value=-10.5, vary=True)

				p1.add('b', value=8.0, vary=True)
				p1.add('c', value=11.5, vary=True)

				mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, 1))
				weighted = np.divide(1, abs(mi.residual))
				weighted = weighted/weighted.sum()
				mi = lmfit.minimize(new_residual, p1, args=(pixel, DOY, weighted))

				a = mi.params['a'].value
				b = mi.params['b'].value
				c = mi.params['c'].value
			#print "a = ", np.abs(a)
			#print "b = ", b
			#print "c = ", c
			amplitude[row, col] = np.abs(a)
			phase[row, col] = b
			mean_t[row, col] = c

		col +=1
	row +=1


out_path_A = "F:\\Chesapeak_bay_files_NEW\\AA1.tif"
writeFloatRaster(out_path_A, amplitude, wst_all.WST[0][2])
out_path_B = "F:\\Chesapeak_bay_files_NEW\\BB1.tif"
writeFloatRaster(out_path_B, phase, wst_all.WST[0][2])
out_path_C = "F:\\Chesapeak_bay_files_NEW\\CC1.tif"
writeFloatRaster(out_path_C, mean_t, wst_all.WST[0][2])
'''

################### linear regression ################################
def basic_linear_regression(x, y):
	length = len(x)
	sum_x = sum(x)
	sum_y = sum(y)
	# sigma_x^2, and sigma_xy,
	sum_x_squared = sum(map(lambda a: a*a, x))
	covariance = sum([x[i]*y[i] for i in range(length)])
	
	# Magic formula
	a = (covariance - (sum_x *sum_y)/length)/(sum_x_squared-((sum_x**2)/length))
	b = (sum_y - a*sum_x)/length
	return a, b
###########################################################################
import pandas as pd
import scipy
#import statsmodels.api as sm
import statsmodels.api as sm
from xlrd import open_workbook
import xlwt


file_list = ["Aberdeen.xls", "Annapolis.xls", "Baltimore.xls", "Cambridge.xls", "Fredericksburg.xls", "Royal.xls", "Princess.xls"]

'''
############ Annapolis site
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Annapolis.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
#res.summary()

intecept = res.params[0]
slope = res.params[1]
print 'intecept is ', intecept

T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')

sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')

print 'slope is ', slope
ax1.set_ylim([12, 16])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]
if	res_T_air_year.params[0]>0:
	ax1.text(1986, 15.3, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1986, 15.3, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1986, 15.15, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1986, 15, 'p-value<0.0001')
else:
	ax1.text(1986, 15, 'p-value = %.4f' % res_T_air_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1986, 14.80, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1986, 14.80, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1986, 14.65, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1986, 14.5, 'p-value<0.0001')
else:
	ax1.text(1986, 14.5, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels)
ax1.set_title("Annapolis")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
ax2.set_ylim([12, 15.0])
if intecept>0:
	ax2.text(13.5, 14.3, 'Water Temp = %.2f*Air Temp+%.2f' %(slope, intecept))
else:
	ax2.text(13.5, 14.3, 'Water Temp = %.2f*Air Temp%.2f' %(slope, intecept))
ax2.text(13.5, 14.15, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(13.5, 14.0, 'p-value<0.0001')
else:
	ax2.text(13.5, 14, 'p-value = %.4f' % res.f_pvalue)
	
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels) 
ax2.set_title("Annapolis")

print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()

plt.show()
'''



'''
### Princess site
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Princess.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

#linear regression between air temperature and year
T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([11, 13.5])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]
if	res_T_air_year.params[0]>0:
	ax1.text(1982, 14.5, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1982, 14.5, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1982, 14.4, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1982, 14.3, 'p-value<0.0001')
else:
	ax1.text(1982, 14.3, 'p-value =' % res_T_air_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1982, 14.15, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1982, 14.15, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1982, 14.05, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1982, 13.95, 'p-value<0.0001')
else:
	ax1.text(1982, 13.95, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels, 'best')
ax1.set_title("Princess")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([11.2, 13])
if intecept > 0:
	ax2.text(12.25, 14.2, 'Water Temp = %.2f*Air Temp + %.2f' %(slope, intecept))
else:
	ax2.text(12.25, 14.2, 'Water Temp = %.2f*Air Temp %.2f' %(slope, intecept))

ax2.text(12.25, 14.1, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(12.25, 14, 'p-value<0.0001')
else:
	ax2.text(12.25, 14, 'p-value =%.4f' %res.f_pvalue)
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels, 'best') 
ax2.set_title("Princess")

print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''



'''
# Royal site
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Royal.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

#linear regression between air temperature and year
T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([11, 15])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]

if	res_T_air_year.params[0]>0:
	ax1.text(1996, 13.3, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1996, 13.3, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1996, 13.2, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1996, 13.1, 'p-value<0.0001')
else:
	ax1.text(1996, 13.1, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1996, 13., 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1996, 13., 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1996, 12.9, '$R^2$= %.2f' %res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1996, 12.8, 'p-value<0.0001')
else:
	ax1.text(1996, 12.8, 'p-value = %.4f' % res_T_water_year.f_pvalue)
ax1.legend(labels, 'best')
ax1.set_title("Royal")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([11.2, 12.5])
if intecept > 0:
	ax2.text(13.5, 14.5, 'Water Temp = %.2f*Air Temp + %.2f' %(slope, intecept))
else:
	ax2.text(13.5, 14.5, 'Water Temp = %.2f*Air Temp %.2f' %(slope, intecept))

ax2.text(13.5, 14.4, '$R^2$= %.2f' % res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(13.5, 14.3, 'p-value<0.0001')
else:
	ax2.text(13.5, 14.3, 'p-value =%.4f' %res.f_pvalue)
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels, 'best') 
ax2.set_title("Royal")


print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''


'''
# Fred
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Fredericksburg.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

#linear regression between air temperature and year
T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([11, 14.5])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]

if	res_T_air_year.params[0]>0:
	ax1.text(1982, 15.2, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1982, 15.2, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1982, 15.1, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1982, 15., 'p-value<0.0001')
else:
	ax1.text(1982, 15., 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1982, 14.9, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1982, 14.9, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1982, 14.8, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1982, 14.7, 'p-value<0.0001')
else:
	ax1.text(1982, 14.7, 'p-value = %.4f' % res_T_water_year.f_pvalue)
   

ax1.legend(labels, 'best')
ax1.set_title("Fredericksburg")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([11.2, 12.5])
if intecept > 0:
	ax2.text(12.9, 14.8, 'Water Temp = %.2f*Air Temp + %.2f' %(slope, intecept))
else:
	ax2.text(12.9, 14.8, 'Water Temp = %.2f*Air Temp %.2f' %(slope, intecept))
ax2.text(12.9, 14.7, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(12.9, 14.6, 'p-value<0.0001')
else:
	ax2.text(12.9, 14.6, 'p-value =%.4f' %res.f_pvalue)
	
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels, 'best') 
ax2.set_title("Fredericksburg")
print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''



'''
# Cambridge site
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Cambridge.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

#linear regression between air temperature and year
T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([10.5, 15])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]

if	res_T_air_year.params[0]>0:
	ax1.text(1995, 13.3, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1995, 13.3, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1995, 13.2, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1995, 13.1, 'p-value<0.0001')
else:
	ax1.text(1995, 13.1, 'p-value = %.4f' % res_T_air_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1995, 13., 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1995, 13., 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1995, 12.9, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1995, 12.8, 'p-value<0.0001')
else:
	ax1.text(1995, 12.8, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels)
ax1.set_title("Cambridge")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([11.2, 12.5])
if intecept >0:
	ax2.text(13.5, 14.3, 'Water Temp = %.2f*Air Temp+%.2f' %(slope, intecept))
else:
	ax2.text(13.5, 14.3, 'Water Temp = %.2f*Air Temp %.2f' %(slope, intecept))
ax2.text(13.5, 14.15, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(13.5, 14.0, 'p-value<0.0001')
else:
	ax2.text(13.5, 14.0, 'p-value = %.4f' % res.f_pvalue)
	
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels) 
ax2.set_title("Cambridge")
print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''


'''
############ Baltimore site
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Baltimore.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

#linear regression between air temperature and year
T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
ax1.set_ylim([11, 17.5])
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature", "Mean Water surface temperature"]

if	res_T_air_year.params[0]>0:
	ax1.text(1995, 12.8, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1995, 12.8, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1995, 12.6, '$R^2$= %.2f'% res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1995, 12.4, 'p-value<0.0001')
else:
	ax1.text(1995, 12.4, 'p-value = %.4f' % res_T_air_year.f_pvalue)
	
	
if res_T_water_year.params[0]>0:
	ax1.text(1995, 12.10, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1995, 12.10, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1995, 11.9, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1995, 11.7, 'p-value<0.0001')
else:
	ax1.text(1995, 11.7, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels)
ax1.set_title("Baltimore")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([11.4, 12.5])
#ax2.set_xlim([15, 16.2])
if intecept>0:
	ax2.text(15.1, 14.5, 'Water Temp = %.2f*Air Temp+%.2f' %(slope, intecept))
else:
	ax2.text(15.1, 14.5, 'Water Temp = %.2f*Air Temp%.2f' %(slope, intecept))
ax2.text(15.1, 14.35, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(15.1, 14.20, 'p-value<0.0001')
else:
	ax2.text(15.1, 14.20, 'p-value = %.4f' % res.f_pvalue)
	
ax2_labels = ["Temperature", "Regression Line"]
ax2.legend(ax2_labels, loc=2) 
ax2.set_title("Baltimore")
print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''





'''
################### plot of Chestertown site 
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Chestertown.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([9.5, 14])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]
if	res_T_air_year.params[0]>0:
	ax1.text(1996, 12.3, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1996, 12.3, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1996, 12.24, '$R^2$= %.2f' % res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1996, 12.18, 'p-value<0.0001')
else:
	ax1.text(1996, 12.18, 'p-value = %.4f' % res_T_air_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1996, 12.1, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1996, 12.1, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1996, 12.04, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1996, 11.98, 'p-value<0.0001')
else:
	ax1.text(1996, 11.98, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels, 'best')
ax1.set_title("Chestertown")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([9.6, 11.8])

if intecept>0:
	ax2.text(12.3, 12.9, 'Water Temp = %.2f*Air Temp + %.2f' %(slope, intecept))
else:
	ax2.text(12.3, 12.9, 'Water Temp =%.2f*Air Temp %.2f' %(slope, intecept))

ax2.text(12.3, 12.8, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(12.3, 12.7, 'p-value<0.0001')
else:
	ax2.text(12.3, 12.7, 'p-value = %.4f' % res.f_pvalue)
	
ax2_labels = ["Temperature", "Linear Regression"]
ax2.legend(ax2_labels, 'best') 
ax2.set_title("Chestertown")

print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()
plt.show()
'''


'''
################### plot of Millington site 
fold_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature"
file = "Millington.xls"
excel_file = os.path.join(fold_name, file)

xl = pd.ExcelFile(excel_file)
xl.sheet_names
df = xl.parse("Parameters")
df.head()
mod = sm.OLS.from_formula('S_Tmean ~ A_Tmean', data=df)	 # y and x are column names in the DataFrame
res = mod.fit()
res.summary()

intecept = res.params[0]
slope = res.params[1]

T_air_year = sm.OLS.from_formula('A_Tmean ~ Year', data=df)
res_T_air_year = T_air_year.fit()
T_water_year = sm.OLS.from_formula('S_Tmean ~ Year', data=df)
res_T_water_year = T_water_year.fit()

ax1 = plt.subplot(121)
ax1.plot(df.Year, df.A_Tmean, 'k^', df.Year, df.S_Tmean, 'ko')
sm.graphics.abline_plot(model_results=res_T_air_year, ax=ax1, color='black')
sm.graphics.abline_plot(model_results=res_T_water_year, ax=ax1, color='black')
#ax1.set_ylim([10.5, 13.5])
ax1.set_xlabel("Year")
ax1.set_ylabel("Temperature ($^\circ$C)")
labels=["Mean Air Temperature","Mean Water surface temperature"]
if	res_T_air_year.params[0]>0:
	ax1.text(1982, 13.4, 'Air Temp=%.2f*Year+%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
else:
	ax1.text(1982, 13.4, 'Air Temp=%.2f*Year%.2f' %(res_T_air_year.params[1], res_T_air_year.params[0]))
ax1.text(1982, 13.34, '$R^2$= %.2f' % res_T_air_year.rsquared)
if res_T_air_year.f_pvalue<0.0001:
	ax1.text(1982, 13.28, 'p-value<0.0001')
else:
	ax1.text(1982, 13.28, 'p-value = %.4f' % res_T_air_year.f_pvalue)
	
if res_T_water_year.params[0]>0:
	ax1.text(1982, 13.22, 'Water Temp=%.2f*Year+%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
else:
	ax1.text(1982, 13.22, 'Water Temp=%.2f*Year%.2f' %(res_T_water_year.params[1], res_T_water_year.params[0]))
ax1.text(1982, 13.16, '$R^2$= %.2f'% res_T_water_year.rsquared)
if res_T_water_year.f_pvalue<0.0001:
	ax1.text(1982, 13.1, 'p-value<0.0001')
else:
	ax1.text(1982, 13.1, 'p-value = %.4f' % res_T_water_year.f_pvalue)
	
ax1.legend(labels, 'best')
ax1.set_title("Millington")

ax2 = plt.subplot(122)
ax2.plot(df.A_Tmean, df.S_Tmean, 'k*')
sm.graphics.abline_plot(model_results=res, ax=ax2, color='black')
ax2.set_xlabel("Mean Air Temperature ($^\circ$C)")
ax2.set_ylabel("Mean Water surface temperature ($^\circ$C)")
#ax2.set_ylim([10.5, 11.8])
if intecept>0:
	ax2.text(12.05, 13.5, 'Water Temp = %.2f*Air Temp + %.2f' %(slope, intecept))
else:
	ax2.text(12.05, 13.5, 'Water Temp = %.2f*Air Temp %.2f' %(slope, intecept))

ax2.text(12.05, 13.4, '$R^2$= %.2f'% res.rsquared)
if res.f_pvalue<0.0001:
	ax2.text(12.05, 13.3, 'p-value<0.0001')
else:
	ax2.text(12.05, 13.3, 'p-value = %.4f' % res.f_pvalue)
	
ax2_labels = ["Temperature", "Linear Regression"]
ax2.legend(ax2_labels, 'best') 
ax2.set_title("Millington")

print 'Information about S_Tmean ~ A_Tmean'
print res.mse_resid
print res.summary()

print 'Information about T_air ~ year'
print res_T_air_year.mse_resid
print res_T_air_year.summary()

print 'Information about T_water ~ year'
print res_T_water_year.mse_resid
print res_T_water_year.summary()

plt.show()
'''



###############################################################################
'''
from xlrd import open_workbook
import xlwt

file = 'F:\\Chesapeak_Bay_water_temperature_data_Atmospheric_corrected'
workbook = xlwt.Workbook()
sheet = workbook.add_sheet("Parameters")
	
for year in arange(1984, 2008):
	wst = water_surface_temperature(file, year, year+5)
	wst.get_dataset()
	
	wst.get_parameters_for_subimage(4855, 3268, 60, 60)
	#print wst.sub_max_mean_min

	
	sheet.write(year-1983, 0, year) # row, column, value
	sheet.write(year-1983, 1, wst.sub_max_mean_min[0])
	sheet.write(year-1983, 2, wst.sub_max_mean_min[1])
	sheet.write(year-1983, 3, wst.sub_max_mean_min[2]) 
	os.chdir("E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature")
	workbook.save("Annapolic_wst.xls")
'''


'''
file = "F:\Chesapeak_Bay_water_temperature_data\B6_1984240.tif"
g = gdal.Open(file)
# cambridge site (5940, 4650), 60*60, done
# royal oak site (5400, 4300), 100*100, done
# solomons site (4890, 5695), 60*60, done
# Fredericksburg site (2474, 5459), 60*60, done
# Baltimore site (4601, 2235), 50*50, done
# Millington site (6160, 1790), 50*50
# Colonial beach site (3400, 5900), 100*100, done
# Indianhead site (2680, 4570), 60*60, done
# Annapolic site (4855, 3268), 60*60, done
# Princess site (5600, 6400), 100*100
# Chestertown site (5680, 2855), 100*100
col_0 = 5680
row_0 = 2855
data = g.ReadAsArray(col_0, row_0, 100, 100)
plt.imshow(data)
plt.show()
'''

######################### get air temperature for different stations
'''
from xlrd import open_workbook
import xlwt
import pandas as pd

file_name = "E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature\\Chestertown_MD_181750_4888.xls"
workbook = open_workbook(file_name)
sheet_0 = workbook.sheet_by_index(0)

workbook = xlwt.Workbook()
sheet = workbook.add_sheet("Parameters")
	

for year in arange(1984, 2008):
	Tmax_0 = list()
	Tmean_0 = list()
	Tmin_0 = list()
	for i in arange(1, sheet_0.nrows):
		if sheet_0.cell_value(i, 1)>=year and sheet_0.cell_value(i, 1)<=year+4:
			
			if sheet_0.cell_value(i, 3)<0:
				pass
			else:
				Tmax_0.append(sheet_0.cell_value(i, 3)/1.0)
			   
			if sheet_0.cell_value(i, 4)<0:
				#print sheet_0.cell_value(i, 6)
				pass
			else:
				Tmean_0.append(sheet_0.cell_value(i, 4)/1.0)
				
			if sheet_0.cell_value(i, 5)<0:
				pass
			else:
				Tmin_0.append(sheet_0.cell_value(i, 5)/10.0)
	#print Tmax_0
	sheet.write(year-1983, 0, year) # row, column, value
	sheet.write(year-1983, 1, (max(np.asarray(Tmax_0))-32)*5/9)
	sheet.write(year-1983, 2, (mean(np.asarray(Tmean_0))-32)*5/9)
	sheet.write(year-1983, 3, (min(np.asarray(Tmin_0))-32)*5/9)
	Tmax_0 = None
	Tmean_0 = None
	Tmin_0 = None	 
	
os.chdir("E:\\Haiyong_File_AL\\p2\\weather_station_air_temperature")
workbook.save("Chestertown.xls")
'''



######## Statistics for the slope larger than zero ##################################
# this code is used to find out the percentage of slopes larger than zero in the Chesapeake Bay
'''
slope_file = "F:\\slope_new_clip.tif"
g = gdal.Open(slope_file)
cols = g.RasterXSize
rows = g.RasterYSize
slope = g.ReadAsArray(0, 0, cols, rows)
sum(slope>0)/sum(slope>-10.0) 
'''
