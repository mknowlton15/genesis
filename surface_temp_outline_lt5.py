import sys, os
import arcpy
import numpy

class Point:
	
	def __init__(self, x, y)
		self.x = doy
		self.y = s
		
def cleanMakeDir(out_dir,mode=0777):
	if os.path.exists(out_dir):
		print "Output Directory exist!"
		removeall(out_dir)   #this erases the contents of the directory
		print "Cleaned Output Directory!"
	else:
		os.mkdir(out_dir)
		os.chmod(out_dir,mode)
		print "Created Output Directory!"
def listFiles(dir):
	rootdir = dir
	for root, subFolders, files in os.walk(rootdir):
		for file in files:
			yield os.path.join(root,file)
	return
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
def readVapor(afile): #returns a dictionary object with scene/values pairs
	f=open(afile)
	lines=f.readlines()
	lines.pop(0)
	adict=dict()
	for line in lines:
		line=line.rstrip()
		words=line.split(',')
		adict[words[0][:16]]=words[1]
	return adict
	
wd=sys.argv[1]
outdir=sys.argv[2]
os.chdir(wd)
cleanMakeDir(outdir)
#set up constants matrix
c = [[.14714, -.15583, 1.1234], [-1.1836, -.37607, -.52894], [-.04554, 1.8719, -.39071]]
dic=readVapor('E:\\Data\\p14r34_Utility\\atmos_water_scene.csv')
#make a list of the .tif's
tif_list=list()
tree=listFiles(wd)
for t in tree:
	if os.path.splitext(t)[1]=='.tif':
		tif_list.append(os.path.basename(t))
tif_list.sort()
#retrieve water data for each scene desired and create water matrix
for tif in tif_list:
	# open the raster
	band_name=os.path.join(tif,'Band_1') #raster object
	arr=arcpy.RasterToNumPyArray(band_name)
	#lowerLeft = arcpy.Point(band_name.extent.XMin,band_name.extent.YMin)
	#cellSize = band_name.meanCellWidth
	outpath=os.path.join(outdir,os.path.splitext(tif)[0]+'_cor.tif')
	w=float(dic[os.path.splitext(tif)[0]]) #get the water content values from the csv
	b = [[w*w], [w], [1]]
	# #do the matrix math to get the phi values
	m = numpy.dot(c, b)	
	# #find Tsen value using landsat constants and the thermal band radiance
	ls = (arr*0.055158)+1.2378 #numbers for before 4/2/2007 (after that date .055376 and 1.18)
	ts = 1260.56/(numpy.log((607.76/(ls))+1))
	# #use Tsen just calculated and thermal band radiance to calculate planck's function factors
	by = 1256
	y = (ts*ts)/(by*ls)
	d = ts-((ts*ts)/by)
	#for each scene calculate the temperature using the calculated phi values and the radiance values at each pixel in the mask
	s = y*((1/.96)*(((m[0])*ls)+m[1])+m[2])+d
	newRaster = arcpy.NumPyArrayToRaster(s)
	newRaster.save(outpath)
	arr=None
	w=None
	b=None
	m=None
	ls=None
	ts=None
	y=None
	d=None
	s=None
	band_name=None