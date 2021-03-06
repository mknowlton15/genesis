import sys, os
import arcpy
import numpy

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
c = [[.04019, .02916, 1.01523], [-.38333, -1.50294, .20324], [.00918, 1.36072, -.27514]]
dic=readVapor('E:\\Data\\p14r34_Utility\\atmos_water_scene_lc8.csv')
#make a list of the .tif's
tif_list=list()
tree=listFiles(wd)
for t in tree:
	if os.path.splitext(t)[1]=='.TIF':
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
	ls = (arr*0.0003342)+.1
	ts = 1201.14/(numpy.log((480.89/(ls))+1))
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
	newRaster=None