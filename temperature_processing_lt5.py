# -*- coding: utf-8 -*-
#Steven Guinn
#02/24/2016
#University of Maryland Center for Environmental Science
#Appalachian Laboratory            http://www.al.umces.edu/ 
#301 689 7146
#301 Braddock Rd    Frostburg MD 215325

import sys,os
# Import arcpy module
import arcpy

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

wd=sys.argv[1]
outdir=sys.argv[2]
os.chdir(wd)

# get the list of compressed archives
gz_list=list()
tree=listFiles(wd)
cleanMakeDir(outdir)
for t in tree:
	if os.path.splitext(t)[1]=='.gz':
		gz_list.append(os.path.basename(t))
gz_list.sort()
_temp_dir = os.path. join(wd,'_temp')
for gz in gz_list:
	cleanMakeDir(_temp_dir)
	_command='E:\\Data\\z\\7z.exe x '+gz+' -o'+_temp_dir
	print _command
	os.system(_command)
	tar_name = os.path.join(_temp_dir,os.path.splitext(gz)[0])
	_command='E:\\Data\\z\\7z.exe x '+tar_name+' -o'+outdir+' *cfmask.tif'
	print _command
	os.system(_command)
	tif_list=list()
	tiffs=listFiles(_temp_dir)
	for t in tiffs:
		if os.path.splitext(t)[1]=='.tif':
			tif_list.append(t)
	scene_name_tif=os.path.join(outdir,(os.path.basename(gz).split('-')[0]+'.tif'))
	# arcpy.CompositeBands_management(tif_list[0]+';'+tif_list[1]+';'+tif_list[2], scene_name_tif)
	