import sys,os
import csv
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
os.chdir(wd)

tif_list=list()
tree=listFiles(wd)
for t in tree:
	if os.path.splitext(t)[1]=='.tif':
		tif_list.append(os.path.basename(t))
tif_list.sort()
doy=list()
for tif in tif_list:
	doy.append(tif[13:16])
doy.sort()
print doy

with open('filename', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(doy)
	