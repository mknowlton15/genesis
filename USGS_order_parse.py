#Steven Guinn
#1/14/2013
#University of Maryland Center for Environmental Science
#Appalachian Laboratory			http://www.al.umces.edu/ 
#301 689 7146
#301 Braddock Rd	Frostburg MD 21532

import os,urllib,sys
file_path = sys.argv[1]
outdir=sys.argv[2]
f=open(file_path,'r')
line_list=f.readlines()
path_list=list()
for line in line_list:
	words=line.split('"')
	if len(words)>1:
		for w in words:
			if w.find('.tar.gz')!= -1:
				path_list.append(w)
print len(path_list)," images to be downloaded"
for path in path_list:
	if not os.path.isfile(os.path.join(outdir,os.path.basename(path.strip()))):
		print os.path.basename(path.strip())
		urllib.urlretrieve(path.strip(),os.path.join(outdir,os.path.basename(path.strip())))
	else:
		print os.path.basename(path.strip()), "exists!"