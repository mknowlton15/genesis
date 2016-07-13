# -*- coding: utf-8 -*-
#Steven Guinn
#05/03/2016
#University of Maryland Center for Environmental Science
#Appalachian Laboratoryhttp://www.al.umces.edu/ 
#301 689 7146
#301 Braddock RdFrostburg MD 215325


from copy import deepcopy
import multiprocessing as mp
from multiprocessing import sharedctypes
import ctypes
import numpy as np
import math

#************************************** Variable Declarations *************************************
########################################    Constants      ########################################
TORAD = math.pi/180
INT16_NODATA = -32768
FLOAT32_NODATA = -3.402823466e+38
ACCESS_MODE = 0777
SCALEF = .0001
MEST=np.mat('0.191943;0.540254;116.499403;16.876135;296.285736;21.686111;0.000').T 		# Average parameters for the 'NBY' population for intial estimates
########################################    Structures      #######################################
# Create a dictionary to convert data types (ctypes to np.dtype)
ctypes_to_numpy = {ctypes.c_char   : np.dtype(np.uint8),ctypes.c_wchar  : np.dtype(np.int16),ctypes.c_byte   : np.dtype(np.int8),ctypes.c_ubyte  : np.dtype(np.uint8),ctypes.c_short  : np.dtype(np.int16),ctypes.c_ushort : np.dtype(np.uint16),ctypes.c_int: np.dtype(np.int32),ctypes.c_uint   : np.dtype(np.uint32),ctypes.c_long   : np.dtype(np.int64),ctypes.c_ulong  : np.dtype(np.uint64),ctypes.c_float  : np.dtype(np.float32),ctypes.c_double : np.dtype(np.float64)}

# Create a dictionary to convert data types (np.dtype to ctypes)
numpy_to_ctypes = dict(zip(ctypes_to_numpy.values(), ctypes_to_numpy.keys()))
#_numpy_dtype_to_unicode={str(np.dtype(np.int8)):'b',str(np.dtype(np.uint8)):'u1',str(np.dtype(np.float64)):'d',str(np.dtype(np.int16)):'i'}
#************************************** Function Definitions **************************************
########################################    OS File operations     ################################
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
	MEST=np.mat('0.10482;0.95465;124.56483;21.26626;310.78949;17.82388; 0.00126').T 	
	mref=MEST   
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
	dresid = (vegmat-fwdlogit(MEST,doymat))
	#set up the stopping values
	misfit = np.sqrt(np.power(Wd*dresid.T,2)).sum()
	oldmisfit=2000
	pert=0
	##this is where the start of the loop goes
	try:
		while (pert < 2) or ((pert < 101) and (not misfitSmall(misfit,dresid)) and (not misfitDelta(misfit,oldmisfit))):
			G=calcG(MEST,doymat)
			if G == None:
				return mdefault.T
			GtWgG = G.T*Wg*G
			GtWgGdmp = GtWgG+C        
			dresid = (vegmat-fwdlogit(MEST,doymat))
			GtWg = G.T*Wg        
			GtWgdresid = (GtWg*dresid.T)+(C*(mref-MEST).T)
			mpert = GtWgGdmp.I * GtWgdresid        
			MEST = MEST+mpert.T     
			Wg=l1_solution(vegmat,MEST,doymat,Wd)
			oldmisfit=misfit
			misfit = np.sqrt(np.power(Wd*dresid.T,2)).sum()
			pert +=1
		return np.hstack([MEST,np.array([misfit]).reshape((1,1)),np.array([pert]).reshape((1,1)),np.array([doymat.shape[1]]).reshape((1,1))])
	except numpy.linalg.linalg.LinAlgError:
		print "whoops"
		return np.mat('0.0;0.0;0.0;0.0;0.0;0.0;0.0; 0.0; 0.0; 0.0').T