
'''
This script implements ***part of the "Model C" analyses from:

	Pataky TC, Koseki M, Cox PG (2016) Probabilistic biomechanical
	finite element simulations: whole-model classical hypothesis testing
	based on upcrossing geometry. PeerJ Computer Science. (in review)

See README.txt for a list of steps required to get started with this
script.

*** Since this model takes a long time to solve (~18 minutes),
this scipts implements only a single iteration. Users wishing to
reproduce the results from the paper should run this script
iteratively using the material values from Table 1, save the
resulting strain / stress fields for each iteraton to disk, then
apply a permutation approach to build a primary probability density
similar to the "ttest_nonparametric" procedure in "modelA.py".
Since this is a two-sample comparison, there are two main differences
between these analyses in the Model A analyses:
1) The t statistic should be computed according to typical two-sample
definition.
2) Permutation should randomly assign the 20 observations to two groups
(10 in each group), and not conduct sign permutation as in "modelA.py".
Python functions implementing these two changes is provided as a
comment at the end of this script.

The results can be visualized using a variety of 3D visualization
packages like VTK. One option with FEB-file compatibility is
FEBlender (https://github.com/0todd0000/feblender) which uses
Blender (blender.org) for 3D finite element field visualizations.

The main procedures implemented in this script include:
	1. FEBio model file manipulation (material parameter)
	2. Simulation (using the FEBio solver version 2.4.2)
	3. FEBio results parsing

Version 0.1   (2016.08.30)
'''



import os
from xml.etree.ElementTree import ElementTree
import numpy as np



#---------------------------------------------------------------#
### USER VARIABLES ###
### Specify the path to the FEBio binary executable:
path2febio = '/Applications/febio/v2.4.2/bin/FEBio2'
### Default executable locations are:
###    Mac OS X:   /Applications/FEBio2.X.X/bin/FEBio2
###    Windows:    C:\\Program Files\\febio-2.X.X\\bin\\FEBio2.exe
### Specify paths to the FEB files:
#   Original model file
fnameFEB   = '/tmp/hip_n10rb.feb'
#   Temporary model file containing adjusted material parameters 
fnameTEMP  = '/tmp/temp.feb'
### Set material parameter:
K          = 1350
#---------------------------------------------------------------#



def parse_logfile(fname):
	'''
	Reads the strain and stress tensor fields from the final data record in an FEBio log file.
	
	Arguments:
	fname -- full path to the log file
	
	Returns:
	A -- an (nElement x 12) array containing the strain and stress tensor fields
	'''
	with open(fnameLOG, 'r') as fid:
		lines = fid.readlines()
	for i,s in enumerate(lines):
		if s.startswith('Data Record'):
			break
	i += 5
	A  = [s.strip().split(' ')[1:]   for s in lines[i:i+149300]]
	return np.asarray(A, dtype=float)


def tensor2effective(Y):
	'''
	Compute effective strain field from a strain tensor field.
	(Or compute von Mises stress field from a stress tensor field)
	
	Arguments:
	Y -- a (101,6) numpy array containing the tensor field
	
	Returns:
	y -- effective strain field (or von Mises stress field): a (101,) numpy array
	'''
	x0,x1,x2, a,b,c = Y.T
	s = (x0-x1)**2 + (x0-x2)**2 + (x1-x2)**2 + 6*(a*a + b*b + c*c)
	return (0.5*s)**0.5


#(0) Adjust material parameter value:
print('Adjusting material parameters...')
### parse template FEB file:
tree        = ElementTree()
tree.parse(fnameFEB)
root        = tree.getroot()
### adjust material:
materials   = root.find('Material')
mat0,mat1   = materials[2:]  #the two cartilage materials
k0,k1       = mat0.find('k'), mat1.find('k')
k0.text     = str(K)
k1.text     = str(K)
### write to disk:
tree.write(fnameTEMP, encoding="ISO-8859-1")


#(1) Simulate:
print('Running simulation...')
dirWORK     = os.path.split(fnameTEMP)[0]
os.chdir(dirWORK)
os.system('%s %s'%(path2febio, fnameTEMP))


#(2) Parse log file:
print('Parsing log file...')
fnameLOG    = os.path.splitext(fnameTEMP)[0] + '.log'
A           = parse_logfile( fnameLOG )
strain      = tensor2effective(A[:,:6])
stress      = tensor2effective(A[:,6:])





'''
Python functions for conducting a two-sample permutation test:


import itertools
from math import sqrt
import numpy as np


def tstat2(YA, YB):
	mA,mB  = YA.mean(axis=0), YB.mean(axis=0)
	sA,sB  = YA.std(ddof=1, axis=0), YB.std(ddof=1, axis=0)
	n      = YA.shape[0]
	s      = np.sqrt(0.5*(sA*sA + sB*sB))
	t      = (mA-mB) / ( s *sqrt(2.0/n))
	return t

def tstat2_labels(Y, ones):
	labels     = np.array([0]*Y.shape[0])
	labels[list(ones)] = 1
	YA,YB      = Y[labels==0], Y[labels==1]
	return tstat2(YA,YB)

def ttest2_nonparam(YA, YB, alpha=0.05, nIterations=-1):
	nA,nB    = YA.shape[0], YB.shape[0]
	n        = nA + nB
	if YA.ndim==1:
		Y    = np.hstack([YA,YB])
	else:
		Y    = np.vstack([YA,YB])
	n        = Y.shape[0]
	t0       = tstat2(YA, YB)
	if nIterations==-1:
		ONES     = list(itertools.combinations(range(n), n/2))
		T        = [tstat2_labels(Y, ones).max()  for ones in ONES]
	else:
		T        = []
		for i in range(nIterations):
			ones = np.random.permutation(range(n))[:nA]
			t    = tstat2_labels(Y, ones)
			T.append( t.max() )
	tstar    = np.percentile(T, 100*(1-alpha))
	return t0,tstar
'''
