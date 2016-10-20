
'''
This script implements the "Model B" analyses from:

	Pataky TC, Koseki M, Cox PG (2016) Probabilistic biomechanical
	finite element simulations: whole-model classical hypothesis testing
	based on upcrossing geometry. PeerJ Computer Science. (in press)

!!!NOTE!!!
In order to run this script you must modify the "path2febio" variable below.


There are three sub models:
	modelB0.feb
	modelB1.feb
	modelB2.feb

Model B0 contains an indenter with a perfectly flat surface.
Models B1 and B2 add random noise to the indenter surface's z coordinates,
with the randomness amplitude biased toward one side (B1) or the other (B2).

The main procedures implemented in this script include:
	1. FEBio model file manipulation (material stiffness distribution)
	2. Simulation (using the FEBio solver version 2.5)
	3. FEBio results parsing

Software dependencies:
	(other versions of the following packages should also work)
	Non-Python software:
		FEBio 2.4     (febio.org)

	Python software:
		Python 2.7       (python.org)
		NumPy 1.10       (scipy.org)
		Matplotlib 1.5   (matplotlib.org)

This script runs in 2.1 minutes on:
     Mac OS X 10.11, 2.7 GHz 12-core Intel Xeon E5, 32 GB 1866 MHz DDR3 ECC
It was also tested on Windows 7 32-bit

Version 0.1   (2016.09.02)
'''




import os
import numpy as np
from matplotlib import pyplot
from xml.etree.ElementTree import ElementTree



#---------------------------------------------------------------#
### USER PARAMETERS ###
### Specify the path to the FEBio binary executable:
path2febio = '/Applications/febio/v2.5.0/bin/FEBio2'
### Default executable locations are:
###    Mac OS X:   /Applications/FEBio2.X.X/bin/FEBio2
###    Windows:    C:\\Program Files\\febio-2.X.X\\bin\\FEBio2.exe
#---------------------------------------------------------------#





def parse_logfile(fnameLOG, nElements=2048):
	'''
	Reads the strain and stress tensor fields from the final data record in an FEBio log file.

	Arguments:
	fname -- full path to the log file

	Returns:
	A -- an (nElement x 12) array containing the strain and stress tensor fields
	'''
	with open(fnameLOG, 'r') as fid:
		lines = fid.readlines()
	ind = []
	for i,s in enumerate(lines):
		if s.startswith('Data Record'):
			ind.append(i)
	i  = ind[-1] + 5
	A  = [s.strip().split(' ')[1:]   for s in lines[i:i+nElements]]
	return np.asarray(A, dtype=float)



def simulate(fname0, k, fname1, silent=False):
	'''
	Simulate the model given a material value "k".
	
	Arguments:
	fname0 -- template FEB file
	E -- stiffness profile:  a (101,) numpy array
	fname1 -- temporary FEB file to be written and simulated
	
	Returns:
	stress -- von Mises stress field of the indented surface (top elements only)
	'''
	write_model(fname0, k, fname1)
	### simulate:
	command   = '"%s" -i %s' %(path2febio, fnameFEB1)
	if silent:
		command += ' -silent'
	os.system( command )
	### parse output:
	fnameLOG  = os.path.splitext(fnameFEB1)[0] + '.log'
	A         = parse_logfile(fnameLOG)
	stress    = tensor2effective(A[:,6:])
	### reshape into an image:
	nex,ney   = 32,32   #numbers of elements in the x and y directions
	n         = nex*ney #total number of elements in the top layer
	return stress[-n:].reshape([nex,ney])




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


def write_model(fnameFEB0, k, fnameFEB1):
	'''
	Write a new FEB file with a new Mooney-Rivlin parameter "k"
	
	Arguments:
	fnameFEB0 : template FEB file
	k : scalar
	fnameFEB1 : FEB file to be written (will be overwritten if it exists)
	'''
	tree      = ElementTree()
	tree.parse(fnameFEB0)
	root      = tree.getroot()
	root      = tree.getroot()
	node      = root.find('Material/material/k')
	node.text = str(k)
	tree.write(fnameFEB1, encoding='ISO-8859-1')





#(0) Run model:
dir0      = os.path.split(__file__)[0]
model     = 1  #0, 1 or 2  (0=flat contact surface,  1&2=jagged contact surfaces)
fnameFEB0 = os.path.join( dir0 , 'modelB%d.feb' %model)
fnameFEB1 = os.path.join( dir0 , 'temp.feb')
k         = 800
S         = simulate(fnameFEB0, k, fnameFEB1, silent=False) #silent=True will silence FEBio output


#(1) Plot the distribution:
pyplot.close('all')
fig   = pyplot.figure(figsize=(6,4))
pyplot.get_current_fig_manager().window.move(0, 0)
ax    = pyplot.axes()
ax.imshow(S, interpolation='nearest')
cb    = pyplot.colorbar(mappable=ax.images[0])
cb.set_label('von Mises stress  (Pa)')
pyplot.show()






