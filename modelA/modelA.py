
'''
This script implements the "Model A" analyses from:

	Pataky TC, Koseki M, Cox PG (2016) Probabilistic biomechanical
	finite element simulations: whole-model classical hypothesis testing
	based on upcrossing geometry. PeerJ Computer Science. (in review)

!!!NOTE!!!
In order to run this script you must modify the "path2febio" variable on
Line 53 below.

The main procedures implemented in this script include:
	1. FEBio model file manipulation (material stiffness distribution)
	2. Simulation (using the FEBio solver version 2.5)
	3. FEBio results parsing
	4. Non-parametric one-sample t test on the results
	5. Results plotting

Software dependencies:
	(other versions of the following packages should also work)
	Non-Python software:
		FEBio 2.5.0      (febio.org)

	Python software:
		Python 2.7       (python.org)
		NumPy 1.10       (scipy.org)
		SciPy 0.17       (scipy.org)
		Matplotlib 1.5   (matplotlib.org)

This script runs in 1.3 seconds on:
     Mac OS X 10.11, 2.7 GHz 12-core Intel Xeon E5, 32 GB 1866 MHz DDR3 ECC
It was also tested on Windows 7 (32-bit) with similar calculation speeds.

Version 0.1   (2016.07.01)
'''






import os,itertools
import numpy as np
from scipy import ndimage
from xml.etree.ElementTree import ElementTree
from matplotlib import pyplot


#---------------------------------------------------------------#
### USER VARIABLES ###
### Specify the path to the FEBio binary executable:
path2febio = '/Applications/febio/v2.5.0/bin/FEBio2'
### Default executable locations are:
###    Mac OS X:   /Applications/FEBio2.X.X/bin/FEBio2
###    Windows:    C:\\Program Files\\febio-2.X.X\\bin\\FEBio2.exe
#---------------------------------------------------------------#


def check_paths(path2febio, fnameCSV, fnameFEB0):
	'''
	Check that the necessary executable and files exist.
	'''
	if not os.path.exists(path2febio):
		raise( IOError('The specified path2febio does not exist:\n   %s' %path2febio)  )
	if not os.path.exists(fnameCSV):
		raise( IOError('The specified fnameCSV does not exist:\n   %s\n\nThis script should be run in the folder containing the "stiffness_profiles.csv" file\n' %fnameCSV)  )
	if not os.path.exists(fnameFEB0):
		raise( IOError('The specified fnameFEB0 does not exist:\n   %s\n\nThis script should be run in the folder containing the "template.feb" file\n' %fnameFEB0)  )
	


def cluster_integral(z, thresh, i):
	'''
	Compute the integral of a supratheshold cluster using a trapezoidal approximation.
	
	Arguments:
	z -- test statistic field
	thresh -- scalar value used to threshold the test statistic field
	i -- binary field specifying a specific cluster location  (multiple clusters may exist above thresh)
	'''
	if i.sum()==1:
		x = z[i] - thresh
	else:
		x = np.trapz(  z[i]-thresh  )
	return x


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
	A  = [s.strip().split(' ')[1:]   for s in lines[i:i+101]]
	return np.asarray(A, dtype=float)


def plot_stats_results(ax, x, t, tCrit):
	'''
	Plot the results of a non-parametric field-wide one-sample t test.
	
	Arguments:
	ax -- a Maplotlib axes object (where the results will be plotted)
	x -- field element labels
	t -- test statistic field
	tCrit -- critical test statistic value
	'''
	ax.plot(x, t, color='b', lw=3)
	ax.axhline(tCrit, color='r', linestyle='--')
	ax.axhline(-tCrit, color='r', linestyle='--')
	ax.axhline(0, color='k', linestyle='-', lw=0.5)


def simulate(fname0, E, fname1):
	'''
	Simulate the model given a stiffness profile "E"
	
	Arguments:
	fname0 -- template FEB file
	E -- stiffness profile:  a (101,) numpy array
	fname1 -- temporary FEB file to be written and simulated
	
	Returns:
	strain -- effective strain field: a (101,) numpy array
	stress -- von Mises stress field: a (101,) numpy array
	'''
	write_model(fname0, E, fname1)
	### simulate:
	dir0      = os.path.split(fname0)[0]
	os.chdir(dir0)
	os.system('"%s" -i %s'%(path2febio, fnameFEB1))
	### parse output:
	A         = parse_logfile(fnameLOG)
	strain    = tensor2effective(A[:,:6])
	stress    = tensor2effective(A[:,6:])
	return strain, stress


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


def ttest_nonparametric(x, mu, alpha=0.05):
	'''
	Conduct a non-parametric field-wide one-sample t test (two-tailed).
	
	Arguments:
	x -- a (101,N) numpy array containing N observations of 101-element scalar fields
	mu -- a (101,) numpy array representing the datum to which the observations in "x" will be compared
	alpha -- type I error rate
	
	Returns:
	t0 -- a (101,) numpy array containing the test statistic field
	tCrit -- a scalar representing the critical threshold (at a type I error rate of alpha)
	p -- probability values for clusters which survive the "tCrit" threshold;  if no regions of "t0" exceed "tCrit" then "p" will be np.nan
	'''
	### preliminaries:
	y          = x.T - mu     #datum-corrected observations
	n          = y.shape[0]   #number of observations
	sqrtN      = n**0.5       #to be used in test statistic computations below
	### compute original test statistic field:
	t0         = y.mean(axis=0) / y.std(ddof=1, axis=0) * sqrtN
	### build primary permutation PDF (using the maximum t value):
	LABELS     = list(itertools.product([0,1], repeat=n))  #all permutations of (+) or (-) labels
	T          = []  #initialize the primary PDF
	for labels in LABELS:
		signs  = -2*np.array(labels) + 1  #sequence of -1 and +1
		yy     = (y.copy().T*signs).T     #sign-permuted observations
		t      = yy.mean(axis=0) / yy.std(ddof=1, axis=0) * sqrtN
		T.append(t.max())
	T          = np.array(T)
	### compute critical threshold:
	tCrit      = np.percentile(T, 100*(1-alpha))
	### build secondary permutation PDF (cluster size):
	M          = []  #initialize the secondary PDF
	for labels in LABELS:
		signs  = -2*np.array(labels) + 1
		yy     = (y.copy().T*signs).T
		t      = yy.mean(axis=0) / yy.std(ddof=1, axis=0) * sqrtN
		### compute maximum cluster integral:
		L,nC   = ndimage.label( np.abs(t) > tCrit )
		m      = 0
		if nC>0:
			m  = [cluster_integral(np.abs(t), tCrit, L==(i+1))   for i in range(nC)]
			m  = max(m)
		M.append(m)
	M          = np.array(M)
	### compute p values:
	L,nC   = ndimage.label( np.abs(t0) > tCrit )
	if nC>0:
		m0 = [cluster_integral(np.abs(t0), tCrit, L==(i+1))   for i in range(nC)]
		p  = [(M > m).mean()   for m in m0]
		p  = [max(pp, 1.0/M.size)  for pp in p]  #correct for the case when the original test statistic field contains the largest clusters from the secondary PDF
	else:
		p = np.nan
	return t0,tCrit,p


def write_model(fnameFEB0, E, fnameFEB1):
	'''
	Write a new FEB file based on a template and the stiffness profile "E"
	
	Arguments:
	fnameFEB0 -- template FEB file
	E -- stiffness profile:  a (101,) numpy array
	fnameFEB1 -- FEB file to be written (will be overwritten if it exists)
	'''
	tree      = ElementTree()
	tree.parse(fnameFEB0)
	root      = tree.getroot()
	materials = root.findall('Material/material')
	for mat,e in zip(materials,E):
		mat.find('E').text = str(e)
	tree.write(fnameFEB1, encoding='ISO-8859-1')








#---------------------------------------------------------------#
# MAIN SCRIPT
#---------------------------------------------------------------#




#(0) Specify file names:
dir0            = os.path.split(__file__)[0]
fnameCSV        = os.path.join(dir0, 'stiffness_profiles.csv')  #stifness profiles from the paper
fnameFEB0       = os.path.join(dir0, 'template.feb')            #template FEB file
fnameFEB1       = os.path.join(dir0, 'temp.feb')                #temporary FEB file (will be overwritten if it exists)
fnameLOG        = os.path.join(dir0, 'temp.log')                #FEBio log file containing strain/stress fields
check_paths(path2febio, fnameCSV, fnameFEB0)


#(1) Simulate the datum:
E0              = 14e9 * np.ones(101)  #constant stiffness for all elements
strain0,stress0 = simulate(fnameFEB0, E0, fnameFEB1)



#(2) Cycle through all stiffness profiles:
EE              = np.loadtxt(fnameCSV, delimiter=',')
STRAIN,STRESS   = [],[]
for E in EE.T:
	strn,strs   = simulate(fnameFEB0, E, fnameFEB1)
	STRAIN.append(strn)
	STRESS.append(strs)
STRAIN,STRESS   = np.asarray(STRAIN).T, np.asarray(STRESS).T



#(3) Constrain the hypotheses to elements 50-90 (the region of local stiffness change)
i               = np.array([False]*101)
i[50:91]        = True
E0,strain0,stress0,EE,STRAIN,STRESS = [x[i]  for x in [E0,strain0,stress0,EE,STRAIN,STRESS]]



#(4) Run non-parametric tests:
t0,tCrit0,p0    = ttest_nonparametric(EE, E0, alpha=0.05)
t1,tCrit1,p1    = ttest_nonparametric(STRAIN, strain0, alpha=0.05)
t2,tCrit2,p2    = ttest_nonparametric(STRESS, stress0, alpha=0.05)
print( "p values (Young's modulus):  %s" %str(p0))
print( "p values (effective strain): %s" %str(p1))
print( "p values (vonMises stress):  %s" %str(p2))



#(5) Plot results:
pyplot.close('all')
x  = np.arange(101)[i]
### plot stiffnesses:
pyplot.figure(figsize=(8,6))
pyplot.get_current_fig_manager().window.move(0, 0)
ax = pyplot.subplot(221);  ax.plot(x, 1e-9*E0, color='r');      ax.plot(x, 1e-9*EE, color='b' );      ax.set_title('Stiffness profile  (GPa)')
ax = pyplot.subplot(223);  ax.plot(x, 1e6*strain0, color='r');  ax.plot(x, 1e6*STRAIN, color='b' );   ax.set_title('Effective stain  (1e-6)')
ax = pyplot.subplot(224);  ax.plot(x, 1e-3*stress0, color='r'); ax.plot(x, 1e-3*STRESS, color='b' );  ax.set_title('von Mises stress  (kPa)')
ax = pyplot.subplot(221);  ax.legend(['Datum', 'Observation'], fontsize=8)  
pyplot.suptitle('All observations', size=20)
pyplot.show()
### plot statistical results:
pyplot.figure(figsize=(8,6))
pyplot.get_current_fig_manager().window.move(50, 50)
ax = pyplot.subplot(221);  plot_stats_results(ax, x, t0, tCrit0);  ax.set_title('Stiffness profile');  ax.set_ylim(-4, 4)
ax = pyplot.subplot(223);  plot_stats_results(ax, x, t1, tCrit1);  ax.set_title('Effective stain');    ax.set_ylim(-4, 4)
ax = pyplot.subplot(224);  plot_stats_results(ax, x, t2, tCrit2);  ax.set_title('von Mises stress');   ax.set_ylim(-4, 4)
ax = pyplot.subplot(223);  ax.legend(['Test statistic', 'Critical threshold'], fontsize=8)
pyplot.suptitle('One-sample t tests', size=20)
pyplot.show()


