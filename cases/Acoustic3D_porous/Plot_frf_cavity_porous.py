from numpy import *
import string
import time
import scipy
#from scipy.sparse import *
#from scipy import *
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen
#from numpy.linalg import solve, norm
import os
#from scipy.linalg import decomp
import pylab as pl
import pickle

f=open('cavity1_results.frf','rb')
frf_cavity1=pickle.load(f)
f.close()

#frf_romain=loadtxt('cavity1_resultats_romain.res')
frf_romain=loadtxt('resultat.res')

prefsquare=20e-6**2

pl.figure(1)
pl.plot(frf_cavity1[0],10*log10(frf_cavity1[1]/prefsquare),'ok-',label='cavity 1', linewidth=2)
pl.plot(frf_romain[:,0],frf_romain[:,1],'or-',label='Romain', linewidth=2)
#pl.axis([100.0, 700.0, 210, 260])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()

