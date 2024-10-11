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

f=open('results/cavity1_fluid_only_results.frf','rb')
frf_classic_no_damping_cavity1=pickle.load(f)
f.close()

f=open('results/cavity2_results.frf','rb')
frf_classic_no_damping_cavity2=pickle.load(f)
f.close()


f=open('results/cavity2_damping_results.frf','rb')
frf_classic_damping_cavity2=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_classic_no_damping_cavity1[0],10*log10(frf_classic_no_damping_cavity1[1]/prefsquare),'k-',label='Classic, no damping, cavity 1', linewidth=2)
pl.plot(frf_classic_no_damping_cavity2[0],10*log10(frf_classic_no_damping_cavity2[1]/prefsquare),'r-',label='Classic, no damping, cavity 2', linewidth=2)
pl.plot(frf_classic_damping_cavity2[0],10*log10(frf_classic_damping_cavity2[1]/prefsquare),'b-',label='Classic, damping, cavity 2', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()

