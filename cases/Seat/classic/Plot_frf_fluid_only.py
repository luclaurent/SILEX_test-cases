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

f=open('results/classic_classic_fluid_cavity_results_no_damping.frf','rb')
frf_classic_no_damping=pickle.load(f)
f.close()

f=open('results/classic_classic_fluid_cavity_results_with_damping.frf','rb')
frf_classic_with_damping=pickle.load(f)
f.close()

f=open('results/classic_classic_fluid_cavity_with_modes_results_with_damping.frf','rb')
frf_classic_with_modes_with_damping=pickle.load(f)
f.close()


f=open('results/classic_pgd_fluid_cavity.frf','rb')
frf_classic_pgd=pickle.load(f)
f.close()
frf_classic_pgd[1]=scipy.array(frf_classic_pgd[1])

f=open('results/classic_pgd_fluid_cavity_with_modes.frf','rb')
frf_classic_pgd_with_modes=pickle.load(f)
f.close()
frf_classic_pgd_with_modes=scipy.array(frf_classic_pgd_with_modes)

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_classic_no_damping[0],10*log10(frf_classic_no_damping[1]/prefsquare),'k-',label='Classic, no damping', linewidth=2)
pl.plot(frf_classic_pgd[0],10*log10(frf_classic_pgd[1]/prefsquare),'og-',label='Classic, PGD', linewidth=2)
pl.plot(frf_classic_with_damping[0],10*log10(frf_classic_with_damping[1]/prefsquare),'b-',label='Classic, with damping', linewidth=2)
pl.plot(frf_classic_pgd_with_modes[0],10*log10(frf_classic_pgd_with_modes[1]/prefsquare),'or-',label='Classic, PGD, modal projection', linewidth=2)
pl.plot(frf_classic_with_modes_with_damping[0],10*log10(frf_classic_with_modes_with_damping[1]/prefsquare),'y-',label='Classic, with damping, modal projection', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()

