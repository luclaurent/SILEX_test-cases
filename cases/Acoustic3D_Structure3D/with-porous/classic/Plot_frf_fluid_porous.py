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

#f=open('results/cavity1_fluid_only_results.frf','rb')
#frf_classic_no_damping_cavity1=pickle.load(f)
#f.close()

f=open('results/cavity1_with_porous_20_1_angle_90_results.frf','rb')
frf_classic_with_porous_20_1_angle_90=pickle.load(f)
f.close()

f=open('results/cavity1_no_porous_20_1_angle_90_results.frf','rb')
frf_classic_no_porous_20_1_angle_90=pickle.load(f)
f.close()

f=open('results/cavity1_with_porous_20_1_angle_30_results.frf','rb')
frf_classic_with_porous_20_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity1_with_porous_20_1_angle_60_results.frf','rb')
frf_classic_with_porous_20_1_angle_60=pickle.load(f)
f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_classic_no_porous_20_1_angle_90[0],10*log10(frf_classic_no_porous_20_1_angle_90[1]/prefsquare),'k-',label='Classic, no porous, 20_1, angle 90', linewidth=2)
pl.plot(frf_classic_with_porous_20_1_angle_90[0],10*log10(frf_classic_with_porous_20_1_angle_90[1]/prefsquare),'b-',label='Classic, with porous, 20_1, angle 90', linewidth=2)
pl.plot(frf_classic_with_porous_20_1_angle_30[0],10*log10(frf_classic_with_porous_20_1_angle_30[1]/prefsquare),'r-',label='Classic, with porous, 20_1, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_20_1_angle_60[0],10*log10(frf_classic_with_porous_20_1_angle_60[1]/prefsquare),'y-',label='Classic, with porous, 20_1, angle 60', linewidth=2)
pl.axis([1.0, 200.0, 40, 160])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()

