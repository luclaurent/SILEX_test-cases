from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

#f=open('results/cavity_acou3D_struc_3D_v3_results.frf','rb')
f=open('results/cavity_acou3D_struc_3D_v3_nowall_results.frf','rb')
#f=open('results/cavity_acou3D_struc_3D_v3_4.00E+00_3.00E+00_2.49E+00_1.00E+00_results.frf','rb')
frf_0=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_0[0],10*log10(frf_0[1]/prefsquare),'r-',label='test 0', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid(True)
pl.legend(loc=4)

pl.show()
