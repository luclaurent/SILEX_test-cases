from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

# reference   
f=open('cube_acou_h_10/cube_acou_tet4_results.frf','rb')
frf_tet4_ref=pickle.load(f)
f.close()

f=open('cube_acou_h_10/cube_acou_tet10_results.frf','rb')
frf_tet10_ref=pickle.load(f)
f.close()

#
f=open('cube_acou_tet4_results.frf','rb')
frf_tet4=pickle.load(f)
f.close()

f=open('cube_acou_tet10_results.frf','rb')
frf_tet10=pickle.load(f)
f.close()


prefsquare=20e-6**2

pl.figure(1)
pl.plot(frf_tet4_ref[0],10*log10(frf_tet4_ref[1]/prefsquare),'r-',label='tet4 reference', linewidth=1)
pl.plot(frf_tet4[0],10*log10(frf_tet4[1]/prefsquare),'b-',label='tet4', linewidth=1)
pl.plot(frf_tet10_ref[0],10*log10(frf_tet10_ref[1]/prefsquare),'g-',label='tet10 reference', linewidth=1)
pl.plot(frf_tet10[0],10*log10(frf_tet10[1]/prefsquare),'m-',label='tet10', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()


