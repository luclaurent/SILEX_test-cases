from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

# xfem + POROUS + flexible structure + reduction
f=open('results/cavity9_with_porous_air_flexible_structure_CB_ang90_gradient_lx3_2000_results.frf','rb')
frf_lx3_2000=pickle.load(f)
f.close()

f=open('results/cavity9_with_porous_air_flexible_structure_CB_ang90_gradient_lx3_2001_results.frf','rb')
frf_lx3_2001=pickle.load(f)
f.close()

f=open('results/cavity9_with_porous_air_flexible_structure_CB_ang90_gradient_lx3_2002_results.frf','rb')
frf_lx3_2002=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_lx3_2000[0],10*log10(frf_lx3_2000[1]/prefsquare),'or-',label='lx3 2.000', linewidth=1)
pl.plot(frf_lx3_2001[0],10*log10(frf_lx3_2001[1]/prefsquare),'og-',label='lx3 2.001', linewidth=1)
pl.plot(frf_lx3_2002[0],10*log10(frf_lx3_2002[1]/prefsquare),'ob-',label='lx3 2.002', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)


#fi=28
#print(str(frf_lx3_2000[0][fi])+' Hz')

Dp_Dtheta_diff = (frf_lx3_2002[1]-frf_lx3_2000[1])/(2.002-2.000)
Dp_Dtheta_grad = frf_lx3_2001[2]

Dp_Dtheta_diff/Dp_Dtheta_grad

pl.figure(2)
pl.plot((Dp_Dtheta_diff-Dp_Dtheta_grad)/Dp_Dtheta_diff)
pl.figure(3)
pl.plot(Dp_Dtheta_diff)
pl.figure(3)
pl.plot(Dp_Dtheta_grad)

pl.show()
