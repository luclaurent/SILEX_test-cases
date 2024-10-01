from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

f=open('results/cavity_acou3D_struc_3D_v3_1.00E+00_1.00E+00_5.00E-01_8.00E-01_results.frf','rb')
frf_no_reduc=pickle.load(f)
f.close()

#f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_1.00E+00_1.00E+00_5.00E-01_8.00E-01_results.frf','rb')
#frf_reduc=pickle.load(f)
#f.close()

f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient_1.00E+00_1.00E+00_NbModesFluid_20_5.00E-01_NbModesFluid_20_8.00E-01_NbModesFluid_20_results.frf','rb')
frf_reduc_gradient=pickle.load(f)
f.close()

#f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_1.00E+00_1.00E+00_5.00E-01_8.00E-01_results.frf','rb')
#frf_reduc_classic=pickle.load(f)
#f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_no_reduc[0],10*log10(frf_no_reduc[1]/prefsquare),'b-',label='no reduction', linewidth=1)
#pl.plot(frf_reduc[0],10*log10(frf_reduc[1]/prefsquare),'r-',label='reduction', linewidth=1)
pl.plot(frf_reduc_gradient[0],10*log10(frf_reduc_gradient[1]/prefsquare),'r-',label='reduction and gradient', linewidth=1)
#pl.plot(frf_reduc_classic[0],10*log10(frf_reduc_classic[1]/prefsquare),'g-',label='reduction', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
pl.plot(frf_no_reduc[0],frf_no_reduc[2],'b-',label='no reduction, gradient pos. x', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[2],'r-',label='reduction, gradient pos. x', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(3)
pl.plot(frf_no_reduc[0],frf_no_reduc[3],'b-',label='no reduction, gradient pos. y', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[3],'r-',label='reduction, gradient pos. y', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(4)
pl.plot(frf_no_reduc[0],frf_no_reduc[4],'b-',label='no reduction, gradient pos. z', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[4],'r-',label='reduction, gradient pos. z', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(5)
pl.plot(frf_no_reduc[0],frf_no_reduc[5],'b-',label='no reduction, gradient radius', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[5],'r-',label='reduction, gradient radius', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(6)
pl.plot(frf_no_reduc[0],frf_no_reduc[2]-frf_reduc_gradient[2],'b-',label='gradient difference: no reduction - reduction, pos. x', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient difference in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.show()

