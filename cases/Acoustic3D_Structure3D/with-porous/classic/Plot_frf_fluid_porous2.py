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


f=open('results/cavity2_with_porous_20_1_angle_90_results.frf','rb')
frf_classic_with_porous_20_1_angle_90=pickle.load(f)
f.close()

#f=open('results/cavity2_no_porous_20_1_angle_90_results.frf','rb')
#frf_classic_no_porous_20_1_angle_90=pickle.load(f)
#f.close()

f=open('results/cavity2_with_porous_20_1_angle_30_results.frf','rb')
frf_classic_with_porous_20_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_2_angle_30_results.frf','rb')
frf_classic_with_porous_20_2_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_30_2_angle_30_results.frf','rb')
frf_classic_with_porous_30_2_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_30_1_angle_30_results.frf','rb')
frf_classic_with_porous_30_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_30_3_angle_30_results.frf','rb')
frf_classic_with_porous_30_3_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_40_1_angle_30_results.frf','rb')
frf_classic_with_porous_40_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_40_2_angle_30_results.frf','rb')
frf_classic_with_porous_40_2_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_50_1_angle_30_results.frf','rb')
frf_classic_with_porous_50_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_60_1_angle_30_results.frf','rb')
frf_classic_with_porous_60_1_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_60_2_angle_30_results.frf','rb')
frf_classic_with_porous_60_2_angle_30=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_60_results.frf','rb')
frf_classic_with_porous_20_1_angle_60=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_00_results.frf','rb')
frf_classic_with_porous_20_1_angle_00=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_120_results.frf','rb')
frf_classic_with_porous_20_1_angle_120=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_150_results.frf','rb')
frf_classic_with_porous_20_1_angle_150=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_180_results.frf','rb')
frf_classic_with_porous_20_1_angle_180=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_210_results.frf','rb')
frf_classic_with_porous_20_1_angle_210=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_240_results.frf','rb')
frf_classic_with_porous_20_1_angle_240=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_270_results.frf','rb')
frf_classic_with_porous_20_1_angle_270=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_300_results.frf','rb')
frf_classic_with_porous_20_1_angle_300=pickle.load(f)
f.close()

f=open('results/cavity2_with_porous_20_1_angle_330_results.frf','rb')
frf_classic_with_porous_20_1_angle_330=pickle.load(f)
f.close()


#f=open('results/cavity1_with_porous_20_1_angle_60_results.frf','rb')
#frf_classic_with_porous_20_1_angle_60=pickle.load(f)
#f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
#pl.plot(frf_classic_no_porous_20_1_angle_90[0],10*log10(frf_classic_no_porous_20_1_angle_90[1]/prefsquare),'k-',label='Classic, no porous, 20_1, angle 90', linewidth=2)
pl.plot(frf_classic_with_porous_20_1_angle_00[0],10*log10(frf_classic_with_porous_20_1_angle_00[1]/prefsquare),'k-',label='Classic, with porous, 20_1, angle 00', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_30[0],10*log10(frf_classic_with_porous_20_1_angle_30[1]/prefsquare),'r-',label='Classic, with porous, 20_1, angle 30', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_60[0],10*log10(frf_classic_with_porous_20_1_angle_60[1]/prefsquare),'y-',label='Classic, with porous, 20_1, angle 60', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_90[0],10*log10(frf_classic_with_porous_20_1_angle_90[1]/prefsquare),'b-',label='Classic, with porous, 20_1, angle 90', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_120[0],10*log10(frf_classic_with_porous_20_1_angle_120[1]/prefsquare),'m-',label='Classic, with porous, 20_1, angle 120', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_150[0],10*log10(frf_classic_with_porous_20_1_angle_150[1]/prefsquare),'c-',label='Classic, with porous, 20_1, angle 150', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_180[0],10*log10(frf_classic_with_porous_20_1_angle_180[1]/prefsquare),'k--',label='Classic, with porous, 20_1, angle 180', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_210[0],10*log10(frf_classic_with_porous_20_1_angle_210[1]/prefsquare),'r--',label='Classic, with porous, 20_1, angle 210', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_240[0],10*log10(frf_classic_with_porous_20_1_angle_240[1]/prefsquare),'y--',label='Classic, with porous, 20_1, angle 240', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_270[0],10*log10(frf_classic_with_porous_20_1_angle_270[1]/prefsquare),'b--',label='Classic, with porous, 20_1, angle 270', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_300[0],10*log10(frf_classic_with_porous_20_1_angle_300[1]/prefsquare),'m--',label='Classic, with porous, 20_1, angle 300', linewidth=1)
pl.plot(frf_classic_with_porous_20_1_angle_330[0],10*log10(frf_classic_with_porous_20_1_angle_330[1]/prefsquare),'c--',label='Classic, with porous, 20_1, angle 330', linewidth=1)
pl.axis([1.0, 200.0, 40, 160])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
#pl.plot(frf_classic_no_porous_20_1_angle_90[0],10*log10(frf_classic_no_porous_20_1_angle_90[1]/prefsquare),'k-',label='Classic, no porous, 20_1, angle 90', linewidth=2)
#pl.plot(frf_classic_with_porous_20_1_angle_90[0],10*log10(frf_classic_with_porous_20_1_angle_90[1]/prefsquare),'b-',label='Classic, with porous, 20_1, angle 90', linewidth=2)
pl.plot(frf_classic_with_porous_20_1_angle_30[0],10*log10(frf_classic_with_porous_20_1_angle_30[1]/prefsquare),'r-',label='Classic, with porous, 20_1, angle 30', linewidth=2)
#pl.plot(frf_classic_with_porous_20_2_angle_30[0],10*log10(frf_classic_with_porous_20_2_angle_30[1]/prefsquare),'y-',label='Classic, with porous, 20_2, angle 30', linewidth=2)
#pl.plot(frf_classic_with_porous_30_1_angle_30[0],10*log10(frf_classic_with_porous_30_1_angle_30[1]/prefsquare),'k-',label='Classic, with porous, 30_1, angle 30', linewidth=2)
#pl.plot(frf_classic_with_porous_30_2_angle_30[0],10*log10(frf_classic_with_porous_30_2_angle_30[1]/prefsquare),'g-',label='Classic, with porous, 30_2, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_30_3_angle_30[0],10*log10(frf_classic_with_porous_30_3_angle_30[1]/prefsquare),'b-',label='Classic, with porous, 30_3, angle 30', linewidth=2)
#pl.plot(frf_classic_with_porous_40_1_angle_30[0],10*log10(frf_classic_with_porous_40_1_angle_30[1]/prefsquare),'c-',label='Classic, with porous, 40_1, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_40_2_angle_30[0],10*log10(frf_classic_with_porous_40_2_angle_30[1]/prefsquare),'k-',label='Classic, with porous, 40_2, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_50_1_angle_30[0],10*log10(frf_classic_with_porous_50_1_angle_30[1]/prefsquare),'m-',label='Classic, with porous, 50_1, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_60_1_angle_30[0],10*log10(frf_classic_with_porous_60_1_angle_30[1]/prefsquare),'g-',label='Classic, with porous, 60_1, angle 30', linewidth=2)
pl.plot(frf_classic_with_porous_60_2_angle_30[0],10*log10(frf_classic_with_porous_60_2_angle_30[1]/prefsquare),'y-',label='Classic, with porous, 60_2, angle 30', linewidth=2)
#pl.plot(frf_classic_with_porous_20_1_angle_60[0],10*log10(frf_classic_with_porous_20_1_angle_60[1]/prefsquare),'y-',label='Classic, with porous, 20_1, angle 60', linewidth=2)
#pl.plot(frf_classic_with_porous_20_1_angle_00[0],10*log10(frf_classic_with_porous_20_1_angle_00[1]/prefsquare),'k-',label='Classic, with porous, 20_1, angle 00', linewidth=2)
pl.axis([1.0, 200.0, 40, 160])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)


pl.show()

