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


f=open('results/cavity4_with_porous_dy_1_0_results.frf','rb')
frf_classic_with_porous_dy_1_0=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_1_5_results.frf','rb')
frf_classic_with_porous_dy_1_5=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_2_0_results.frf','rb')
frf_classic_with_porous_dy_2_0=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_2_5_results.frf','rb')
frf_classic_with_porous_dy_2_5=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_3_0_results.frf','rb')
frf_classic_with_porous_dy_3_0=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_3_5_results.frf','rb')
frf_classic_with_porous_dy_3_5=pickle.load(f)
f.close()
f=open('results/cavity4_with_porous_dy_4_0_results.frf','rb')
frf_classic_with_porous_dy_4_0=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle0_results.frf','rb')
frf_classic_with_porous_angle_0=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle0_30_2_results.frf','rb')
frf_classic_with_porous_angle_0_30_2=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle0_25_3_results.frf','rb')
frf_classic_with_porous_angle_0_25_3=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle0_30_3_results.frf','rb')
frf_classic_with_porous_angle_0_30_3=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle0_25_4_results.frf','rb')
frf_classic_with_porous_angle_0_25_4=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle90_results.frf','rb')
frf_classic_with_porous_angle_90=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle180_results.frf','rb')
frf_classic_with_porous_angle_180=pickle.load(f)
f.close()

f=open('results/cavity4_with_porous_angle270_results.frf','rb')
frf_classic_with_porous_angle_270=pickle.load(f)
f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
#pl.plot(frf_classic_with_porous_dy_1_0[0],10*log10(frf_classic_with_porous_dy_1_0[1]/prefsquare),'k-',label='Classic, with porous, dy=1.0', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_1_5[0],10*log10(frf_classic_with_porous_dy_1_5[1]/prefsquare),'r-',label='Classic, with porous, dy=1.5', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_2_0[0],10*log10(frf_classic_with_porous_dy_2_0[1]/prefsquare),'g-',label='Classic, with porous, dy=2.0', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_2_5[0],10*log10(frf_classic_with_porous_dy_2_5[1]/prefsquare),'b-',label='Classic, with porous, dy=2.5', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_3_0[0],10*log10(frf_classic_with_porous_dy_3_0[1]/prefsquare),'c-',label='Classic, with porous, dy=3.0', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_3_5[0],10*log10(frf_classic_with_porous_dy_3_5[1]/prefsquare),'y-',label='Classic, with porous, dy=3.5', linewidth=1)
#pl.plot(frf_classic_with_porous_dy_4_0[0],10*log10(frf_classic_with_porous_dy_4_0[1]/prefsquare),'m-',label='Classic, with porous, dy=4.0', linewidth=1)
pl.plot(frf_classic_with_porous_angle_0[0],10*log10(frf_classic_with_porous_angle_0[1]/prefsquare),'k-',label='Classic, with porous, angle 0, 20_2', linewidth=1)
pl.plot(frf_classic_with_porous_angle_0_30_2[0],10*log10(frf_classic_with_porous_angle_0_30_2[1]/prefsquare),'m-',label='Classic, with porous, angle 0, 30_2', linewidth=1)
pl.plot(frf_classic_with_porous_angle_0_25_3[0],10*log10(frf_classic_with_porous_angle_0_25_3[1]/prefsquare),'y-',label='Classic, with porous, angle 0, 25_3', linewidth=1)
pl.plot(frf_classic_with_porous_angle_0_30_3[0],10*log10(frf_classic_with_porous_angle_0_30_3[1]/prefsquare),'g-',label='Classic, with porous, angle 0, 30_3', linewidth=1)
pl.plot(frf_classic_with_porous_angle_0_25_4[0],10*log10(frf_classic_with_porous_angle_0_25_4[1]/prefsquare),'b-',label='Classic, with porous, angle 0, 25_4', linewidth=1)
#pl.plot(frf_classic_with_porous_angle_90[0],10*log10(frf_classic_with_porous_angle_90[1]/prefsquare),'r-',label='Classic, with porous, angle 90', linewidth=1)
#pl.plot(frf_classic_with_porous_angle_180[0],10*log10(frf_classic_with_porous_angle_180[1]/prefsquare),'g-',label='Classic, with porous, angle 180', linewidth=1)
#pl.plot(frf_classic_with_porous_angle_270[0],10*log10(frf_classic_with_porous_angle_270[1]/prefsquare),'b-',label='Classic, with porous, angle 270', linewidth=1)
#pl.axis([1.0, 200.0, 40, 160])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)



pl.show()

