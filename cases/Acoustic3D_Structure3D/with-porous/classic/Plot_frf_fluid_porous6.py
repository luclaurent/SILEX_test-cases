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

f=open('results/cavity6_with_porous_angle_0_20_2_results.frf','rb')
frf_angle_0_20_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_30_20_2_results.frf','rb')
frf_angle_30_20_2=pickle.load(f)
f.close()


f=open('results/cavity6_with_porous_angle_60_15_2_results.frf','rb')
frf_angle_60_15_2=pickle.load(f)
f.close()

#f=open('results/cavity6_with_porous_angle_60_15_2bis_results.frf','rb')
#frf_angle_60_15_2bis=pickle.load(f)
#f.close()

f=open('results/cavity6_with_porous_angle_60_10_2_results.frf','rb')
frf_angle_60_10_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_60_20_1_results.frf','rb')
frf_angle_60_20_1=pickle.load(f)
f.close()



f=open('results/cavity6_with_porous_angle_60_20_3_results.frf','rb')
frf_angle_60_20_3=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_60_20_2_results.frf','rb')
frf_angle_60_20_2=pickle.load(f)
f.close()


f=open('results/cavity6_with_porous_angle_60_30_2_results.frf','rb')
frf_angle_60_30_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_60_40_4_results.frf','rb')
frf_angle_60_40_4=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_90_20_2_results.frf','rb')
frf_angle_90_20_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_120_20_2_results.frf','rb')
frf_angle_120_20_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_150_20_2_results.frf','rb')
frf_angle_150_20_2=pickle.load(f)
f.close()

f=open('results/cavity6_with_porous_angle_150_30_3_results.frf','rb')
frf_angle_150_30_3=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
#pl.plot(frf_angle_60_15_2[0],10*log10(frf_angle_60_15_2[1]/prefsquare),'k-',label='Classic, with porous, angle 60, 15_2', linewidth=1)
#pl.plot(frf_angle_60_15_2bis[0],10*log10(frf_angle_60_15_2bis[1]/prefsquare),'r-',label='Classic, with porous, angle 60, 15_2bis', linewidth=1)
#pl.plot(frf_angle_60_10_2[0],10*log10(frf_angle_60_10_2[1]/prefsquare),'b-',label='Classic, with porous, angle 60, 10_2', linewidth=1)
#pl.plot(frf_angle_60_20_3[0],10*log10(frf_angle_60_20_3[1]/prefsquare),'r-',label='Classic, with porous, angle 60, 20_3', linewidth=1)
pl.plot(frf_angle_0_20_2[0],10*log10(frf_angle_0_20_2[1]/prefsquare),'k-',label='Classic, with porous, angle 0, 20_2', linewidth=1)
pl.plot(frf_angle_30_20_2[0],10*log10(frf_angle_30_20_2[1]/prefsquare),'r-',label='Classic, with porous, angle 30, 20_2', linewidth=1)
pl.plot(frf_angle_60_20_2[0],10*log10(frf_angle_60_20_2[1]/prefsquare),'g-',label='Classic, with porous, angle 60, 20_2', linewidth=1)
pl.plot(frf_angle_60_20_1[0],10*log10(frf_angle_60_20_1[1]/prefsquare),'y-',label='Classic, with porous, angle 60, 20_1', linewidth=1)
pl.plot(frf_angle_90_20_2[0],10*log10(frf_angle_90_20_2[1]/prefsquare),'m-',label='Classic, with porous, angle 90, 20_2', linewidth=1)
#pl.plot(frf_angle_60_30_2[0],10*log10(frf_angle_60_30_2[1]/prefsquare),'m-',label='Classic, with porous, angle 60, 30_2', linewidth=1)
#pl.plot(frf_angle_60_40_4[0],10*log10(frf_angle_60_40_4[1]/prefsquare),'c-',label='Classic, with porous, angle 60, 40_4', linewidth=1)
pl.plot(frf_angle_120_20_2[0],10*log10(frf_angle_120_20_2[1]/prefsquare),'c-',label='Classic, with porous, angle 120, 20_2', linewidth=1)
pl.plot(frf_angle_150_20_2[0],10*log10(frf_angle_150_20_2[1]/prefsquare),'b-',label='Classic, with porous, angle 150, 20_2', linewidth=1)
#pl.plot(frf_angle_150_30_3[0],10*log10(frf_angle_150_30_3[1]/prefsquare),'k-',label='Classic, with porous, angle 150, 30_3', linewidth=1)
pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

i_min=161
i_max=199
f_min=frf_angle_0_20_2[0][i_min]
f_max=frf_angle_0_20_2[0][i_max]
print("freq min=",f_min)
print("freq max=",f_max)
moyenne=[average(10*log10(frf_angle_0_20_2[1][i_min:i_max]/prefsquare)),
         average(10*log10(frf_angle_30_20_2[1][i_min:i_max]/prefsquare)),
         average(10*log10(frf_angle_60_20_2[1][i_min:i_max]/prefsquare)),
         average(10*log10(frf_angle_90_20_2[1][i_min:i_max]/prefsquare)),
         average(10*log10(frf_angle_120_20_2[1][i_min:i_max]/prefsquare)),
         average(10*log10(frf_angle_150_20_2[1][i_min:i_max]/prefsquare))]

angles=[0,30,60,90,120,150]

pl.figure(2)
pl.plot(angles,moyenne,'bo-',label='with porous')
pl.xlabel('angle')
pl.ylabel('Mean of Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)


pl.show()

