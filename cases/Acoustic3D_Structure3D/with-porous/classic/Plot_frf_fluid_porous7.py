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

f=open('results/cavity7_with_porous_results.frf','rb')
frf=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf[0],10*log10(frf[1]/prefsquare),'k-',label='Classic, with porous', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

##i_min=161
##i_max=199
##f_min=frf_angle_0_20_2[0][i_min]
##f_max=frf_angle_0_20_2[0][i_max]
##print("freq min=",f_min)
##print("freq max=",f_max)
##moyenne=[average(10*log10(frf_angle_0_20_2[1][i_min:i_max]/prefsquare)),
##         average(10*log10(frf_angle_30_20_2[1][i_min:i_max]/prefsquare)),
##         average(10*log10(frf_angle_60_20_2[1][i_min:i_max]/prefsquare)),
##         average(10*log10(frf_angle_90_20_2[1][i_min:i_max]/prefsquare)),
##         average(10*log10(frf_angle_120_20_2[1][i_min:i_max]/prefsquare)),
##         average(10*log10(frf_angle_150_20_2[1][i_min:i_max]/prefsquare))]
##
##angles=[0,30,60,90,120,150]
##
##pl.figure(2)
##pl.plot(angles,moyenne,'bo-',label='with porous')
##pl.xlabel('angle')
##pl.ylabel('Mean of Mean quadratic pressure (dB)')
##pl.grid('on')
##pl.legend(loc=4)


pl.show()

