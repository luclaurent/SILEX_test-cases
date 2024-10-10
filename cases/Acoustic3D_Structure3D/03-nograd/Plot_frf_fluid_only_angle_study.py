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

f=open('results/Cavity3_frf_results_angle_study.frf','rb')
[frf,angles]=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6
pl.figure(1)
i=0
moyenne=[]
print(moyenne)
for frf_i in frf:
    pl.plot(frf_i[0],10*log10(frf_i[1]/prefsquare),label='angle '+str(angles[i]), linewidth=1)
    i=i+1
    moyenne.append(average(10*log10(frf_i[1]/prefsquare)))

#pl.axis([10.0, 200.0, 20, 90])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
pl.plot(angles,moyenne)
pl.xlabel('angle')
pl.ylabel('Mean of Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

