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

f=open('results/CavityPorous3_frf_results_angle_study.frf','rb')
[frf,angles]=pickle.load(f)
f.close()
freq_min=850
freq_max=900
weight_function=[]
i=0
for freq in frf[0][0]:
    weight=1e-3
    if freq>freq_min:
        if freq<freq_max:
            i_max=i
            weight=1.0
    else:
        i_min=i+1
    i=i+1
    weight_function.append(weight)
weight_function=scipy.array(weight_function)


prefsquare=20e-6*20e-6
pl.figure(1)
i=0
moyenne=[]
moyenne_weighted=[]
angles_plotted=[]
print(moyenne)
for frf_i in frf:
    pl.plot(frf_i[0],10*log10(frf_i[1]/prefsquare),'-',label='angle '+str(angles[i]), linewidth=1)
    moyenne.append(average(10*log10(frf_i[1]/prefsquare)))
    moyenne_weighted.append(average(10*log10(weight_function[list(range(i_min,i_max,1))]*frf_i[1][list(range(i_min,i_max,1))]/prefsquare)))
    angles_plotted.append(angles[i])
    i=i+1

##f=open('results/CavityPorous3_no_porous_frf_results_angle_study.frf','rb')
##[frf_no_porous,angles_no_porous]=pickle.load(f)
##f.close()
##
##i=0
##moyenne_no_porous=[]
##angles_plotted_no_porous=[]
##print(moyenne)
##for frf_i in frf_no_porous:
##    pl.plot(frf_i[0],10*log10(frf_i[1]/prefsquare),'--',label='no porous, angle '+str(angles_no_porous[i]), linewidth=1)
##    moyenne_no_porous.append(average(10*log10(frf_i[1]/prefsquare)))
##    angles_plotted_no_porous.append(angles_no_porous[i])
##    i=i+1



#pl.axis([10.0, 200.0, 20, 90])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
#pl.plot(angles_plotted,moyenne,'k-',label='with porous')
pl.plot(angles_plotted,moyenne_weighted,'b-',label='with porous, weighted')
#pl.plot(angles_plotted_no_porous,moyenne_no_porous,'b-',label='no porous')
pl.xlabel('angle')
pl.ylabel('Mean of Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

