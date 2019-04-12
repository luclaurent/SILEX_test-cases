from Main_acou3D_struc3D_v3_grad import *

freqMin=10.
freqMax=200.
paraVal=scipy.array([1.0,1.0])
nbStep=200


#load info from MPI
nbProc,rank,comm=mpiInfo()  

dataFRFgrad=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal)

#finite differences (CD2)
dd=1e-6
ddV=scipy.array([dd,0.])
paraValB=paraVal-ddV
paraValF=paraVal+ddV

dataFRFB=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValB)
dataFRFF=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValF)


f=open('debug'+'_results.frf','wb')
pickle.dump([dataFRFgrad,dataFRFB,dataFRFF], f)
f.close()

FF=pickle.load(open('debug_results.frf','rb'))
dataFRFgrad=FF[0]
dataFRFB=FF[1]
dataFRFF=FF[2]

import pylab as pl
from numpy import *

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFB[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFF[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

pl.figure(2)
pl.plot(dataFRFgrad[0],dataFRFgrad[2],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFF[1]-dataFRFB[1])/2*dd,'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

