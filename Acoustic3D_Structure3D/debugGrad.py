from Main_acou3D_struc3D_v3_grad import *

freqMin=10.
freqMax=100.
paraVal=scipy.array([1.0,1.0])
nbStep=400


#load info from MPI
nbProc,rank,comm=mpiInfo()  

dataFRFgrad=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal)

#finite differences (CD2)
dd=1e-1
ddVX=scipy.array([dd,0.])
paraValBX=paraVal-ddVX
paraValFX=paraVal+ddVX

dataFRFBX=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBX)
dataFRFFX=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFX)

ddVY=scipy.array([0.,dd])
paraValBY=paraVal-ddVY
paraValFY=paraVal+ddVY

dataFRFBY=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBY)
dataFRFFY=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFY)


f=open('debug'+'_results.frf','wb')
pickle.dump([dataFRFgrad,dataFRFBX,dataFRFFX,dataFRFBY,dataFRFFY], f)
f.close()

FF=pickle.load(open('debug_results.frf','rb'))
dataFRFgrad=FF[0]
dataFRFBX=FF[1]
dataFRFFX=FF[2]
dataFRFBY=FF[3]
dataFRFFY=FF[4]

import pylab as pl
from numpy import *

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFBX[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFFX[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

pl.figure(2)
pl.plot(dataFRFgrad[0],dataFRFgrad[2],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFX[1]-dataFRFBX[1])/2*dd,'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad X')
pl.grid('on')
pl.legend(loc=4)

pl.show()

pl.figure(3)
pl.plot(dataFRFgrad[0],dataFRFgrad[2],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFY[1]-dataFRFBY[1])/2*dd,'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad Y')
pl.grid('on')
pl.legend(loc=4)

pl.show()

