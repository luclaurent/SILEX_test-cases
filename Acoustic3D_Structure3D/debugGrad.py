### OMP_NUM_THREADS
### OPENBLAS_NUM_THREADS
### MKL_NUM_THREADS


import cProfile
cp = cProfile.Profile()


from Main_acou3D_struc3D_v3_grad import *

freqMin=10.
freqMax=200.
paraVal=scipy.array([2.0,2.0,1.0,1.0])
nbStep=800


#load info from MPI
nbProc,rank,comm=mpiInfo() 

#cp.enable()
dataFRFgrad=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,[0,1,2,3],0)
#cp.disable()
#cp.print_stats()

#finite differences (CD2)
dd=1e-4
ddVX=scipy.array([dd,0.,0.,0.])
paraValBX=paraVal-ddVX
paraValFX=paraVal+ddVX

dataFRFBX=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBX,saveResults=0)
dataFRFFX=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFX,saveResults=0)

ddVY=scipy.array([0.,dd,0.,0.])
paraValBY=paraVal-ddVY
paraValFY=paraVal+ddVY

dataFRFBY=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBY,saveResults=0)
dataFRFFY=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFY,saveResults=0)

ddVZ=scipy.array([0.,0.,dd,0.])
paraValBZ=paraVal-ddVZ
paraValFZ=paraVal+ddVZ

dataFRFBZ=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBZ,saveResults=0)
dataFRFFZ=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFZ,saveResults=0)

ddVR=scipy.array([0.,0.,0.,dd])
paraValBR=paraVal-ddVR
paraValFR=paraVal+ddVR

dataFRFBR=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValBR,saveResults=0)
dataFRFFR=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValFR,saveResults=0)


import pickle

f=open('debug'+'_results.frf','wb')
pickle.dump([dataFRFgrad,dataFRFBX,dataFRFFX,dataFRFBY,dataFRFFY,dataFRFBZ,dataFRFFZ,dataFRFBR,dataFRFFR], f)
#pickle.dump([dataFRFgrad,dataFRFBX,dataFRFFX], f)
f.close()

FF=pickle.load(open('debug_results.frf','rb'))
dataFRFgrad=FF[0]
dataFRFBX=FF[1]
dataFRFFX=FF[2]
dataFRFBY=FF[3]
dataFRFFY=FF[4]
dataFRFBZ=FF[5]
dataFRFFZ=FF[6]
dataFRFBR=FF[7]
dataFRFFR=FF[8]

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

pl.draw()

pl.figure(2)
pl.plot(dataFRFgrad[0],dataFRFgrad[2],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFX[1]-dataFRFBX[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad X')
pl.grid('on')
pl.legend(loc=4)

pl.show()

pl.figure(3)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFBY[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFFY[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.draw()

pl.figure(4)
pl.plot(dataFRFgrad[0],dataFRFgrad[3],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFY[1]-dataFRFBY[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad Y')
pl.grid('on')
pl.legend(loc=4)

pl.show()

pl.figure(5)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFBZ[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFFZ[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.draw()

pl.figure(6)
pl.plot(dataFRFgrad[0],dataFRFgrad[4],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFZ[1]-dataFRFBZ[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad Z')
pl.grid('on')
pl.legend(loc=4)

pl.show()
print("  ")
pl.figure(5)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFBZ[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFFZ[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.draw()

pl.figure(6)
pl.plot(dataFRFgrad[0],dataFRFgrad[4],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFZ[1]-dataFRFBZ[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad Z')
pl.grid('on')
pl.legend(loc=4)

pl.show()
print("  ")

pl.figure(7)
pl.plot(dataFRFgrad[0],10*log10(dataFRFgrad[1]/prefsquare),'r-',label='test 0', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFBR[1]/prefsquare),'b-',label='B', linewidth=1)
pl.plot(dataFRFgrad[0],10*log10(dataFRFFR[1]/prefsquare),'b-',label='F', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.draw()

pl.figure(8)
pl.plot(dataFRFgrad[0],dataFRFgrad[5],'r-',label='grad (full)', linewidth=1)
pl.plot(dataFRFgrad[0],(dataFRFFR[1]-dataFRFBR[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Grad R')
pl.grid('on')
pl.legend(loc=4)

pl.show()
print("  ")

print("  ")

# pl.figure(3)
# pl.plot(dataFRFgrad[0],dataFRFgrad[2],'r-',label='grad (full)', linewidth=1)
# pl.plot(dataFRFgrad[0],(dataFRFFY[1]-dataFRFBY[1])/(2*dd),'b-',label='grad (FD)', linewidth=1)
# #pl.axis([1.0, 120.0, 70, 105])
# pl.xlabel('Frequency (Hz)')
# pl.ylabel('Grad Y')
# pl.grid('on')
# pl.legend(loc=4)

# pl.show()


#cp.print_stats()
