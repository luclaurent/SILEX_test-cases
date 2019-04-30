### OMP_NUM_THREADS
### OPENBLAS_NUM_THREADS
### MKL_NUM_THREADS


import cProfile
cp = cProfile.Profile()


from Main_acou3D_struc3D_v3_grad import *

freqMin=10.
freqMax=600.
paraValN=scipy.array([1.0,1.0,0.0,1.0])
nbStep=1600

import scipy

#prepare parameters evolution
nbVal=40
Zmin=0.1
Zmax=2.
Rmin=0.5
Rmax=3.

Zl=scipy.linspace(Zmin,Zmax,nbVal)
Rl=scipy.linspace(Rmin,Rmax,nbVal)
Zm,Rm=scipy.meshgrid(Zl,Rl)
ZZm=list()

#load info from MPI
nbProc,rank,comm=mpiInfo() 

#along values
it=0.
for key,val in scipy.ndenumerate(Rm):
    paraValRun=paraValN.copy()
    paraValRun[2]=Zm[key]
    paraValRun[3]=Rm[key]
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(it)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    Zval=RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraValRun,saveResults=0)
    if rank==0:
        ZZm.append(Zval.copy())
    it=it+1

if rank ==0:
    import pickle
    
    f=open('paraZR'+str(nbStep)+'-'+str(nbVal**2)+'.pck','wb')
    pickle.dump([Zm,Rm,ZZm,[freqMin,freqMax,nbStep,paraValN,nbVal]], f)
    f.close()
