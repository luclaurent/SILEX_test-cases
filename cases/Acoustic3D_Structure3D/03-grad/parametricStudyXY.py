### OMP_NUM_THREADS
### OPENBLAS_NUM_THREADS
### MKL_NUM_THREADS


import cProfile
cp = cProfile.Profile()


from Main_acou3D_struc3D_v3_grad import *

freqMin=10.
freqMax=200.
paraValN=scipy.array([1.0,1.0,1.,1.0])
nbStep=800

import scipy

#prepare parameters evolution
nbVal=30
Xmin=0.
Xmax=4.
Ymin=0.
Ymax=3.

Xl=scipy.linspace(Xmin,Xmax,nbVal)
Yl=scipy.linspace(Ymin,Ymax,nbVal)
Xm,Ym=scipy.meshgrid(Xl,Yl)
Zm=list()

#load info from MPI
nbProc,rank,comm=mpiInfo() 

#along values
it=0.
for key,val in scipy.ndenumerate(Xm):
    paraValRun=paraValN.copy()
    paraValRun[0]=Xm[key]
    paraValRun[1]=Ym[key]
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
        Zm.append(Zval.copy())
    it=it+1

if rank ==0:
    import pickle
    
    f=open('paraXY_'+str(freqMin)+'-'+str(freqMax)+'_'+str(nbStep)+'-'+str(nbVal**2)+'.pck','wb')
    pickle.dump([Xm,Ym,Zm,[freqMin,freqMax,nbStep,paraValN,nbVal]], f)
    f.close()

    import scipy.io
    scipy.io.savemat('paraXY_'+str(freqMin)+'-'+str(freqMax)+'_'+str(nbStep)+'-'+str(nbVal**2)+'.mat',mdict={'Xm':Xm,'Ym':Ym,'FRF':FRF,'freqMin':freqMin,'freqMax':freqMax,'nbStep':nbStep,'paraValN':paraValN,'nbVal':nbVal})


