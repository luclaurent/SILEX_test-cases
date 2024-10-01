import numpy,pickle,scipy.io


pp=pickle.load(open('paraZR800-800.pck','rb'))

#pickle.dump([Xm,Ym,Zm,[freqMin,freqMax,nbStep,paraValN,nbVal]], f)
Zm=pp[0]
Rm=pp[1]
FRF=pp[2]
freqMin=pp[3][0]
freqMax=pp[3][1]
nbStep=pp[3][2]
paraValN=pp[3][3]
nbVal=pp[3][4]

scipy.io.savemat('paraZR800-800.mat',mdict={'Zm':Zm,'Rm':Rm,'FRF':FRF,'freqMin':freqMin,'freqMax':freqMax,'nbStep':nbStep,'paraValN':paraValN,'nbVal':nbVal})