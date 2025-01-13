import string
import time
import os
import pylab as pl
import pickle
import scipy

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/classic-test1_pressTip.frf','r')
AllpressTipRef=pickle.load(f)
f.close()
pressTipRef=AllpressTipRef[0]
nodes=AllpressTipRef[1]
freq_ref=AllpressTipRef[2]
nnodes=len(pressTipRef)

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test1-with-edge_pressTip.frf','r')
AllpressTipXfem1=pickle.load(f)
f.close()
pressTipXfem1=AllpressTipXfem1[0]
thetaXfem1=AllpressTipXfem1[1]
freq_Xfem1=AllpressTipXfem1[2]

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test2-with-edge_pressTip.frf','r')
AllpressTipXfem2=pickle.load(f)
f.close()
pressTipXfem2=AllpressTipXfem2[0]
thetaXfem2=AllpressTipXfem2[1]
freq_Xfem2=AllpressTipXfem2[2]

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test3-with-edge_pressTip.frf','r')
AllpressTipXfem3=pickle.load(f)
f.close()
pressTipXfem3=AllpressTipXfem3[0]
thetaXfem3=AllpressTipXfem3[1]
freq_Xfem3=AllpressTipXfem3[2]

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test4-with-edge_pressTip.frf','r')
AllpressTipXfem4=pickle.load(f)
f.close()
pressTipXfem4=AllpressTipXfem4[0]
thetaXfem4=AllpressTipXfem4[1]
freq_Xfem4=AllpressTipXfem4[2]

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test5-with-edge_pressTip.frf','r')
AllpressTipXfem5=pickle.load(f)
f.close()
pressTipXfem5=AllpressTipXfem5[0]
thetaXfem5=AllpressTipXfem5[1]
freq_Xfem5=AllpressTipXfem5[2]

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test6-with-edge_pressTip.frf','r')
AllpressTipXfem6=pickle.load(f)
f.close()
pressTipXfem6=AllpressTipXfem6[0]
thetaXfem6=AllpressTipXfem6[1]
freq_Xfem6=AllpressTipXfem6[2]

freq_Xfem=[freq_Xfem1,freq_Xfem2,freq_Xfem3,freq_Xfem4,freq_Xfem5,freq_Xfem6]

print 'Frequence at which it has been computed: reference =',freq_ref,'Hz   ////       Xfem=',freq_Xfem,'Hz'

XC=0.6
YC=0.65

thetaRef=[]
pref=20e-6

for i in range(nnodes):
    x=nodes[i][0]-XC
    y=nodes[i][1]-YC
    thetaRef.append(scipy.arctan(x/y)*180.0/scipy.pi)

pl.figure(1)
pl.plot(thetaRef,20*scipy.log10(scipy.real(abs(pressTipRef))/pref),'ko-',label='Reference', linewidth=2)
pl.plot(thetaXfem1,20*scipy.log10(scipy.real(abs(pressTipXfem1))/pref),'y-', linewidth=2)
pl.plot(thetaXfem1[range(0,len(thetaXfem1),10)],20*scipy.log10(scipy.real(abs(pressTipXfem1[range(0,len(thetaXfem1),10)]))/pref),'*y-',label='Mesh 1')
pl.plot(thetaXfem2,20*scipy.log10(scipy.real(abs(pressTipXfem2))/pref),'m-', linewidth=2)
pl.plot(thetaXfem2[range(0,len(thetaXfem2),10)],20*scipy.log10(scipy.real(abs(pressTipXfem2[range(0,len(thetaXfem2),10)]))/pref),'xm-',label='Mesh 2')
pl.plot(thetaXfem3,20*scipy.log10(scipy.real(abs(pressTipXfem3))/pref),'g-', linewidth=2)
pl.plot(thetaXfem3[range(0,len(thetaXfem3),10)],20*scipy.log10(scipy.real(abs(pressTipXfem3[range(0,len(thetaXfem3),10)]))/pref),'>g-',label='Mesh 3')
pl.plot(thetaXfem4,20*scipy.log10(scipy.real(abs(pressTipXfem4)/pref)),'b-', linewidth=2)
pl.plot(thetaXfem4[range(0,len(thetaXfem4),10)],20*scipy.log10(scipy.real(abs(pressTipXfem4[range(0,len(thetaXfem4),10)]))/pref),'<b-',label='Mesh 4')
pl.plot(thetaXfem5,20*scipy.log10(scipy.real(abs(pressTipXfem5))/pref),'c-', linewidth=2)
pl.plot(thetaXfem5[range(0,len(thetaXfem5),10)],20*scipy.log10(scipy.real(abs(pressTipXfem5[range(0,len(thetaXfem5),10)]))/pref),'^c-',label='Mesh 5')
pl.plot(thetaXfem5,20*scipy.log10(scipy.real(abs(pressTipXfem6))/pref),'r-', linewidth=2)
pl.plot(thetaXfem6[range(0,len(thetaXfem6),10)],20*scipy.log10(scipy.real(abs(pressTipXfem6[range(0,len(thetaXfem6),10)]))/pref),'sr-',label='Mesh 6')

pl.axis([-90.0, 90.0, 70, 110])

pl.xlabel('Angle (degree)')
pl.ylabel('Pressure (dB)')
pl.grid('on')
pl.legend()
pl.show()
