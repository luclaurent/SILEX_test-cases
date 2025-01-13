from numpy import *
import string
import time
import os
import matplotlib
matplotlib.use('TkAgg')
import pylab as pl
import pickle

##f=open('xfem_3_62000_results.frf','rb')
##frf_1=pickle.load(f)
##f.close()
##
##f=open('xfem_3_62001_results.frf','rb')
##frf_2=pickle.load(f)
##f.close()
##
##f=open('xfem_3_62002_results.frf','rb')
##frf_3=pickle.load(f)
##f.close()
ii=0
val=[0]*7
valP=[0]*7
frf=[0]*7
for itP in [0,1,5,10,50,100,500]:
    val[ii]=620000+itP
    valP[ii]=(val[ii]-620000)/1000000
    filename='xfem_3_'+str(val[ii])+'_results.frf'
    f=open(filename,'rb')
    frf[ii]=pickle.load(f)
    f.close()
    ii=ii+1
    #frf_1=pickle.load(f)
    #f.close()

    #f=open('xfem_3_6201_results.frf','rb')
    #frf_2=pickle.load(f)
    #f.close()

    #f=open('xfem_3_6202_results.frf','rb')
    #frf_3=pickle.load(f)
    #f.close()

#id_node=7
#id_freq=20
#print('freq=',frf_6200[0][id_freq])
#press_node_6200=frf_6200[3][id_freq][id_node]
#press_node_6202=frf_6202[3][id_freq][id_node]
#dpress_node_diff=(press_node_6202-press_node_6200)/(0.6202-0.6200)
#dpress_node_gradient=frf_6201[4][id_freq][id_node]

prefsquare=20e-6*20e-6
Dptmp=[0]*6
pl.figure(1)
for it in range(1,6):
    Dptmp[it]=(frf[it][1]-frf[0][1])/valP[1]
    pl.plot(frf[0][0],Dptmp[it],'k-', linewidth=1)


# Dp_Dtheta_diff = (frf_3[1]-frf_1[1])/(0.6202-0.6200)
# Dp_Dtheta_grad = frf_2[2]

# pl.figure(1)
# pl.plot(frf_1[0],10*log10(frf_1[1]/prefsquare),'k-',label='0.62000', linewidth=1)
# pl.plot(frf_2[0],10*log10(frf_2[1]/prefsquare),'g-',label='0.62001', linewidth=1)
# pl.plot(frf_3[0],10*log10(frf_3[1]/prefsquare),'b-',label='0.62002', linewidth=1)
# #pl.axis([100.0, 300.0, 70, 130])
# pl.xlabel('Frequency (Hz)')
# pl.ylabel('Mean quadratic pressure (dB)')
# pl.grid('on')
# pl.legend(loc=3)

# pl.figure(2)
# pl.plot(frf_2[0],Dp_Dtheta_diff,'ok-',label='finite diff.', linewidth=1)
# pl.plot(frf_2[0],Dp_Dtheta_grad,'og-',label='analytic', linewidth=1)

# pl.figure(4)
# pl.plot(frf_2[0],Dp_Dtheta_diff/Dp_Dtheta_grad,'ok-',label='finite diff./analytic', linewidth=1)


# pl.figure(3)
# for id_freq in range(len(frf_1[0])):
#     id_node=5
#     #print('freq=',frf_6200[0][id_freq])
#     press_node_1=frf_1[3][id_freq][id_node-1]
#     press_node_3=frf_3[3][id_freq][id_node-1]
#     dpress_node_diff=(press_node_3-press_node_1)/(0.6202-0.6200)
#     dpress_node_gradient=frf_2[4][id_freq][id_node]
#     pl.plot(frf_2[0][id_freq],dpress_node_diff,'ok')
#     pl.plot(frf_2[0][id_freq],dpress_node_gradient,'og')
# #pl.axis([100.0, 200.0, -1, 10])

# pl.figure(5)
# for id_freq in range(len(frf_1[0])):
#     id_node=5
#     #print('freq=',frf_6200[0][id_freq])
#     press_node_1=frf_1[3][id_freq][id_node-1]
#     press_node_3=frf_3[3][id_freq][id_node-1]
#     dpress_node_diff=(press_node_3-press_node_1)/(0.6202-0.6200)
#     dpress_node_gradient=frf_2[4][id_freq][id_node]
#     pl.plot(frf_2[0][id_freq],dpress_node_diff/dpress_node_gradient,'ok')
#     #pl.plot(frf_2[0][id_freq],dpress_node_gradient,'og')

# pl.show()

