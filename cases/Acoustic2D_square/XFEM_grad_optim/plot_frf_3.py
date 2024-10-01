from numpy import *
import string
import time
import os
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
f=open('xfem_3_6199_results.frf','rb')
frf_1=pickle.load(f)
f.close()

f=open('xfem_3_6200_results.frf','rb')
frf_2=pickle.load(f)
f.close()

f=open('xfem_3_6201_results.frf','rb')
frf_3=pickle.load(f)
f.close()

#id_node=7
#id_freq=20
#print('freq=',frf_6200[0][id_freq])
#press_node_6200=frf_6200[3][id_freq][id_node]
#press_node_6202=frf_6202[3][id_freq][id_node]
#dpress_node_diff=(press_node_6202-press_node_6200)/(0.6202-0.6200)
#dpress_node_gradient=frf_6201[4][id_freq][id_node]

prefsquare=20e-6*20e-6

Dp_Dtheta_diff = (frf_3[1]-frf_1[1])/(0.6202-0.6200)
Dp_Dtheta_grad = frf_2[2]

pl.figure(1)
pl.plot(frf_1[0],10*log10(frf_1[1]/prefsquare),'k-',label='0.6199', linewidth=1)
pl.plot(frf_2[0],10*log10(frf_2[1]/prefsquare),'g-',label='0.6200', linewidth=1)
pl.plot(frf_3[0],10*log10(frf_3[1]/prefsquare),'b-',label='0.6201', linewidth=1)
#pl.axis([100.0, 300.0, 70, 130])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=3)

pl.figure(2)
pl.plot(frf_2[0],Dp_Dtheta_diff,'ok-',label='finite diff.', linewidth=1)
pl.plot(frf_2[0],Dp_Dtheta_grad,'og-',label='analytic', linewidth=1)

pl.figure(4)
pl.plot(frf_2[0],Dp_Dtheta_diff/Dp_Dtheta_grad,'ok',label='finite diff./analytic', linewidth=1)
pl.axis([100.0, 500.0, 0.95, 1.05])

pl.figure(3)
for id_freq in range(len(frf_1[0])):
    id_node=5
    #print('freq=',frf_6200[0][id_freq])
    press_node_1=frf_1[3][id_freq][id_node-1]
    press_node_3=frf_3[3][id_freq][id_node-1]
    #print(press_node_1,press_node_3)
    dpress_node_diff=(press_node_3-press_node_1)/(0.6201-0.6199)
    dpress_node_gradient=frf_2[4][id_freq][id_node-1]
    pl.plot(frf_2[0][id_freq],dpress_node_diff,'ok')
    pl.plot(frf_2[0][id_freq],dpress_node_gradient,'og')
#    pl.axis([100.0, 500.0, 0.9, 1.1])

pl.figure(5)
for id_freq in range(len(frf_1[0])):
    id_node=5
    #print('freq=',frf_6200[0][id_freq])
    press_node_1=frf_1[3][id_freq][id_node-1]
    press_node_3=frf_3[3][id_freq][id_node-1]
    dpress_node_diff=(press_node_3-press_node_1)/(0.6201-0.6199)
    dpress_node_gradient=frf_2[4][id_freq][id_node-1]
    pl.plot(frf_2[0][id_freq],dpress_node_diff/dpress_node_gradient,'ok')
    pl.axis([100.0, 500.0, 0.95, 1.05])
    #pl.plot(frf_2[0][id_freq],dpress_node_gradient,'og')

pl.show()

