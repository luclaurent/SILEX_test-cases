from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

f=open('results/cavity7_with_porous_air_results.frf','rb')
frf_xfem=pickle.load(f)
f.close()

f=open('../classic/results/cavity7_with_porous_results.frf','rb')
frf_classic=pickle.load(f)
f.close()


## Wrong results because nothing is done for the air-porous interface
f=open('results/cavity7_CB_with_porous_air_results.frf','rb')
frf_CB_xfem=pickle.load(f)
f.close()

f=open('results/cavity7_CB_with_porous_air_v2_results.frf','rb')
frf_CB2_xfem=pickle.load(f)
f.close()

f=open('results/cavity7_CB_with_porous_air_v3_results.frf','rb')
frf_CB3_xfem=pickle.load(f)
f.close()

f=open('results/cavity7_CB_with_porous_air_v4_results.frf','rb')
frf_CB4_xfem=pickle.load(f)
f.close()

f=open('results/cavity7_with_porous_air_flexible_structure_results.frf','rb')
frf_CB5_xfem=pickle.load(f)
f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
#pl.plot(frf_classic[0],10*log10(frf_classic[1]/prefsquare),'k-',label='Classic', linewidth=1)
#pl.plot(frf_xfem[0],10*log10(frf_xfem[1]/prefsquare),'r-',label='xfem', linewidth=1)
#pl.plot(frf_CB_xfem[0],10*log10(frf_CB_xfem[1]/prefsquare),'g-',label='CB + xfem', linewidth=1)
#pl.plot(frf_CB2_xfem[0],10*log10(frf_CB2_xfem[1]/prefsquare),'b-',label='CB full + xfem', linewidth=1)
#pl.plot(frf_CB3_xfem[0],10*log10(frf_CB3_xfem[1]/prefsquare),'g-',label='CB full + xfem + flexible structure', linewidth=1)
pl.plot(frf_CB4_xfem[0],10*log10(frf_CB4_xfem[1]/prefsquare),'m-',label='CB full + xfem + modal flex. struc.', linewidth=1)
pl.plot(frf_CB5_xfem[0],10*log10(frf_CB5_xfem[1]/prefsquare),'y-',label='xfem + flex. struc.', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()

##Number of nodes: 10604
##Number of elements in air: 47915
##Number of elements in porous: 6813
##Number of nodes in air: 9356
##Number of nodes in porous: 2359
##Number of nodes at interface: 1111
##nnodes for structure= 345
##nelem for structure= 624

##Number of nodes: 38631
##Number of elements in air: 176079
##Number of elements in porous: 35655
##Number of nodes in air: 32604
##Number of nodes in porous: 10305
##Number of nodes at interface: 4278
##nnodes for structure= 4243
##nelem for structure= 8242


