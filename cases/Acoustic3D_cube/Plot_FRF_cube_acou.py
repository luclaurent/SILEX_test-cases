from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path

filename = Path(__file__).parent / 'results/cube_acou_tet4_results.frf'
f=open(filename,'rb')
frf_tet4=pickle.load(f)
f.close()

#
#filename = Path(__file__).parent / 'h_30/cube_acou_tet4_results.frf'
#f=open(filename,'rb')
#frf_tet4_h_30=pickle.load(f)
#f.close()

#f=open('h_20/cube_acou_tet4_results.frf','rb')
#frf_tet4_h_20=pickle.load(f)
#f.close()
#
#f=open('h_20/cube_acou_tet10_results.frf','rb')
#frf_tet10_h_20=pickle.load(f)
#f.close()
#
#f=open('h_10/cube_acou_tet4_results.frf','rb')
#frf_tet4_h_10=pickle.load(f)
#f.close()
#
#f=open('h_10/cube_acou_tet10_results.frf','rb')
#frf_tet10_h_10=pickle.load(f)
#f.close()
#
#f=open('cube_acou_tet4_results.frf','rb')
#frf_tet4=pickle.load(f)
#f.close()

prefsquare=20e-6**2

pl.figure(1)
#pl.plot(frf_tet4_h_30[0],10*log10(frf_tet4_h_30[1]/prefsquare),'g-',label='tet4 / h 30', linewidth=1)
#pl.plot(frf_tet4_h_20[0],10*log10(frf_tet4_h_20[1]/prefsquare),'b-',label='tet4 / h 20', linewidth=1)
pl.plot(frf_tet4[0],10*log10(frf_tet4[1]/prefsquare),'m-',label='tet4', linewidth=1)

#pl.plot(frf_tet10_h_20[0],10*log10(frf_tet10_h_20[1]/prefsquare),'ro-',label='tet10 / h 20', linewidth=1)
#pl.plot(frf_tet4_h_10[0],10*log10(frf_tet4_h_10[1]/prefsquare),'mo-',label='tet4 / h 10', linewidth=1)
#pl.plot(frf_tet10_h_10[0],10*log10(frf_tet10_h_10[1]/prefsquare),'go-',label='tet10 / h 10', linewidth=1)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()

#pl.figure(2)
#pl.plot(frf_tet4[0],10*log10(frf_tet4[1])-10*log(frf_tet10[1]),'bo-',label='tet4/tet10', linewidth=2)
##pl.axis([1.0, 120.0, 70, 105])
#pl.xlabel('Frequency (Hz)')
#pl.ylabel('Mean quadratic pressure (dB)')
#pl.grid('on')
#pl.legend(loc=4)

