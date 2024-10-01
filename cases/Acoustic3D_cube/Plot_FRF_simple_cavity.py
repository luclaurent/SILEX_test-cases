from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

#   
##f=open('simple_cavity_h_10/simple_cavity_tet4_results.frf','rb')
##frf_tet4_h_10=pickle.load(f)
##f.close()
##
##f=open('simple_cavity_h_10/simple_cavity_tet10_results.frf','rb')
##frf_tet10_h_10=pickle.load(f)
##f.close()
##

f=open('simple_cavity_h_10/simple_cavity_tet4_results.frf','rb')
frf_tet4_h_10=pickle.load(f)
f.close()

f=open('simple_cavity_h_10/simple_cavity_tet10_results.frf','rb')
frf_tet10_h_10=pickle.load(f)
f.close()

f=open('simple_cavity_h_20/simple_cavity_tet4_results.frf','rb')
frf_tet4_h_20=pickle.load(f)
f.close()

f=open('simple_cavity_h_20/simple_cavity_tet10_results.frf','rb')
frf_tet10_h_20=pickle.load(f)
f.close()

f=open('simple_cavity_h_30/simple_cavity_tet4_results.frf','rb')
frf_tet4_h_30=pickle.load(f)
f.close()
##
####f=open('simple_cavity_h_30/simple_cavity_tet10_results.frf','rb')
####frf_tet10_h_30=pickle.load(f)
####f.close()
##
##f=open('simple_cavity_h_40/simple_cavity_tet4_results.frf','rb')
##frf_tet4_h_40=pickle.load(f)
##f.close()
##
##f=open('simple_cavity_h_50/simple_cavity_tet4_results.frf','rb')
##frf_tet4_h_50=pickle.load(f)
##f.close()
##
##f=open('simple_cavity_h_60/simple_cavity_tet4_results.frf','rb')
##frf_tet4_h_60=pickle.load(f)
##f.close()
##
##f=open('simple_cavity_h_25/simple_cavity_tet10_results.frf','rb')
##frf_tet10_h_25=pickle.load(f)
##f.close()

prefsquare=20e-6**2
pl.figure(1)
##pl.plot(frf_tet4[0],10*log10(frf_tet4[1]/prefsquare),'bo-',label='tet4', linewidth=1)
##pl.plot(frf_tet10[0],10*log10(frf_tet10[1]/prefsquare),'go-',label='tet10', linewidth=1)
pl.plot(frf_tet4_h_10[0],10*log10(frf_tet4_h_10[1]/prefsquare),'b-',label='tet4 / h 10 / 1050 nodes', linewidth=1)
pl.plot(frf_tet10_h_10[0],10*log10(frf_tet10_h_10[1]/prefsquare),'m-',label='tet10 / h 10 / 6850 nodes', linewidth=1)
pl.plot(frf_tet4_h_20[0],10*log10(frf_tet4_h_20[1]/prefsquare),'k-',label='tet4 / h 20 / 5042 nodes', linewidth=1)
pl.plot(frf_tet10_h_20[0],10*log10(frf_tet10_h_20[1]/prefsquare),'r-',label='tet10 / h 20 / 34985 nodes', linewidth=1)
pl.plot(frf_tet4_h_30[0],10*log10(frf_tet4_h_30[1]/prefsquare),'g-',label='tet4 / h 30 / 13280 nodes', linewidth=1)



##pl.figure(1)
###pl.plot(frf_tet4_h_10[0],10*log10(frf_tet4_h_10[1]/prefsquare),'bo-',label='tet4 / h 10 / 1050 dofs', linewidth=1)
###pl.plot(frf_tet10_h_10[0],10*log10(frf_tet10_h_10[1]/prefsquare),'mo-',label='tet10 / h 10 / 6850 dofs', linewidth=1)
##pl.plot(frf_tet4_h_30[0],10*log10(frf_tet4_h_30[1]/prefsquare),'go-',label='tet4 / h 30 / 13280 nodes', linewidth=1)
##pl.plot(frf_tet4_h_40[0],10*log10(frf_tet4_h_40[1]/prefsquare),'yo-',label='tet4 / h 40 / 25210 dofs', linewidth=1)
##pl.plot(frf_tet4_h_50[0],10*log10(frf_tet4_h_50[1]/prefsquare),'ro-',label='tet4 / h 50 / 42478 dofs', linewidth=1)
##pl.plot(frf_tet4_h_60[0],10*log10(frf_tet4_h_60[1]/prefsquare),'bo-',label='tet4 / h 60 / 84102 dofs', linewidth=1)
##pl.plot(frf_tet10_h_25[0],10*log10(frf_tet10_h_25[1]/prefsquare),'c-',label='tet10 / h 25 / 61950 dofs', linewidth=1)
####pl.plot(frf_tet10_h_30[0],10*log10(frf_tet10_h_30[1]/prefsquare),'ro-',label='tet10 / h 30', linewidth=1)
####pl.axis([1.0, 120.0, 70, 105])
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Mean quadratic pressure (dB)')
##pl.grid('on')
##pl.legend(loc=4)
##pl.figure(2)
###pl.plot(frf_tet4_h_10[0],10*log10(frf_tet4_h_10[1]/prefsquare),'bo-',label='tet4 / h 10 / 1050 dofs', linewidth=1)
##pl.plot(frf_tet10_h_10[0],10*log10(frf_tet10_h_10[1]/prefsquare),'m-',label='tet10 / h 10 / 6850 dofs', linewidth=1)
###pl.plot(frf_tet4_h_20[0],10*log10(frf_tet4_h_20[1]/prefsquare),'ko-',label='tet4 / h 20 / 5042 dofs', linewidth=1)
##pl.plot(frf_tet10_h_20[0],10*log10(frf_tet10_h_20[1]/prefsquare),'r-',label='tet10 / h 20 / 34985 dofs', linewidth=1)
###pl.plot(frf_tet4_h_30[0],10*log10(frf_tet4_h_30[1]/prefsquare),'go-',label='tet4 / h 30 / 13280 dofs', linewidth=1)
###pl.plot(frf_tet4_h_40[0],10*log10(frf_tet4_h_40[1]/prefsquare),'yo-',label='tet4 / h 40 / 25210 dofs', linewidth=1)
###pl.plot(frf_tet4_h_50[0],10*log10(frf_tet4_h_50[1]/prefsquare),'ro-',label='tet4 / h 50 / 42478 dofs', linewidth=1)
###pl.plot(frf_tet4_h_60[0],10*log10(frf_tet4_h_60[1]/prefsquare),'bo-',label='tet4 / h 60 / 84102 dofs', linewidth=1)
##pl.plot(frf_tet10_h_25[0],10*log10(frf_tet10_h_25[1]/prefsquare),'c-',label='tet10 / h 25 / 61950 dofs', linewidth=1)
###pl.plot(frf_tet10_h_30[0],10*log10(frf_tet10_h_30[1]/prefsquare),'ro-',label='tet10 / h 30', linewidth=1)
##pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()
