from numpy import *
import string
import time
import scipy
#from scipy.sparse import *
#from scipy import *
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen
#from numpy.linalg import solve, norm
import os
#from scipy.linalg import decomp
import pylab as pl
import pickle

f=open('results/classic_results_no_damping.frf','rb')
frf_classic_no_damping=pickle.load(f)
f.close()

f=open('results/classic_results_with_damping.frf','rb')
frf_classic_with_damping=pickle.load(f)
f.close()

f=open('results/classic_pgd.frf','rb')
frf_classic_pgd=pickle.load(f)
f.close()
frf_classic_pgd[1]=scipy.array(frf_classic_pgd[1])
prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_classic_no_damping[0],10*log10(frf_classic_no_damping[1]/prefsquare),'k-',label='Classic, no damping', linewidth=2)
pl.plot(frf_classic_pgd[0],10*log10(frf_classic_pgd[1]/prefsquare),'og-',label='Classic, PGD', linewidth=2)
#pl.plot(frf_classic_with_damping[0],10*log10(frf_classic_with_damping[1]/prefsquare),'b-',label='Classic, with damping', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)
pl.show()


##f=open('/home/legay/Codes/XFEM-Acoustique/V2013/Sieges/results/xfem_struc_projection_fluid_damping_results.frf','r')
##frf_xfem_reference=pickle.load(f)
##f.close()
##
##f=open('/home/legay/Codes/XFEM-Acoustique/V2013/Sieges/results/xfem_struc_CB_projection_fluid_damping_250F-44S_results.frf','r')
##frf_CB_xfem_1=pickle.load(f)
##f.close()
##
##f=open('/home/legay/Codes/XFEM-Acoustique/V2013/Sieges/results/xfem_struc_CB_projection_fluid_damping_102F-32S_results.frf','r')
##frf_CB_xfem_2=pickle.load(f)
##f.close()
##
##f=open('/home/legay/Codes/XFEM-Acoustique/V2013/Sieges/results/xfem_struc_CB_projection_fluid_damping_38F-24S_results.frf','r')
##frf_CB_xfem_3=pickle.load(f)
##f.close()


##prefsquare=20e-6*20e-6
##
##dBerror_1 =abs(10*log10(frf_CB_xfem_1[1]/prefsquare)-10*log10(frf_xfem_reference[1]/prefsquare))
##dBerror_2 =abs(10*log10(frf_CB_xfem_2[1]/prefsquare)-10*log10(frf_xfem_reference[1]/prefsquare))
##dBerror_3 =abs(10*log10(frf_CB_xfem_3[1]/prefsquare)-10*log10(frf_xfem_reference[1]/prefsquare))
##
##globalerror=[(sum(dBerror_1))/400,(sum(dBerror_2))/400,(sum(dBerror_3))/400]
##
##lastfreqF=[ 437.33080439115651,300.8908154499278,200.63124208219497]
##lastfreqS=[ 468.71121478270476,298.11692698746924,202.25256304638719]
##
##pl.figure(1)
##pl.plot(frf_xfem_reference[0],10*log10(frf_xfem_reference[1]/prefsquare),'k-',label='Reference', linewidth=2)
##pl.plot(frf_CB_xfem_3[0],10*log10(frf_CB_xfem_3[1]/prefsquare),'g:',label='38F-24S', linewidth=2)
##pl.plot(frf_CB_xfem_2[0],10*log10(frf_CB_xfem_2[1]/prefsquare),'b-.',label='102F-32S', linewidth=2)
##pl.plot(frf_CB_xfem_1[0],10*log10(frf_CB_xfem_1[1]/prefsquare),'r--',label='250F-44S', linewidth=2)
##pl.axis([10.0, 200.0, 20, 90])
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Mean quadratic pressure (dB)')
##pl.grid('on')
##pl.legend(loc=4)
##
##pl.figure(2)
##pl.plot(array(lastfreqF)/200.0,globalerror,'rs--',label='Fluid basis', linewidth=1)
##pl.plot(array(lastfreqS)/200.0,globalerror,'bo--',label='Structure basis', linewidth=1)
##pl.figtext(0.225,0.8,'38F-24S')
##pl.figtext(0.45,0.29,'102F-32S')
##pl.figtext(0.76,0.15,'250F-44S')
##pl.axis([0.5, 2.5, 0.0, 1.3])
##pl.xlabel('Highest frequency / Highest observed frequency')
##pl.ylabel('Mean dB difference')
##pl.grid('on')
##pl.legend()
##
##
##fixedtime_ref=10*60-9*60-25  # 15:10:00 - 15:09:25
##totaltime_ref=fixedtime_ref+500.0*(16*3600.0+37*60.0+45.0-15*3600-10*60-0)/400.0 # 16:37:45 - 15:10:00 for 400 steps
##
###250F-44S
##fixedtime_xfem_1=172.86
##totaltime_xfem_1=fixedtime_xfem_1+500.0*(48*60+43-42*60-18)/400.0 # 15:48:43 - 15:42:18 for 400 steps
##
###102F-32S
##fixedtime_xfem_2=135.14
##totaltime_xfem_2=fixedtime_xfem_1+500.0*(36*60+12-31*60-34)/400.0 # 15:36:12 - 15:31:34 for 400 steps
##
###38F-24S
##fixedtime_xfem_3=119.5
##totaltime_xfem_3=fixedtime_xfem_3+500.0*(16*3600+0*60+54-15*3600-57*60-1)/400.0 # 16:00:54 - 15:57:01 for 400 steps
##
##time0=[fixedtime_ref,totaltime_ref]
##
##time1=[fixedtime_xfem_1,totaltime_xfem_1]
##time2=[fixedtime_xfem_2,totaltime_xfem_2]
##time3=[fixedtime_xfem_3,totaltime_xfem_3]
##
##nbsteps=[0,500]
##
##pl.figure(3)
##pl.plot(array(nbsteps),array(time1)/totaltime_ref,'r--',label='250F-44S', linewidth=3)
##pl.plot(array(nbsteps),array(time2)/totaltime_ref,'b-.',label='102F-32S', linewidth=3)
##pl.plot(array(nbsteps),array(time3)/totaltime_ref,'g:',label='38F-24S', linewidth=3)
##pl.plot(array(nbsteps),array(time0)/totaltime_ref,'k-',label='ref', linewidth=2)
##pl.axis([0, 500, 0.0, 0.1])
##pl.grid('on')
##pl.legend(loc=4)
##pl.xlabel('Nb. frequency steps')
##pl.ylabel('Cpu time / Reference cpu time')
##
###pl.figure(4)
###pl.plot(array([771+250+44+1,771+102+32+1,771+38+24+1]),array([totaltime_xfem_1,totaltime_xfem_2,totaltime_xfem_3])/totaltime_ref,'or-',label='250F-44S', linewidth=2)
###pl.axis([0, 1100, 0.0, 0.1])
##
##pl.show()
##
