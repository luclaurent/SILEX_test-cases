from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

# "fine":
# // size of elements
# h = lx1/50; air
# h2 = ee/2.5;
# h = lx1/70; structure
# nb_mode_F = 200
# nb_mode_S = 100
# freq_ini     = 10.0
# freq_end     = 140.0
# nb_freq_step_per_proc=80 / 20 proc.
# Last eigenfrequency in fluid basis:  268.614567554
# Last eigenfrequency in structure basis:  310.121842914
# nb of total dofs:  6760 / CB+XFEM
# Number of nodes: 38496
# Number of elements in air: 198102
# Number of elements in porous: 6416
# Number of nodes in air: 35677
# Number of nodes in porous: 1913
# Number of nodes at interface: 702
# nnodes for structure= 3294
# nelem for structure= 6373

##CB+XFEM:
##time at the beginning of the computation: Mon Apr 24 18:52:05 2017
##time for computing structure: 0.33467
##time for computing the structure modal basis: 13.516658
##time to compute level set: 1.0327529999999996
##time to compute tangent level set: 0.15114099999999908
##time to find surface enriched elements: 6.504829000000001
##time to find edge enriched elements: 0.8139910000000015
##time to compute coupling matrices: 4.290914999999998
##time to compute Heaviside enrichment: 0.40980300000000014
##time to compute edge enrichment: 2.410295999999999
##LAST fluid eigen frequencies :  268.614567554
##time for computing the fluid modes: 113.59602900000002
##Compute PSI_IA
##LU decomposition has been made
##End of PSI_IA computing, send to other proc.
##time to compute PSI_IA: 29.478661999999986
##Compute PSI_IB
##End of PSI_IB computing, send to other proc.
##time to compute PSI_FB: 11.035554000000019
##time to compute Khat_hat: 0.5784010000000137
##time to compute Mstar: 6.7018170000000055
##time to compute Mhat: 4.536851000000013
##Proc.  1  / time at the beginning of the FRF: Mon Apr 24 18:55:59 2017
##Proc.  1  / time at the end of the FRF: Mon Apr 24 19:00:37 2017

##NO reduction / NO porous
##time at the beginning of the computation: Mon Apr 24 19:00:42 2017
##Proc.  1  / time at the beginning of the FRF: Mon Apr 24 19:01:02 2017
##Proc.  1  / time at the end of the FRF: Mon Apr 24 19:09:06 2017

##NO REDUCTION / With POROUS
##time at the beginning of the computation: Mon Apr 24 18:41:44 2017
##Proc.  1  / time at the beginning of the FRF: Mon Apr 24 18:42:04 2017
##Proc.  1  / time at the end of the FRF: Mon Apr 24 18:51:56 2017
##nb of total dofs:  60742


# xfem + NO POROUS + flexible structure / no reduction
f=open('results/cavity8_no_porous_air_flexible_structure_ang30_results.frf','rb')
frf_xfem_no_porous_30=pickle.load(f)
f.close()
f=open('results/cavity8_no_porous_air_flexible_structure_ang120_results.frf','rb')
frf_xfem_no_porous_120=pickle.load(f)
f.close()
f=open('results/cavity8_no_porous_air_flexible_structure_ang30_fine_results.frf','rb')
frf_xfem_no_porous_30_fine=pickle.load(f)
f.close()
f=open('results/cavity8_no_porous_air_flexible_structure_ang120_fine_results.frf','rb')
frf_xfem_no_porous_120_fine=pickle.load(f)
f.close()

# xfem + porous + flexible structure / no reduction
f=open('results/cavity8_with_porous_air_flexible_structure_ang120_results.frf','rb')
frf_xfem_120=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_ang30_results.frf','rb')
frf_xfem_30=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_ang120_fine_results.frf','rb')
frf_xfem_120_fine=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_ang30_fine_results.frf','rb')
frf_xfem_30_fine=pickle.load(f)
f.close()

##f=open('../classic/results/cavity7_with_porous_results.frf','rb')
##frf_classic=pickle.load(f)
##f.close()

f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang120_results.frf','rb')
frf_xfem_CB_120=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang30_results.frf','rb')
frf_xfem_CB_30=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang120_fine_results.frf','rb')
frf_xfem_CB_120_fine=pickle.load(f)
f.close()
f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang30_fine_results.frf','rb')
frf_xfem_CB_30_fine=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_xfem_no_porous_30_fine[0],10*log10(frf_xfem_no_porous_30_fine[1]/prefsquare),'r-',label='FINE / xfem + flex. struc. NO POROUS / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_30_fine[0],10*log10(frf_xfem_30_fine[1]/prefsquare),'b-',label='FINE / xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_no_porous_30[0],10*log10(frf_xfem_no_porous_30[1]/prefsquare),'k-',label='xfem + flex. struc. NO POROUS / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_30[0],10*log10(frf_xfem_30[1]/prefsquare),'g-',label='xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
pl.plot(frf_xfem_no_porous_120_fine[0],10*log10(frf_xfem_no_porous_120_fine[1]/prefsquare),'r-',label='FINE / xfem + flex. struc. NO POROUS / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_120_fine[0],10*log10(frf_xfem_120_fine[1]/prefsquare),'b-',label='FINE / xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_no_porous_120[0],10*log10(frf_xfem_no_porous_120[1]/prefsquare),'k-',label='xfem + flex. struc. NO POROUS / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_120[0],10*log10(frf_xfem_120[1]/prefsquare),'g-',label='xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(3)
pl.plot(frf_xfem_120_fine[0],10*log10(frf_xfem_120_fine[1]/prefsquare),'r-',label='FINE / xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_30_fine[0],10*log10(frf_xfem_30_fine[1]/prefsquare),'b-',label='FINE / xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_120[0],10*log10(frf_xfem_120[1]/prefsquare),'k-',label='xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_30[0],10*log10(frf_xfem_30[1]/prefsquare),'g-',label='xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(4)
pl.plot(frf_xfem_30_fine[0],10*log10(frf_xfem_30_fine[1]/prefsquare),'r-',label='FINE / xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_CB_30_fine[0],10*log10(frf_xfem_CB_30_fine[1]/prefsquare),'b-',label='FINE / xfem + flex. struc. / 30 deg / CB reduction', linewidth=2)
pl.plot(frf_xfem_30[0],10*log10(frf_xfem_30[1]/prefsquare),'k-',label='xfem + flex. struc. / 30 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_CB_30[0],10*log10(frf_xfem_CB_30[1]/prefsquare),'g-',label='xfem + flex. struc. / 30 deg / CB reduction', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(5)
pl.plot(frf_xfem_120_fine[0],10*log10(frf_xfem_120_fine[1]/prefsquare),'r-',label='FINE / xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_CB_120_fine[0],10*log10(frf_xfem_CB_120_fine[1]/prefsquare),'b-',label='FINE / xfem + flex. struc. / 120 deg / CB reduction', linewidth=2)
pl.plot(frf_xfem_120[0],10*log10(frf_xfem_120[1]/prefsquare),'k-',label='xfem + flex. struc. / 120 deg/ no reduction', linewidth=2)
pl.plot(frf_xfem_CB_120[0],10*log10(frf_xfem_CB_120[1]/prefsquare),'g-',label='xfem + flex. struc. / 120 deg / CB reduction', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(6)
pl.plot(frf_xfem_30_fine[0],10*log10(frf_xfem_30_fine[1]/prefsquare),'g-', linewidth=2)
#pl.plot(frf_xfem_CB_30_fine[0],10*log10(frf_xfem_CB_30_fine[1]/prefsquare),'g-', linewidth=2)
pl.axis([20.0, 140.0, 70, 130])
pl.xlabel('Frequence (Hz)')
pl.ylabel('Pression quadratique moyenne (dB)')
pl.grid('on')
pl.legend(loc=4)


pl.show()



