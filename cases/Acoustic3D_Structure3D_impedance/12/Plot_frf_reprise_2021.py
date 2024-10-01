from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle


# xfem + POROUS + flexible structure / no reduction
f=open('results/cavity8_with_porous_air_flexible_structure_ang30_pour_article_results.frf','rb')
frf_xfem=pickle.load(f)
f.close()

# xfem + POROUS + flexible structure /  reduction air + structure
#f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang30_results.frf','rb')
#frf_xfem_reduc_air_struc=pickle.load(f)
#f.close()

# xfem + POROUS + flexible structure /  no reduction air + new porous
#f=open('results/cavity10_with_new_porous_air_flexible_structure_ang30_results.frf','rb')
#frf_xfem_new_porous=pickle.load(f)
#f.close()

# xfem + impedance paroi + flexible structure /  no reduction 
f=open('results/cavity12_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_imp_no_reduc=pickle.load(f)
f.close()

# xfem + impedance paroi + flexible structure /  no reduction 
f=open('results/cavity12_with_impedance_air_flexible_structure_CB_reduction_results.frf','rb')
frf_xfem_imp_CB_reduc=pickle.load(f)
f.close()

# xfem + NO impedance paroi + flexible structure /  no reduction 
f=open('results/cavity12_no_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_no_imp_no_reduc=pickle.load(f)
f.close()


prefsquare=20e-6*20e-6
print('moyenne sur la bande de frequence -   no  impedance [dB]',sum(10*log10(frf_xfem_no_imp_no_reduc[1]/prefsquare))/len(frf_xfem_no_imp_no_reduc[1]))
print('moyenne sur la bande de frequence - with impedance [dB]',sum(10*log10(frf_xfem_imp_no_reduc[1]/prefsquare))/len(frf_xfem_imp_no_reduc[1]))

pl.figure(1)
#pl.plot(frf_xfem[0],10*log10(frf_xfem[1]/prefsquare),'b-',label='xfem + flex. struc. + POROUS Biot-Allard/ 30 deg/ no reduction', linewidth=2)
#pl.plot(frf_xfem_reduc_air_struc[0],10*log10(frf_xfem_reduc_air_struc[1]/prefsquare),'ob--',label='xfem + flex. struc. + POROUS / 30 deg/ reduction', linewidth=2)
pl.plot(frf_xfem_imp_no_reduc[0],10*log10(frf_xfem_imp_no_reduc[1]/prefsquare),'g-',label='xfem + flex. struc. + imp. paroi / no reduc', linewidth=2)
pl.plot(frf_xfem_imp_CB_reduc[0],10*log10(frf_xfem_imp_CB_reduc[1]/prefsquare),'m-',label='xfem + flex. struc. + imp. paroi / CB reduc', linewidth=2)
#pl.plot(frf_xfem_no_imp_no_reduc[0],10*log10(frf_xfem_no_imp_no_reduc[1]/prefsquare),'r-',label='xfem + flex. struc. + NO imp. paroi / no reduc', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.show()



