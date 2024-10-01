from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

import os

# xfem + POROUS + flexible structure / no reduction
#f=open('results/cavity8_with_porous_air_flexible_structure_ang30_pour_article_results.frf','rb')
#frf_xfem=pickle.load(f)
#f.close()

# xfem + POROUS + flexible structure /  reduction air + structure
#f=open('results/cavity8_with_porous_air_flexible_structure_CB_ang30_results.frf','rb')
#frf_xfem_reduc_air_struc=pickle.load(f)
#f.close()

# xfem + POROUS + flexible structure /  no reduction air + new porous
#f=open('results/cavity10_with_new_porous_air_flexible_structure_ang30_results.frf','rb')
#frf_xfem_new_porous=pickle.load(f)
#f.close()

# xfem + impedance paroi + flexible structure /  no reduction 
f=open(os.path.join(os.path.dirname(__file__),'results/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc=pickle.load(f)
f.close()
# f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_40/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
# frf_xfem_imp_no_reduc_mesh_40=pickle.load(f)
# f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_20_angle_0/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc_mesh_20_angle_0=pickle.load(f)
f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_20_angle_60/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc_mesh_20_angle_60=pickle.load(f)
f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_30_angle_60/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc_mesh_30_angle_60=pickle.load(f)
f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_40_angle_60/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc_mesh_40_angle_60=pickle.load(f)
f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_40_angle_0/cavity_with_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_imp_no_reduc_mesh_40_angle_0=pickle.load(f)
f.close()
f=open('results/Mesh_80_angle_60/cavity_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_imp_no_reduc_mesh_80_angle_60=pickle.load(f)
f.close()

### xfem + impedance paroi + flexible structure /  no reduction 
# f=open(os.path.join(os.path.dirname(__file__),'results/cavity_with_impedance_air_flexible_structure_CB_reduction_results.frf'),'rb')
# frf_xfem_imp_CB_reduc=pickle.load(f)
# f.close()

# xfem + NO impedance paroi + flexible structure /  no reduction 
f=open(os.path.join(os.path.dirname(__file__),'results/cavity_no_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_no_imp_no_reduc=pickle.load(f)
f.close()
# f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_40/cavity_no_impedance_air_flexible_structure_results.frf'),'rb')
# frf_xfem_no_imp_no_reduc_mesh_40=pickle.load(f)
# f.close()
#f=open('results/Mesh_20/cavity_no_impedance_air_flexible_structure_results.frf','rb')
#frf_xfem_no_imp_no_reduc_mesh_20=pickle.load(f)
#f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_20_angle_0/cavity_no_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_no_imp_no_reduc_mesh_20_angle_0=pickle.load(f)
f.close()
f=open(os.path.join(os.path.dirname(__file__),'results/Mesh_20_angle_60/cavity_no_impedance_air_flexible_structure_results.frf'),'rb')
frf_xfem_no_imp_no_reduc_mesh_20_angle_60=pickle.load(f)
f.close()

# xfem + NO impedance paroi + rigid structure / no reduction
##f=open('results/cavity13_no_imp_rigid_sphere_structure_no_reduc_results.frf','rb')
##frf_xfem_13=pickle.load(f)
##f.close()




f=open('results/Tet10_flexwood_foamWalid_angle_60_h_25/cavity_tet10_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet10_angle_60_h_25=pickle.load(f)
f.close()
f=open('results/Tet10_flexwood_foamWalid_angle_60_h_20/cavity_tet10_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet10_angle_60_h_20=pickle.load(f)
f.close()
f=open('results/Tet10_flexwood_foamWalid_angle_60_h_15/cavity_tet10_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet10_angle_60_h_15=pickle.load(f)
f.close()
f=open('results/Tet10_flexwood_foamWalid_angle_60_h_10/cavity_tet10_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet10_angle_60_h_10=pickle.load(f)
f.close()

f=open('results/cavity_with_impedance_air_rigid_structure_tet4_results.frf','rb')
frf_xfem_rigid_imp_tet4=pickle.load(f)
f.close()

f=open('results/cavity_tet10_with_impedance_air_rigid_structure_results.frf','rb')
frf_xfem_rigid_imp_tet10=pickle.load(f)
f.close()

f=open('results/Tet4_flexwood_foamWalid_angle_60_h_20/cavity_tet4_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet4_angle_60_h_20=pickle.load(f)
f.close()

f=open('results/Tet4_flexwood_foamWalid_angle_60_h_30/cavity_tet4_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet4_angle_60_h_30=pickle.load(f)
f.close()

f=open('results/cavity_tet10_with_impedance_air_flexible_structure_results.frf','rb')
frf_xfem_tet10=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6
n1=300
n2=399
print('moyenne sur la bande de frequence entre le pas ',n1,' et le pas ',n2,' -   no  impedance [dB]',sum(10*log10(frf_xfem_no_imp_no_reduc[1][n1:n2]/prefsquare))/len(frf_xfem_no_imp_no_reduc[1][n1:n2]))
print('moyenne sur la bande de frequence entre le pas ',n1,' et le pas ',n2,' - with impedance [dB]',sum(10*log10(frf_xfem_imp_no_reduc[1][n1:n2]/prefsquare))/len(frf_xfem_imp_no_reduc[1][n1:n2]))

pl.figure(1)
#pl.plot(frf_xfem_imp_no_reduc_mesh_40[0],10*log10(frf_xfem_imp_no_reduc_mesh_40[1]/prefsquare),'r-',label='xfem + flex. struc. + imp. paroi / no reduc / mesh 40', linewidth=2)
#pl.plot(frf_xfem_imp_no_reduc_mesh_20[0],10*log10(frf_xfem_imp_no_reduc_mesh_20[1]/prefsquare),'g-',label='xfem + flex. struc. + imp. paroi / no reduc / mesh 20', linewidth=2)
#pl.plot(frf_xfem_no_imp_no_reduc_mesh_40[0],10*log10(frf_xfem_no_imp_no_reduc_mesh_40[1]/prefsquare),'b-',label='xfem + flex. struc. + NO imp. paroi / no reduc / mesh 40', linewidth=2)
#pl.plot(frf_xfem_no_imp_no_reduc_fine_mesh[0],10*log10(frf_xfem_no_imp_no_reduc_fine_mesh[1]/prefsquare),'r-',label='xfem + flex. struc. + NO imp. paroi / no reduc / fine mesh', linewidth=2)
#pl.plot(frf_xfem_imp_no_reduc[0],10*log10(frf_xfem_imp_no_reduc[1]/prefsquare),'b-',label='xfem + flex. struc. + imp. paroi / no reduc', linewidth=2)
#pl.plot(frf_xfem_no_imp_no_reduc[0],10*log10(frf_xfem_no_imp_no_reduc[1]/prefsquare),'r--',label='xfem + flex. struc. + NO imp. paroi / no reduc', linewidth=2)
#pl.plot(frf_xfem[0],10*log10(frf_xfem[1]/prefsquare),'b-',label='xfem + flex. struc. + POROUS Biot-Allard/ 30 deg/ no reduction', linewidth=2)
#pl.plot(frf_xfem_reduc_air_struc[0],10*log10(frf_xfem_reduc_air_struc[1]/prefsquare),'ob--',label='xfem + flex. struc. + POROUS / 30 deg/ reduction', linewidth=2)
#pl.plot(frf_xfem_imp_CB_reduc[0],10*log10(frf_xfem_imp_CB_reduc[1]/prefsquare),'m-',label='xfem + flex. struc. + imp. paroi / CB reduc', linewidth=2)
#pl.plot(frf_xfem_13[0],10*log10(frf_xfem_13[1]/prefsquare),'y-',label='xfem + rigid struc. + NO imp. / no reduc', linewidth=2)
#pl.axis([1.0, 120.0, 70, 105])

#pl.plot(frf_xfem_imp_no_reduc_mesh_20_angle_0[0],10*log10(frf_xfem_imp_no_reduc_mesh_20_angle_0[1]/prefsquare),'b-',label='xfem + flex. struc. + imp. paroi / no reduc / angle 0', linewidth=1)
#pl.plot(frf_xfem_no_imp_no_reduc_mesh_20_angle_0[0],10*log10(frf_xfem_no_imp_no_reduc_mesh_20_angle_0[1]/prefsquare),'b--',label='xfem + flex. struc. + NO imp. paroi / no reduc / angle 0', linewidth=1)
#pl.plot(frf_xfem_imp_no_reduc_mesh_40_angle_0[0],10*log10(frf_xfem_imp_no_reduc_mesh_40_angle_0[1]/prefsquare),'y-',label='xfem + flex. struc. + imp. paroi / no reduc / angle 0 / mesh 40', linewidth=1)

#pl.plot(frf_xfem_no_imp_no_reduc_mesh_20_angle_60[0],10*log10(frf_xfem_no_imp_no_reduc_mesh_20_angle_60[1]/prefsquare),'r--',label='xfem + flex. struc. + NO imp. paroi / no reduc/ angle 60', linewidth=1)

#pl.plot(frf_xfem_imp_no_reduc_mesh_20_angle_60[0],10*log10(frf_xfem_imp_no_reduc_mesh_20_angle_60[1]/prefsquare),'r-',label='xfem + flex. struc. + imp. paroi / no reduc/ angle 60 / mesh 20', linewidth=1)
#pl.plot(frf_xfem_imp_no_reduc_mesh_30_angle_60[0],10*log10(frf_xfem_imp_no_reduc_mesh_30_angle_60[1]/prefsquare),'m-',label='xfem + flex. struc. + imp. paroi / no reduc/ angle 60 / mesh 30', linewidth=1)
#pl.plot(frf_xfem_imp_no_reduc_mesh_40_angle_60[0],10*log10(frf_xfem_imp_no_reduc_mesh_40_angle_60[1]/prefsquare),'g-',label='xfem + flex. struc. + imp. paroi / no reduc/ angle 60 / mesh 40', linewidth=1)

#pl.plot(frf_xfem_imp_no_reduc_mesh_80_angle_60[0],10*log10(frf_xfem_imp_no_reduc_mesh_80_angle_60[1]/prefsquare),'k-',label='xfem + flex. struc. + imp. paroi / no reduc/ angle 60 / mesh 80', linewidth=1)

pl.plot(frf_xfem_tet10[0],10*log10(frf_xfem_tet10[1]/prefsquare),'k-',label='xfem + flex. struc. + imp. paroi / TET10', linewidth=1)

#pl.plot(frf_xfem_rigid_imp_tet4[0],10*log10(frf_xfem_rigid_imp_tet4[1]/prefsquare),'g-',label='xfem + rigid. struc. + imp. paroi / TET4', linewidth=1)
#pl.plot(frf_xfem_rigid_imp_tet10[0],10*log10(frf_xfem_rigid_imp_tet10[1]/prefsquare),'b-',label='xfem + rigid. struc. + imp. paroi / TET10', linewidth=1)
pl.plot(frf_xfem_tet10_angle_60_h_25[0],10*log10(frf_xfem_tet10_angle_60_h_25[1]/prefsquare),'r-',label='xfem + flex. wood struc. + imp. paroi Walid / TET10 / h 25 / 75496 nodes', linewidth=1)
pl.plot(frf_xfem_tet10_angle_60_h_20[0],10*log10(frf_xfem_tet10_angle_60_h_20[1]/prefsquare),'b-',label='xfem + flex. wood struc. + imp. paroi Walid / TET10 / h 20 / 41602 nodes', linewidth=1)
pl.plot(frf_xfem_tet10_angle_60_h_15[0],10*log10(frf_xfem_tet10_angle_60_h_15[1]/prefsquare),'m-',label='xfem + flex. wood struc. + imp. paroi Walid / TET10 / h 15 / 19983 nodes', linewidth=1)
#pl.plot(frf_xfem_tet10_angle_60_h_10[0],10*log10(frf_xfem_tet10_angle_60_h_10[1]/prefsquare),'g-',label='xfem + flex. wood struc. + imp. paroi Walid / TET10 / h 10 / 7881 nodes', linewidth=1)
#pl.plot(frf_xfem_tet4_angle_60_h_20[0],10*log10(frf_xfem_tet4_angle_60_h_20[1]/prefsquare),'y-',label='xfem + flex. wood struc. + imp. paroi Walid / TET4 / h 20 /  5657 nodes', linewidth=1)
#pl.plot(frf_xfem_tet4_angle_60_h_30[0],10*log10(frf_xfem_tet4_angle_60_h_30[1]/prefsquare),'c-',label='xfem + flex. wood struc. + imp. paroi Walid / TET4 / h 30 /  16689 nodes', linewidth=1)

pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

##pl.figure(2)
###pl.plot(frf_xfem[0],10*log10(frf_xfem[1]/prefsquare),'b-',label='xfem + flex. struc. + POROUS Biot-Allard/ 30 deg/ no reduction', linewidth=2)
###pl.plot(frf_xfem_reduc_air_struc[0],10*log10(frf_xfem_reduc_air_struc[1]/prefsquare),'ob--',label='xfem + flex. struc. + POROUS / 30 deg/ reduction', linewidth=2)
###pl.plot(frf_xfem_imp_no_reduc[0],10*log10(frf_xfem_imp_no_reduc[1]/prefsquare),'g-',label='xfem + flex. struc. + imp. paroi / no reduc', linewidth=2)
###pl.plot(frf_xfem_no_imp_no_reduc[0],10*log10(frf_xfem_no_imp_no_reduc[1]/prefsquare),'r-',label='xfem + flex. struc. + NO imp. paroi / no reduc', linewidth=2)
###pl.axis([1.0, 120.0, 70, 105])
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Mean quadratic pressure (dB)')
##pl.grid('on')
##pl.legend(loc=4)

pl.show()


##pl.figure(2)
###pl.plot(frf_xfem[0],10*log10(frf_xfem[1]/prefsquare),'b-',label='xfem + flex. struc. + POROUS Biot-Allard/ 30 deg/ no reduction', linewidth=2)
###pl.plot(frf_xfem_reduc_air_struc[0],10*log10(frf_xfem_reduc_air_struc[1]/prefsquare),'ob--',label='xfem + flex. struc. + POROUS / 30 deg/ reduction', linewidth=2)
###pl.plot(frf_xfem_imp_no_reduc[0],10*log10(frf_xfem_imp_no_reduc[1]/prefsquare),'g-',label='xfem + flex. struc. + imp. paroi / no reduc', linewidth=2)
###pl.plot(frf_xfem_imp_CB_reduc[0],10*log10(frf_xfem_imp_CB_reduc[1]/prefsquare),'m-',label='xfem + flex. struc. + imp. paroi / CB reduc', linewidth=2)
##pl.plot(frf_xfem_no_imp_no_reduc[0],10*log10(frf_xfem_imp_CB_reduc[1]/prefsquare)-10*log10(frf_xfem_imp_no_reduc[1]/prefsquare),'r-',label='xfem + flex. struc. + NO imp. paroi / no reduc', linewidth=2)
###pl.axis([1.0, 120.0, 70, 105])
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Mean quadratic pressure (dB)')
##pl.grid('on')
##pl.legend(loc=4)
##
##pl.show()



