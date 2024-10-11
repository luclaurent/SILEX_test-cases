import scipy
import pylab
import pickle

f=open('Results_assemblage_FRF_linear_no_pre_stress_no_damping.pkl','rb')
frf_no_pre_stress_no_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_no_pre_stress_no_damping_no_mass.pkl','rb')
frf_no_pre_stress_no_damping_no_mass=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_no_pre_stress_no_damping_rigid_faces.pkl','rb')
frf_no_pre_stress_no_damping_rigid_faces=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_no_pre_stress_with_fractional_damping.pkl','rb')
frf_no_pre_stress_with_fractional_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_no_pre_stress_with_beta_M_damping.pkl','rb')
frf_no_pre_stress_with_beta_M_damping=pickle.load(f)
f.close()

##f=open('Results_assemblage_FRF_linear_with_pre_stress1000_no_damping.pkl','rb')
##frf_with_pre_stress_1000_no_damping=pickle.load(f)
##f.close()

f=open('Results_assemblage_FRF_linear_with_pre_stress1000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_1000_with_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_with_pre_stress2000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_2000_with_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_with_pre_stress4000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_4000_with_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_with_pre_stress6000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_6000_with_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_with_pre_stress8000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_8000_with_damping=pickle.load(f)
f.close()


f=open('Results_assemblage_FRF_linear_beam1_no_damping.pkl','rb')
frf_beam1_no_pre_stress_no_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_beam2_no_damping.pkl','rb')
frf_beam2_no_pre_stress_no_damping=pickle.load(f)
f.close()


f=open('Results_assemblage_FRF_linear_beam3_no_damping.pkl','rb')
frf_beam3_no_pre_stress_no_damping=pickle.load(f)
f.close()

f=open('Results_assemblage_FRF_linear_beam3_with_damping.pkl','rb')
frf_beam3_no_pre_stress_with_damping=pickle.load(f)
f.close()


f=open('Results_assemblage_FRF_linear_rigid_faces_no_damping.pkl','rb')
frf_rigid_faces_no_pre_stress_no_damping=pickle.load(f)
f.close()


RefDisp=scipy.sqrt(2e-4**2)

pylab.figure(1)
#pylab.plot(frf_no_pre_stress_no_damping[0],scipy.log10(frf_no_pre_stress_no_damping[1]),'k-',label='no damping, no pre-stress', linewidth=2)
pylab.plot(frf_no_pre_stress_with_fractional_damping[0],scipy.log10(frf_no_pre_stress_with_fractional_damping[1]/RefDisp),'k-',label='with fractional derivative damping, no pre-stress', linewidth=2)
#pylab.plot(frf_with_pre_stress_1000_no_damping[0],scipy.log10(frf_with_pre_stress_1000_no_damping[1]),'g-',label='no damping, with pre-stress 1000N', linewidth=2)
pylab.plot(frf_with_pre_stress_1000_with_damping[0],scipy.log10(frf_with_pre_stress_1000_with_damping[1]/RefDisp),'y-',label='with damping, with pre-stress 1000N', linewidth=2)
pylab.plot(frf_with_pre_stress_2000_with_damping[0],scipy.log10(frf_with_pre_stress_2000_with_damping[1]/RefDisp),'g-',label='with damping, with pre-stress 2000N', linewidth=2)
pylab.plot(frf_with_pre_stress_4000_with_damping[0],scipy.log10(frf_with_pre_stress_4000_with_damping[1]/RefDisp),'b-',label='with damping, with pre-stress 4000N', linewidth=2)
pylab.plot(frf_with_pre_stress_6000_with_damping[0],scipy.log10(frf_with_pre_stress_6000_with_damping[1]/RefDisp),'r-',label='with damping, with pre-stress 6000N', linewidth=2)
pylab.plot(frf_with_pre_stress_8000_with_damping[0],scipy.log10(frf_with_pre_stress_8000_with_damping[1]/RefDisp),'m-',label='with damping, with pre-stress 8000N', linewidth=2)
#pylab.plot(frf_rigid_faces_no_pre_stress_no_damping[0],scipy.log10(frf_rigid_faces_no_pre_stress_no_damping[1]),'b-',label='rigid faces, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_no_pre_stress_no_damping_no_mass[0],scipy.log10(frf_no_pre_stress_no_damping_no_mass[1]),'ob-',label='no damping, no pre-stress, no link mass', linewidth=2)
#pylab.plot(frf_no_pre_stress_no_damping_rigid_faces[0],scipy.log10(frf_no_pre_stress_no_damping_rigid_faces[1]),'or-',label='no damping, no pre-stress, rigid faces', linewidth=2)
#pylab.plot(frf_no_pre_stress_with_beta_M_damping[0],scipy.log10(frf_no_pre_stress_with_beta_M_damping[1]),'or-',label='with beta-M damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam1_no_pre_stress_no_damping[0],scipy.log10(frf_beam1_no_pre_stress_no_damping[1]),'oy-',label='BEAM-1, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam2_no_pre_stress_no_damping[0],scipy.log10(frf_beam2_no_pre_stress_no_damping[1]),'og-',label='BEAM-2, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam3_no_pre_stress_no_damping[0],scipy.log10(frf_beam3_no_pre_stress_no_damping[1]),'r-',label='BEAM-CB-3, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam3_no_pre_stress_with_damping[0],scipy.log10(frf_beam3_no_pre_stress_with_damping[1]),'y-',label='BEAM-CB-3, with damping, no pre-stress', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Extremity displacement')
pylab.grid('on')
pylab.legend(loc='upper right')
pylab.show()
