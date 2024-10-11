import scipy
import pylab
import pickle

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_no_pre_stress_no_damping.pkl','rb')
frf_no_pre_stress_no_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_no_pre_stress_with_fractional_damping.pkl','rb')
frf_no_pre_stress_with_fractional_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_with_pre_stress1000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_1000_with_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_with_pre_stress2000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_2000_with_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_with_pre_stress3000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_3000_with_damping=pickle.load(f)
f.close()


f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_with_pre_stress4000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_4000_with_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_with_pre_stress5000_with_fractional_damping.pkl','rb')
frf_with_pre_stress_5000_with_damping=pickle.load(f)
f.close()


f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_rigid_faces_no_damping.pkl','rb')
frf_rigid_faces_no_pre_stress_no_damping=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_beam3_no_damping_10modes.pkl','rb')
frf_beam3_no_pre_stress_no_damping_10modes=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_beam3_no_damping_50modes.pkl','rb')
frf_beam3_no_pre_stress_no_damping_50modes=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_beam3_no_damping_100modes.pkl','rb')
frf_beam3_no_pre_stress_no_damping_100modes=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_beam3_no_damping_500modes.pkl','rb')
frf_beam3_no_pre_stress_no_damping_500modes=pickle.load(f)
f.close()

f=open('Resultats_pre-contrainte-maillage2/Results_assemblage_FRF_linear_beam3_no_damping_3000modes.pkl','rb')
frf_beam3_no_pre_stress_no_damping_3000modes=pickle.load(f)
f.close()


RefDisp=scipy.sqrt(2*(1e-4)**2)

pylab.figure(1)

pylab.plot(frf_no_pre_stress_no_damping[0],20*scipy.log10(frf_no_pre_stress_no_damping[1]/RefDisp),'k--',label='sans amortissement, sans precontrainte', linewidth=2)
pylab.plot(frf_no_pre_stress_with_fractional_damping[0],20*scipy.log10(frf_no_pre_stress_with_fractional_damping[1]/RefDisp),'k-',label='avec amortissement, sans precontrainte', linewidth=2)
#pylab.plot(frf_no_pre_stress_with_fractional_damping[0],scipy.log10(frf_no_pre_stress_with_fractional_damping[1]/RefDisp),'k-',label='with fractional derivative damping, no pre-stress', linewidth=2)
#pylab.plot(frf_with_pre_stress_1000_no_damping[0],scipy.log10(frf_with_pre_stress_1000_no_damping[1]),'g-',label='no damping, with pre-stress 1000N', linewidth=2)
#pylab.plot(frf_with_pre_stress_0_with_damping[0],20*scipy.log10(frf_with_pre_stress_0_with_damping[1]/RefDisp),'k-',label='with damping, with pre-stress 0 N', linewidth=2)
pylab.plot(frf_with_pre_stress_1000_with_damping[0],20*scipy.log10(frf_with_pre_stress_1000_with_damping[1]/RefDisp),'y-',label='avec amortissement, precontrainte de 1000 N', linewidth=2)
pylab.plot(frf_with_pre_stress_2000_with_damping[0],20*scipy.log10(frf_with_pre_stress_2000_with_damping[1]/RefDisp),'g-',label='avec amortissement, precontrainte de 2000 N', linewidth=2)
pylab.plot(frf_with_pre_stress_3000_with_damping[0],20*scipy.log10(frf_with_pre_stress_3000_with_damping[1]/RefDisp),'m-',label='avec amortissement, precontrainte de 3000 N', linewidth=2)
pylab.plot(frf_with_pre_stress_4000_with_damping[0],20*scipy.log10(frf_with_pre_stress_4000_with_damping[1]/RefDisp),'b-',label='avec amortissement, precontrainte de 4000 N', linewidth=2)
pylab.plot(frf_with_pre_stress_5000_with_damping[0],20*scipy.log10(frf_with_pre_stress_5000_with_damping[1]/RefDisp),'r-',label='avec amortissement, precontrainte de 5000 N', linewidth=2)
#pylab.plot(frf_with_pre_stress_6000_with_damping[0],scipy.log10(frf_with_pre_stress_6000_with_damping[1]/RefDisp),'r-',label='with damping, with pre-stress 6000 N', linewidth=2)
#pylab.plot(frf_with_pre_stress_8000_with_damping[0],scipy.log10(frf_with_pre_stress_8000_with_damping[1]/RefDisp),'m-',label='with damping, with pre-stress 8000N', linewidth=2)
#pylab.plot(frf_rigid_faces_no_pre_stress_no_damping[0],scipy.log10(frf_rigid_faces_no_pre_stress_no_damping[1]),'b-',label='rigid faces, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_no_pre_stress_no_damping_no_mass[0],scipy.log10(frf_no_pre_stress_no_damping_no_mass[1]),'ob-',label='no damping, no pre-stress, no link mass', linewidth=2)
#pylab.plot(frf_no_pre_stress_no_damping_rigid_faces[0],scipy.log10(frf_no_pre_stress_no_damping_rigid_faces[1]),'or-',label='no damping, no pre-stress, rigid faces', linewidth=2)
#pylab.plot(frf_no_pre_stress_with_beta_M_damping[0],scipy.log10(frf_no_pre_stress_with_beta_M_damping[1]),'or-',label='with beta-M damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam1_no_pre_stress_no_damping[0],scipy.log10(frf_beam1_no_pre_stress_no_damping[1]),'oy-',label='BEAM-1, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam2_no_pre_stress_no_damping[0],scipy.log10(frf_beam2_no_pre_stress_no_damping[1]),'og-',label='BEAM-2, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam3_no_pre_stress_no_damping[0],scipy.log10(frf_beam3_no_pre_stress_no_damping[1]),'r-',label='BEAM-CB-3, no damping, no pre-stress', linewidth=2)
#pylab.plot(frf_beam3_no_pre_stress_with_damping[0],scipy.log10(frf_beam3_no_pre_stress_with_damping[1]),'y-',label='BEAM-CB-3, with damping, no pre-stress', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequence (Hz)')
pylab.ylabel('Deplacement point P (dB)')
pylab.grid('on')
pylab.legend(loc='upper right')

pylab.figure(2)
pylab.plot(frf_no_pre_stress_no_damping[0],20*scipy.log10(frf_no_pre_stress_no_damping[1]/RefDisp),'k--',label='sans amortissement, sans precontrainte', linewidth=2)
pylab.plot(frf_no_pre_stress_with_fractional_damping[0],20*scipy.log10(frf_no_pre_stress_with_fractional_damping[1]/RefDisp),'k-',label='with fractional derivative damping, no pre-stress', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequence (Hz)')
pylab.ylabel('Deplacement point P (dB)')
pylab.grid('on')
pylab.legend(loc='upper right')

pylab.figure(3)
pylab.plot(frf_no_pre_stress_no_damping[0],20*scipy.log10(frf_no_pre_stress_no_damping[1]/RefDisp),'k-',label='sans amortissement, calcul complet', linewidth=2)
pylab.plot(frf_rigid_faces_no_pre_stress_no_damping[0],20*scipy.log10(frf_rigid_faces_no_pre_stress_no_damping[1]/RefDisp),'b-',label='sans amortissement, faces rigides', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequence (Hz)')
pylab.ylabel('Deplacement point P (dB)')
pylab.grid('on')
pylab.legend(loc='upper right')

pylab.figure(4)
pylab.plot(frf_rigid_faces_no_pre_stress_no_damping[0],20*scipy.log10(frf_rigid_faces_no_pre_stress_no_damping[1]/RefDisp),'k-',label='sans amortissement, faces rigides', linewidth=2)
pylab.plot(frf_beam3_no_pre_stress_no_damping_10modes[0],20*scipy.log10(frf_beam3_no_pre_stress_no_damping_10modes[1]/RefDisp),'b--',label='sans amortissement, 4 liaisons CB, 10 modes (1472 Hz)', linewidth=2)
pylab.plot(frf_beam3_no_pre_stress_no_damping_50modes[0],20*scipy.log10(frf_beam3_no_pre_stress_no_damping_50modes[1]/RefDisp),'g--',label='sans amortissement, 4 liaisons CB, 50 modes (2572 Hz)', linewidth=2)
pylab.plot(frf_beam3_no_pre_stress_no_damping_100modes[0],20*scipy.log10(frf_beam3_no_pre_stress_no_damping_100modes[1]/RefDisp),'m--',label='sans amortissement, 4 liaisons CB, 100 modes (3355 Hz)', linewidth=2)
pylab.plot(frf_beam3_no_pre_stress_no_damping_500modes[0],20*scipy.log10(frf_beam3_no_pre_stress_no_damping_500modes[1]/RefDisp),'y--',label='sans amortissement, 4 liaisons CB, 500 modes (6409 Hz)', linewidth=2)
#pylab.plot(frf_beam3_no_pre_stress_no_damping_3000modes[0],20*scipy.log10(frf_beam3_no_pre_stress_no_damping_3000modes[1]/RefDisp),'r-',label='sans amortissement, 4 liaisons CB, 3000 modes (15825 Hz)', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequence (Hz)')
pylab.ylabel('Deplacement point P (dB)')
pylab.grid('on')
pylab.legend(loc='upper right')


# pour le poster / septembre 2015
pylab.figure(5)
pylab.plot(frf_no_pre_stress_no_damping[0],20*scipy.log10(frf_no_pre_stress_no_damping[1]/RefDisp),'k-',label='liaison non amortie, non pre-contrainte', linewidth=3)
pylab.plot(frf_no_pre_stress_with_fractional_damping[0],20*scipy.log10(frf_no_pre_stress_with_fractional_damping[1]/RefDisp),'r-',label='liaison amortie, non pre-contrainte', linewidth=3)
pylab.plot(frf_with_pre_stress_3000_with_damping[0],20*scipy.log10(frf_with_pre_stress_4000_with_damping[1]/RefDisp),'b-',label='liaison amortie, avec pre-contrainte de 4000 N', linewidth=3)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequence (Hz)')
pylab.ylabel('Deplacement point P (dB)')
pylab.grid('on')
pylab.legend(loc='upper right')

pylab.show()
