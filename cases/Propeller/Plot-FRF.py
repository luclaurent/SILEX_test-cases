import scipy
import pylab
import pickle

f=open('Results_propeller_tet10_no_damping.frf','rb')
frf_tet10_no_damping=pickle.load(f)
f.close()

f=open('Results_propeller_tet4_no_damping.frf','rb')
frf_tet4_no_damping=pickle.load(f)
f.close()

pylab.figure(1)
pylab.plot(frf_tet10_no_damping[0],scipy.log10(frf_tet10_no_damping[1]),'k-',label='no damping, tet10', linewidth=2)
pylab.plot(frf_tet4_no_damping[0],scipy.log10(frf_tet4_no_damping[1]),'b-',label='no damping, tet4', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Extremity displacement')
pylab.grid('on')
pylab.legend(loc=4)
pylab.show()
