import scipy
import pylab
import pickle

f=open('Results_truss_2_XE_1_NO_damping_static_no_damping.frf','rb')
frf_no_damping_XE_1=pickle.load(f)
f.close()

#f=open('Results_truss_damping.frf','rb')
#frf_damping=pickle.load(f)
#f.close()

f=open('Results_truss_2_XE_optim_NO_damping_static_no_damping.frf','rb')
frf_no_damping_XE_static_optim=pickle.load(f)
f.close()

#f=open('Results_truss_static_optim_damping.frf','rb')
#frf_damping_static_optim=pickle.load(f)
#f.close()

pylab.figure(1)
pylab.plot(frf_no_damping_XE_1[0],scipy.log10(frf_no_damping_XE_1[1]),'k-',label='no damping, XE=1', linewidth=2)
#pylab.plot(frf_damping[0],scipy.log10(frf_damping[1]),'b-',label='damping', linewidth=2)
pylab.plot(frf_no_damping_XE_static_optim[0],scipy.log10(frf_no_damping_XE_static_optim[1]),'g-',label='no damping, XE static optim', linewidth=2)
#pylab.plot(frf_damping_static_optim[0],scipy.log10(frf_damping_static_optim[1]),'r-',label='damping, static optim', linewidth=2)
#pl.axis([10.0, 200.0, 20, 90])
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Extremity displacement')
pylab.grid('on')
pylab.legend(loc=4)
pylab.show()
