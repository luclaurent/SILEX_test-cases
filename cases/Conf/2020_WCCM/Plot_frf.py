from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

f=open('results/cavity_acou3D_struc_3D_v3_1.00E+00_1.00E+00_5.00E-01_8.00E-01_results.frf','rb')
frf_no_reduc=pickle.load(f)
f.close()

f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient_1.00E+00_1.00E+00_NbModesFluid_30_5.00E-01_NbModesFluid_30_8.00E-01_NbModesFluid_30_results.frf','rb')
frf_reduc_gradient1=pickle.load(f)
f.close()

f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient_1.00E+00_1.00E+00_NbModesFluid_150_5.00E-01_NbModesFluid_150_8.00E-01_NbModesFluid_150_results.frf','rb')
frf_reduc_gradient2=pickle.load(f)
f.close()

f=open('results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient_1.00E+00_1.00E+00_NbModesFluid_300_5.00E-01_NbModesFluid_300_8.00E-01_NbModesFluid_300_results.frf','rb')
frf_reduc_gradient3=pickle.load(f)
f.close()


prefsquare=20e-6*20e-6

pl.figure(1)
pl.plot(frf_no_reduc[0],10*log10(frf_no_reduc[1]/prefsquare),'k-',label='No reduction', linewidth=1)
pl.plot(frf_reduc_gradient1[0],10*log10(frf_reduc_gradient1[1]/prefsquare),'r-',label='CB reduction, 30 modes, f_max=128 Hz', linewidth=1)
pl.plot(frf_reduc_gradient2[0],10*log10(frf_reduc_gradient2[1]/prefsquare),'g-',label='CB reduction, 150 modes, f_max=267 Hz', linewidth=1)
pl.plot(frf_reduc_gradient3[0],10*log10(frf_reduc_gradient3[1]/prefsquare),'b-',label='CB reduction, 300 modes, f_max=387 Hz', linewidth=1)
pl.axis([1.0, 130.0, 70, 125])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
pl.plot(frf_no_reduc[0],frf_no_reduc[2],'k-',label='No reduction', linewidth=1)
pl.plot(frf_reduc_gradient1[0],frf_reduc_gradient1[2],'r-',label='CB reduction, 30 modes, f_max=128 Hz', linewidth=1)
pl.plot(frf_reduc_gradient2[0],frf_reduc_gradient2[2],'g-',label='CB reduction, 150 modes, f_max=267 Hz', linewidth=1)
pl.plot(frf_reduc_gradient3[0],frf_reduc_gradient3[2],'b-',label='CB reduction, 300 modes, f_max=387 Hz', linewidth=1)
pl.axis([1.0, 130.0, -800, 800])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume for X pos. (dB/m)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(3)
pl.plot(frf_no_reduc[0],frf_no_reduc[5],'k-',label='No reduction', linewidth=1)
pl.plot(frf_reduc_gradient1[0],frf_reduc_gradient1[5],'r-',label='CB reduction, 30 modes, f_max=128 Hz', linewidth=1)
pl.plot(frf_reduc_gradient2[0],frf_reduc_gradient2[5],'g-',label='CB reduction, 150 modes, f_max=267 Hz', linewidth=1)
pl.plot(frf_reduc_gradient3[0],frf_reduc_gradient3[5],'b-',label='CB reduction, 300 modes, f_max=387 Hz', linewidth=1)
pl.axis([1.0, 130.0, -800, 800])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume for R (dB/m)')
pl.grid('on')
pl.legend(loc=4)

Results=[[30,128.101,160.08985700000002,0.5372142852],
         [150,267.038,283.821646,0.9524],
         [300,386.695,438.306018,1.470825563]
         ]

# time for one time step with no reduction
TimeRef_OneFreq=3.1450526107382

NbFluidModes=[]
LastFreq=[]
CPUtime=[]
CPUtimePerStep=[]
for ii in Results:
    NbFluidModes.append(ii[0])
    LastFreq.append(ii[1])
    CPUtime.append(ii[2])
    CPUtimePerStep.append(ii[3])


pl.figure(4)
pl.plot(NbFluidModes,scipy.array(CPUtimePerStep)/TimeRef_OneFreq,'-o')
pl.axis([0, 300, 0, 1])
pl.grid()
pl.xlabel('Nb Fluid Modes')
pl.ylabel('Relative CPU time per time step')

pl.show()

##pl.figure(4)
##pl.plot(frf_no_reduc[0],frf_no_reduc[3],'k-',label='no reduction, gradient pos. y', linewidth=1)
##pl.plot(frf_reduc_gradient1[0],frf_reduc_gradient1[3],'r-',label='reduction, gradient pos. y, 30 modes', linewidth=1)
##pl.plot(frf_reduc_gradient2[0],frf_reduc_gradient2[3],'g-',label='reduction, gradient pos. y, 150 modes', linewidth=1)
##pl.plot(frf_reduc_gradient3[0],frf_reduc_gradient3[3],'b-',label='reduction, gradient pos. y, 300 modes', linewidth=1)
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Pressure gradient in control volume')
##pl.grid('on')
##pl.legend(loc=4)

##pl.figure(5)
##pl.plot(frf_no_reduc[0],frf_no_reduc[4],'k-',label='no reduction, gradient pos. z', linewidth=1)
##pl.plot(frf_reduc_gradient1[0],frf_reduc_gradient1[4],'r-',label='reduction, gradient pos. z, 30 modes', linewidth=1)
##pl.plot(frf_reduc_gradient2[0],frf_reduc_gradient2[4],'g-',label='reduction, gradient pos. z, 150 modes', linewidth=1)
##pl.plot(frf_reduc_gradient3[0],frf_reduc_gradient3[4],'b-',label='reduction, gradient pos. z, 300 modes', linewidth=1)
##pl.xlabel('Frequency (Hz)')
##pl.ylabel('Pressure gradient in control volume')
##pl.grid('on')
##pl.legend(loc=4)

