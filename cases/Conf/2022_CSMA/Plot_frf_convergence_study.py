from numpy import *
import string
import time
import scipy

import pylab as pl
import pickle

f=open('results/cavity_acou3D_struc_3D_v3_1.00E+00_1.00E+00_5.00E-01_8.00E-01_results.frf','rb')
frf_no_reduc=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

Results=[[10,74.229,144.03,0.4833],
         [20,103.199,154.67612499,0.51904739],
         [30,128.101,160.08985700000002,0.5372142852],
         [40,147.007,178.416178,0.598712],
         [50,160.30,185.939574,0.62395830],
         [60,173.222,196.314997,0.6587751577],
         [70,185.086,200.66522,0.6733732],
         [80,198.387,209.382375,0.702625],
         [90,211.838,217.814705,0.7309218288],
         [100,222.6170,233.786786,0.78451941],
         [110,234.292,243.58025,0.81738338],
         [120,242.798,257.784,0.8650],
         [130,251.3540,267.251,0.8968],
         [140,259.4565,276.147115,0.9266],
         [150,267.038,283.821646,0.9524],
         [200,312.398,336.869284,1.130433],
         [250,352.592,387.319332,1.299729],
         [300,386.695,438.306018,1.470825563],
         [350,423.0174,484.57760,1.626099],
         [400,455.75815,524.42785,1.7598250],
         [500,518.8540,735.5804,2.4683909],
         [600,580.3953,803.62085,2.6967142],
         [800,700.89346,908.089094,3.047278838926],
         [1000,822.35672,1006.81347,3.3785686]]

TimeRef_total=937.225678 # temps problème non réduit 
TimeRef_OneFreq=3.1450526107382 # temps un pas de fréquence
NbFreqSteps=len(frf_no_reduc[0])    # 

NbFluidModes=[]
LastFreq=[]
CPUtime=[]
CPUtimePerStep=[]
frf_reduc_gradient=[]
MeandBdifferenceFRF=[]
MeandBdifferenceGradientPosX=[]
MeandBdifferenceGradientPosY=[]
MeandBdifferenceGradientPosZ=[]
MeandBdifferenceGradientR=[]

for ii in Results:
    NbFluidModes.append(ii[0])
    LastFreq.append(ii[1])
    CPUtime.append(ii[2])
    CPUtimePerStep.append(ii[3])
    nb_mode_F=ii[0]
    results_file='results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient_1.00E+00_1.00E+00_NbModesFluid_'+str(nb_mode_F)+'_5.00E-01_NbModesFluid_'+str(nb_mode_F)+'_8.00E-01_NbModesFluid_'+str(nb_mode_F)+'_results.frf'

    f=open(results_file,'rb')
    frf_reduc_gradient.append(pickle.load(f))
    f.close()

for i in range(len(NbFluidModes)):
    MeandBdifferenceFRF.append(abs(sum(frf_no_reduc[1]-frf_reduc_gradient[i][1]))/NbFreqSteps)
    MeandBdifferenceGradientPosX.append(abs(sum(frf_no_reduc[2]-frf_reduc_gradient[i][2]))/NbFreqSteps)
    MeandBdifferenceGradientPosY.append(abs(sum(frf_no_reduc[3]-frf_reduc_gradient[i][3]))/NbFreqSteps)
    MeandBdifferenceGradientPosZ.append(abs(sum(frf_no_reduc[4]-frf_reduc_gradient[i][4]))/NbFreqSteps)
    MeandBdifferenceGradientR.append(abs(sum(frf_no_reduc[5]-frf_reduc_gradient[i][5]))/NbFreqSteps)

pl.figure(1)
pl.plot(NbFluidModes,MeandBdifferenceFRF,'-o')
pl.xlabel('Nb Fluid Modes')
pl.ylabel('Mean quadratic pressure difference with reference (dB)')

pl.figure(2)
pl.plot(LastFreq,log10(MeandBdifferenceFRF),'-o')
pl.xlabel('Last eigen frequency in reduced basis')
pl.ylabel('Mean quadratic pressure difference with reference (dB)')

pl.figure(3)
pl.plot(NbFluidModes,CPUtimePerStep,'-o')
pl.xlabel('Nb Fluid Modes')
pl.ylabel('CPU time per time step')

pl.figure(4)
pl.plot(CPUtimePerStep,log10(MeandBdifferenceFRF),'-o')
pl.xlabel('CPU time per time step')
pl.ylabel('Mean quadratic pressure difference with reference (dB)')


pl.figure(5)
pl.plot(LastFreq,MeandBdifferenceGradientPosX,'-o')
pl.xlabel('Last eigen frequency in reduced basis')
pl.ylabel('Mean quadratic pressure gradient difference (Pos. X) with reference (dB/m)')

pl.figure(6)
pl.plot(LastFreq,MeandBdifferenceGradientPosY,'-o')
pl.xlabel('Last eigen frequency in reduced basis')
pl.ylabel('Mean quadratic pressure gradient difference (Pos. Y) with reference (dB/m)')


pl.figure(7)
pl.plot(LastFreq,MeandBdifferenceGradientPosZ,'-o')
pl.xlabel('Last eigen frequency in reduced basis')
pl.ylabel('Mean quadratic pressure gradient difference (Pos. Z) with reference (dB/m)')


pl.figure(8)
pl.plot(LastFreq,MeandBdifferenceGradientR,'-o')
pl.xlabel('Last eigen frequency in reduced basis')
pl.ylabel('Mean quadratic pressure gradient difference (R) with reference (dB/m)')





for i in range(len(NbFluidModes)):
    pl.figure(i+10)
    pl.plot(frf_no_reduc[0],10*log10(frf_no_reduc[1]/prefsquare),'b-',label='no reduction', linewidth=1)
    pl.plot(frf_reduc_gradient[i][0],10*log10(frf_reduc_gradient[i][1]/prefsquare),'r-',label='reduction and gradient, nb modes '+str(NbFluidModes[i]), linewidth=1)
    #pl.axis([1.0, 120.0, 70, 105])
    pl.xlabel('Frequency (Hz)')
    pl.ylabel('Mean quadratic pressure (dB)')
    pl.grid('on')
    pl.legend(loc=4)

pl.show()


pl.figure(2)
pl.plot(frf_no_reduc[0],frf_no_reduc[2],'b-',label='no reduction, gradient pos. x', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[2],'r-',label='reduction, gradient pos. x', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(3)
pl.plot(frf_no_reduc[0],frf_no_reduc[3],'b-',label='no reduction, gradient pos. y', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[3],'r-',label='reduction, gradient pos. y', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(4)
pl.plot(frf_no_reduc[0],frf_no_reduc[4],'b-',label='no reduction, gradient pos. z', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[4],'r-',label='reduction, gradient pos. z', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(5)
pl.plot(frf_no_reduc[0],frf_no_reduc[5],'b-',label='no reduction, gradient radius', linewidth=1)
pl.plot(frf_reduc_gradient[0],frf_reduc_gradient[5],'r-',label='reduction, gradient radius', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient in control volume')
pl.grid('on')
pl.legend(loc=4)

pl.figure(6)
pl.plot(frf_no_reduc[0],frf_no_reduc[2]-frf_reduc_gradient[2],'b-',label='gradient difference: no reduction - reduction, pos. x', linewidth=1)
pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure gradient difference in control volume')
pl.grid('on')
pl.legend(loc=4)

#print(len(frf_no_reduc[0]))
#print(abs(sum(frf_no_reduc[2]-frf_reduc_gradient[2]))/len(frf_no_reduc[0]))


#pl.show()


