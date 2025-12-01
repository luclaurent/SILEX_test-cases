import numpy as np
import string
import time
import scipy
import pandas as pd

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')



## With no stiffener
#filename = Path(__file__).parent / 'cube_ballotement_tet4_results.frf'
#f=open(filename,'rb')
#frf_tet4_no_stiff=pickle.load(f)
#f.close()
#
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_results.frf'
#f=open(filename,'rb')
#frf_tet4=pickle.load(f)
#f.close()
#
##
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_results.frf'
#f=open(filename,'rb')
#frf_tet10=pickle.load(f)
#f.close()

elem = 'tet10'

filename = Path(__file__).parent / ('classic_' + elem) / ('cube_sloshing_with_stiffener_'+elem+'_up_p000_down_p000_h20_results.frf')
f=open(filename,'rb')
frf_tet_classic_up_0_down_0=pickle.load(f)
f.close()


#
filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_p000_down_p000_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_0_down_0=pickle.load(f)
f.close()

filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_p200_down_p200_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_2_down_2=pickle.load(f)
f.close()

filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_m200_down_p200_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_m2_down_2=pickle.load(f)
f.close()

filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_m200_down_m200_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_m2_down_m2=pickle.load(f)
f.close()

filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_p200_down_m200_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_2_down_m2=pickle.load(f)
f.close()

filename = Path(__file__).parent / ('xfem_'+elem) / ('cube_xfem_sloshing_with_stiffener_'+elem+'_up_p200_down_m200_h20_results.frf')
f=open(filename,'rb')
frf_tet_xfem_up_m0118_down_m0200_optimum=pickle.load(f)
f.close()

#
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_results.frf'
#f=open(filename,'rb')
#frf_tet10_xfem=pickle.load(f)
#f.close()
#
#frf_tet4_tab=[]
#frf_tet10_tab=[]
#frf_tet4_xfem_tab=[]
#frf_tet10_xfem_tab=[]
#h_tab=[]
#
######################################
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h20_results.frf'
#f=open(filename,'rb')
#frf_tet4_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h20_results.frf'
#f=open(filename,'rb')
#frf_tet10_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h20_results.frf'
#f=open(filename,'rb')
#frf_tet4_xfem_tab.append(pickle.load(f))
#f.close()
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h20_results.frf'
#f=open(filename,'rb')
##
#frf_tet10_xfem_tab.append(pickle.load(f))
#f.close()
#h_tab.append('h 20')
#######################################
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h25_results.frf'
#f=open(filename,'rb')
#frf_tet4_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h25_results.frf'
#f=open(filename,'rb')
#frf_tet10_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h25_results.frf'
#f=open(filename,'rb')
#frf_tet4_xfem_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h25_results.frf'
#f=open(filename,'rb')
#frf_tet10_xfem_tab.append(pickle.load(f))
#f.close()
#h_tab.append('h 25')
#######################################
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h30_results.frf'
#f=open(filename,'rb')
#frf_tet4_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h30_results.frf'
#f=open(filename,'rb')
#frf_tet10_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h30_results.frf'
#f=open(filename,'rb')
#frf_tet4_xfem_tab.append(pickle.load(f))
#f.close()
##
#filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h30_results.frf'
#f=open(filename,'rb')
#frf_tet10_xfem_tab.append(pickle.load(f))
#f.close()
#h_tab.append('h 30')
######################################


## Analytic frequencies
#lx1 = 1.0;
#ly1 = 0.8;
#lz1 = 0.6;
#fmn=[]
#mm=[]
#nn=[]
#for m in range(3):
#    for n in range(3):
#        kmn=np.pi*np.sqrt((m/lx1)**2+(n/ly1)**2)
#        mm.append(m)
#        nn.append(n)
#        fmn.append(np.sqrt(9.81*kmn*np.tanh(kmn*lz1))/(2*np.pi))
#

# store data in DataFrame
df = pd.DataFrame(columns=["freq_0_0", "FRF_0_0", 
                           "freq_0.2_0.2", "FRF_0.2_0.2",
                           "freq_-0.2_0.2", "FRF_-0.2_0.2",
                           "freq_-0.2_-0.2", "FRF_-0.2_-0.2",
                           "freq_0.2_-0.2", "FRF_0.2_-0.2"])

df['freq_0_0'] = frf_tet_xfem_up_0_down_0[0]
df['FRF_0_0'] = abs(frf_tet_xfem_up_0_down_0[1])  /(1000*9.81)
df['freq_0.2_0.2'] = frf_tet_xfem_up_2_down_2[0]  
df['FRF_0.2_0.2'] = abs(frf_tet_xfem_up_2_down_2[1])  /(1000*9.81)
df['freq_-0.2_0.2'] = frf_tet_xfem_up_m2_down_2[0] 
df['FRF_-0.2_0.2'] = abs(frf_tet_xfem_up_m2_down_2[1]) /(1000*9.81)
df['freq_-0.2_-0.2'] = frf_tet_xfem_up_m2_down_m2[0]
df['FRF_-0.2_-0.2'] = abs(frf_tet_xfem_up_m2_down_m2[1])/(1000*9.81)
df['freq_0.2_-0.2'] = frf_tet_xfem_up_2_down_m2[0]
df['FRF_0.2_-0.2'] = abs(frf_tet_xfem_up_2_down_m2[1]) /(1000*9.81)
df.to_csv('FRF_2_param_aachen.csv')

pl.figure(1)
#pl.plot(frf_tet4_no_stiff[0],10*np.log10(abs(frf_tet4_no_stiff[1])),'g-',label='no stiffener / tet4', linewidth=1)
#pl.plot(frf_tet4[0],10*np.log10(abs(frf_tet4[1])),'m-',label='with stiffener / tet4', linewidth=1)
pl.semilogy(frf_tet_classic_up_0_down_0[0]  ,(abs(frf_tet_classic_up_0_down_0[1])  /(1000*9.81)),'y-',label='REF$x_{up}=0$, $x_{down}=0$', linewidth=2)
pl.semilogy(frf_tet_xfem_up_0_down_0[0]  ,(abs(frf_tet_xfem_up_0_down_0[1])  /(1000*9.81)),'r-',label='$x_{up}=0$, $x_{down}=0$', linewidth=2)
pl.semilogy(frf_tet_xfem_up_2_down_2[0]  ,(abs(frf_tet_xfem_up_2_down_2[1])  /(1000*9.81)),'g-',label='$x_{up}=+0.2$, $x_{down}=+0.2$', linewidth=1)
pl.semilogy(frf_tet_xfem_up_m2_down_2[0] ,(abs(frf_tet_xfem_up_m2_down_2[1]) /(1000*9.81)),'c-',label='$x_{up}=-0.2$, $x_{down}=+0.2$', linewidth=1)
pl.semilogy(frf_tet_xfem_up_m2_down_m2[0],(abs(frf_tet_xfem_up_m2_down_m2[1])/(1000*9.81)),'b-',label='$x_{up}=-0.2$, $x_{down}=-0.2$', linewidth=1)
#pl.plot(frf_tet10[0],10*np.log10(abs(frf_tet10[1])),'b-',label='with stiffener / tet10', linewidth=1)
#pl.plot(frf_tet10_xfem[0],10*np.log10(abs(frf_tet10_xfem[1])),'c-',label='XFEM, with stiffener / tet10', linewidth=1)

pl.semilogy([1.15,1.15] ,[1e-5,1],'k--',label='', linewidth=1)

pl.xlabel('Frequency (Hz)')
pl.ylabel('Point A elevation [m]')
pl.grid('on')
pl.legend(loc=3)
pl.xlim(0.4,1.5)
pl.ylim(1e-5,1)

for i in range(len(frf_tet_xfem_up_2_down_m2[0])):
    print(i,frf_tet_xfem_up_2_down_m2[0][i])

pl.figure(2)
#pl.plot(frf_tet4_no_stiff[0],10*np.log10(abs(frf_tet4_no_stiff[1])),'g-',label='no stiffener / tet4', linewidth=1)
#pl.plot(frf_tet4[0],10*np.log10(abs(frf_tet4[1])),'m-',label='with stiffener / tet4', linewidth=1)
pl.semilogy(frf_tet_xfem_up_0_down_0[0]  ,(abs(frf_tet_xfem_up_0_down_0[1])  /(1000*9.81)),'r-',label='$x_{up}=0$, $x_{down}=0$', linewidth=2)
pl.semilogy(frf_tet_xfem_up_2_down_2[0]  ,(abs(frf_tet_xfem_up_2_down_2[1])  /(1000*9.81)),'g-',label='$x_{up}=+0.2$, $x_{down}=+0.2$', linewidth=1)
pl.semilogy(frf_tet_xfem_up_m2_down_2[0] ,(abs(frf_tet_xfem_up_m2_down_2[1]) /(1000*9.81)),'c-',label='$x_{up}=-0.2$, $x_{down}=+0.2$', linewidth=1)
pl.semilogy(frf_tet_xfem_up_m2_down_m2[0],(abs(frf_tet_xfem_up_m2_down_m2[1])/(1000*9.81)),'b-',label='$x_{up}=-0.2$, $x_{down}=-0.2$', linewidth=1)
pl.semilogy(frf_tet_xfem_up_2_down_m2[0] ,(abs(frf_tet_xfem_up_2_down_m2[1]) /(1000*9.81)),'m-',label='$x_{up}=+0.2$, $x_{down}=-0.2$', linewidth=1)
#pl.plot(frf_tet10[0],10*np.log10(abs(frf_tet10[1])),'b-',label='with stiffener / tet10', linewidth=1)
#pl.plot(frf_tet10_xfem[0],10*np.log10(abs(frf_tet10_xfem[1])),'c-',label='XFEM, with stiffener / tet10', linewidth=1)
pl.semilogy(frf_tet_xfem_up_m0118_down_m0200_optimum[0] ,(abs(frf_tet_xfem_up_m0118_down_m0200_optimum[1]) /(1000*9.81)),'y-',label='$x_{up}=-0.118$, $x_{down}=-0.2$', linewidth=1)

pl.semilogy([1.15,1.15] ,[1e-5,1],'k--',label='', linewidth=1)

pl.xlabel('Frequency (Hz)')
pl.ylabel('Point A elevation [m]')
pl.grid('on')
pl.legend(loc=3)
pl.xlim(0.4,1.5)
pl.ylim(1e-5,1)

for i in range(len(frf_tet_xfem_up_2_down_m2[0])):
    print(i,frf_tet_xfem_up_2_down_m2[0][i])


#pl.savefig('img_FRF_2_param_aachen.png', dpi=200)  


pl.show()