import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')



# With no stiffener
filename = Path(__file__).parent / 'cube_ballotement_tet4_results.frf'
f=open(filename,'rb')
frf_tet4_no_stiff=pickle.load(f)
f.close()

#
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_results.frf'
f=open(filename,'rb')
frf_tet4=pickle.load(f)
f.close()

#
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_results.frf'
f=open(filename,'rb')
frf_tet10=pickle.load(f)
f.close()


#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_results.frf'
f=open(filename,'rb')
frf_tet4_xfem=pickle.load(f)
f.close()

filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_results.frf'
f=open(filename,'rb')
frf_tet10_xfem=pickle.load(f)
f.close()

frf_tet4_tab=[]
frf_tet10_tab=[]
frf_tet4_xfem_tab=[]
frf_tet10_xfem_tab=[]
h_tab=[]

######################################
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h20_results.frf'
f=open(filename,'rb')
frf_tet4_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h20_results.frf'
f=open(filename,'rb')
frf_tet10_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h20_results.frf'
f=open(filename,'rb')
frf_tet4_xfem_tab.append(pickle.load(f))
f.close()
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h20_results.frf'
f=open(filename,'rb')
#
frf_tet10_xfem_tab.append(pickle.load(f))
f.close()
h_tab.append('h 20')
######################################
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h25_results.frf'
f=open(filename,'rb')
frf_tet4_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h25_results.frf'
f=open(filename,'rb')
frf_tet10_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h25_results.frf'
f=open(filename,'rb')
frf_tet4_xfem_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h25_results.frf'
f=open(filename,'rb')
frf_tet10_xfem_tab.append(pickle.load(f))
f.close()
h_tab.append('h 25')
######################################
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h30_results.frf'
f=open(filename,'rb')
frf_tet4_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h30_results.frf'
f=open(filename,'rb')
frf_tet10_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h30_results.frf'
f=open(filename,'rb')
frf_tet4_xfem_tab.append(pickle.load(f))
f.close()
#
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h30_results.frf'
f=open(filename,'rb')
frf_tet10_xfem_tab.append(pickle.load(f))
f.close()
h_tab.append('h 30')
######################################


# Analytic frequencies
lx1 = 1.0;
ly1 = 0.8;
lz1 = 0.6;
fmn=[]
mm=[]
nn=[]
for m in range(3):
    for n in range(3):
        kmn=np.pi*np.sqrt((m/lx1)**2+(n/ly1)**2)
        mm.append(m)
        nn.append(n)
        fmn.append(np.sqrt(9.81*kmn*np.tanh(kmn*lz1))/(2*np.pi))



pl.figure(1)
#pl.plot(frf_tet4_no_stiff[0],10*np.log10(abs(frf_tet4_no_stiff[1])),'g-',label='no stiffener / tet4', linewidth=1)
pl.plot(frf_tet4[0],10*np.log10(abs(frf_tet4[1])),'m-',label='with stiffener / tet4', linewidth=1)
pl.plot(frf_tet4_xfem[0],10*np.log10(abs(frf_tet4_xfem[1])),'r-',label='XFEM, with stiffener / tet4', linewidth=1)
pl.plot(frf_tet10[0],10*np.log10(abs(frf_tet10[1])),'b-',label='with stiffener / tet10', linewidth=1)
pl.plot(frf_tet10_xfem[0],10*np.log10(abs(frf_tet10_xfem[1])),'c-',label='XFEM, with stiffener / tet10', linewidth=1)

for i in range(len(fmn)):
    pl.plot(fmn[i],min(10*np.log10(abs(frf_tet4[1]))),'ro')


pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)

#############################
pl.figure(2)
for i in range(len(h_tab)):
    pl.plot(frf_tet4_tab[i][0],10*np.log10(abs(frf_tet4_tab[i][1])),'m-',label='with stiffener / tet4 '+h_tab[i], linewidth=1)
    pl.plot(frf_tet4_xfem_tab[i][0],10*np.log10(abs(frf_tet4_xfem_tab[i][1])),'r-',label='XFEM, with stiffener / tet4 '+h_tab[i], linewidth=1)
    pl.plot(frf_tet10_tab[i][0],10*np.log10(abs(frf_tet10_tab[i][1])),'b-',label='with stiffener / tet10 '+h_tab[i], linewidth=1)
    pl.plot(frf_tet10_xfem_tab[i][0],10*np.log10(abs(frf_tet10_xfem_tab[i][1])),'c-',label='XFEM, with stiffener / tet10 '+h_tab[i], linewidth=1)

pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)

pl.show()


