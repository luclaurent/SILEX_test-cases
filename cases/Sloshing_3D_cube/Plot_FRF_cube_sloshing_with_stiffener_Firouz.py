import numpy as np

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')

frf_tet10_xfem_tab=[]
lx_tab=[]

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab = tmp[1]

f.close()
#

param_range = frf_tet10_xfem_tab[0]


p_mean=[]
nb_param_steps=len(param_range)
for k in range(nb_param_steps):
    k=k+1
    tmp=0
    titi=frf_tet10_xfem_tab[k][0]
    for i in range(len(frf_tet10_xfem_tab[k][0])):
        tmp=tmp+abs(frf_tet10_xfem_tab[k][1][i])
    
    p_mean.append(tmp/len(frf_tet10_xfem_tab[k][0]))

pl.figure(10)
for k in range(nb_param_steps):
    param=eigen_freq_xfem_tab[0][k]
    eigen_freq=eigen_freq_xfem_tab[k+1]
    freq1=eigen_freq[1]
    freq2=eigen_freq[2]

    pl.plot(param,freq1,'ob')
    pl.plot(param,freq2,'og')

pl.xlabel('Param')
pl.ylabel('Freq')
pl.grid('on')
pl.legend(loc=4)
pl.show()
    
STOP

#############################
pl.figure(1)
for k in range(len(param_range)):
    k=k+1
    freq  = frf_tet10_xfem_tab[k][0]
    #print(freq)
    level = frf_tet10_xfem_tab[k][1]
    #print(level)
    pl.plot(freq,10*np.log10(abs(level)),'b-', linewidth=1)

pl.xlabel('Fequency [Hz]')
pl.ylabel('Mean pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
pl.plot(param_range,10*np.log10(p_mean),'m-', linewidth=1)
pl.xlabel('Param [m]')
pl.ylabel('Mean pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)

pl.figure(3)
pl.style.use('_mpl-gallery')
fig, ax = pl.subplots(subplot_kw={"projection": "3d"})
for k in range(len(param_range)):
    param = param_range[k]
    k=k+1
    freq  = frf_tet10_xfem_tab[k][0]
    level = frf_tet10_xfem_tab[k][1]
    #print(level)
    #pl.plot(freq,10*np.log10(abs(level)),'b-', linewidth=1)
    ax.plot(param, freq, 10*np.log10(abs(level)))
    #ax.set(xticklabels=['fff'],yticklabels=['ddf'],zticklabels=['gg'])
pl.xlabel('Param [m]')
pl.ylabel('freq [Hz]')
#pl.zlabel('Param [m]')

pl.grid('on')

pl.show()


