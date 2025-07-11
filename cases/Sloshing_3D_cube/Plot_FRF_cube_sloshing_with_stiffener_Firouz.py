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

Firouz_fig_11_a_f1_x = [0.1, 0.14057, 0.17872, 0.21929, 0.25986, 0.29983, 0.33919, 0.37915, 0.41912, 0.45908, 0.49905, 0.53962, 0.57958, 0.61955, 0.65952, 0.69948, 0.73945, 0.77881, 0.81998, 0.85934, 0.89931]
Firouz_fig_11_a_f1_y = [1.69747, 1.69494, 1.68987, 1.68481, 1.67722, 1.66835, 1.66076, 1.64937, 1.63544, 1.62152, 1.60253, 1.58354, 1.55823, 1.53544, 1.50506, 1.47215, 1.43418, 1.39367, 1.34557, 1.29494, 1.23291]
Firouz_fig_11_a_f1 = [Firouz_fig_11_a_f1_x,Firouz_fig_11_a_f1_y]


Firouz_fig_11_a_f2_x = [0.0992, 0.21896, 0.37944, 0.53912, 0.65888, 0.73952, 0.77944, 0.81936, 0.86008, 0.8992]
Firouz_fig_11_a_f2_y = [2.51302, 2.51302, 2.51302, 2.51302, 2.51302, 2.51637, 2.51804, 2.51804, 2.51971, 2.52473]
Firouz_fig_11_a_f2 = [Firouz_fig_11_a_f2_x,Firouz_fig_11_a_f2_y]

Firouz_fig_11_a_f3 = [[0.1,0.9],[2.66,2.66]]




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
    freq3=eigen_freq[3]
    freq4=eigen_freq[4]
    freq5=eigen_freq[5]


    pl.plot(param,freq1,'ob')
    pl.plot(param,freq2,'og')
    pl.plot(param,freq3,'om')
    pl.plot(param,freq4,'ok')
    pl.plot(param,freq5,'oy')


pl.plot(Firouz_fig_11_a_f1[0],np.array(Firouz_fig_11_a_f1[1])/(2*np.pi*(2/9.81)**(0.5)),'b')
pl.plot(Firouz_fig_11_a_f2[0],np.array(Firouz_fig_11_a_f2[1])/(2*np.pi*(2/9.81)**(0.5)),'g')
pl.plot(Firouz_fig_11_a_f3[0],np.array(Firouz_fig_11_a_f3[1])/(2*np.pi*(2/9.81)**(0.5)),'m')

pl.xlabel('Param')
pl.ylabel('Freq')
pl.grid('on')
pl.legend(loc=4)


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

#pl.figure(3)
#pl.style.use('_mpl-gallery')
#fig, ax = pl.subplots(subplot_kw={"projection": "3d"})
#for k in range(len(param_range)):
#    param = param_range[k]
#    k=k+1
#    freq  = frf_tet10_xfem_tab[k][0]
#    level = frf_tet10_xfem_tab[k][1]
#    #print(level)
#    #pl.plot(freq,10*np.log10(abs(level)),'b-', linewidth=1)
#    ax.plot(param, freq, 10*np.log10(abs(level)))
#    #ax.set(xticklabels=['fff'],yticklabels=['ddf'],zticklabels=['gg'])
#pl.xlabel('Param [m]')
#pl.ylabel('freq [Hz]')
##pl.zlabel('Param [m]')
#
#pl.grid('on')

pl.show()


