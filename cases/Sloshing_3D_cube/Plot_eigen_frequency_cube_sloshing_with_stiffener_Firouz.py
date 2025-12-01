import numpy as np

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')

#frf_tet10_xfem_tab=[]
lx_tab=[]

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h15.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_e_h15 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h25.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_e_h25 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h35.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_e_h35 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h45.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_e_h45 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h55.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_e_h55 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_e_h25.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_e_h25 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_e_h35.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_e_h35 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_e_h45.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_e_h45 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_e_h55.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_e_h55 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_d_h15.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_d_h15 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_d_h25.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_d_h25 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_d_h35.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_d_h35 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_xfem_Firouz_tet10_d_h45.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_xfem_tab_tet10_d_h45 = tmp[0]
f.close()
                                   
######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_d_h15.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_d_h15 = tmp[0]
f.close()

######################################
filename = Path(__file__).parent / 'cube_classic_Firouz_tet10_d_h25.pkl'
f=open(filename,'rb')
tmp=pickle.load(f)
#frf_tet10_xfem_tab  = tmp[0]
eigen_freq_classic_tab_tet10_d_h25 = tmp[0]
f.close()


nb_param_steps = len(eigen_freq_xfem_tab_tet10_e_h15)-1

Firouz_fig_11_a_f1_x = [0.1, 0.14057, 0.17872, 0.21929, 0.25986, 0.29983, 0.33919, 0.37915, 0.41912, 0.45908, 0.49905, 0.53962, 0.57958, 0.61955, 0.65952, 0.69948, 0.73945, 0.77881, 0.81998, 0.85934, 0.89931]
Firouz_fig_11_a_f1_y = [1.69747, 1.69494, 1.68987, 1.68481, 1.67722, 1.66835, 1.66076, 1.64937, 1.63544, 1.62152, 1.60253, 1.58354, 1.55823, 1.53544, 1.50506, 1.47215, 1.43418, 1.39367, 1.34557, 1.29494, 1.23291]
Firouz_fig_11_a_f1 = [Firouz_fig_11_a_f1_x,Firouz_fig_11_a_f1_y]


Firouz_fig_11_a_f2_x = [0.0992, 0.21896, 0.37944, 0.53912, 0.65888, 0.73952, 0.77944, 0.81936, 0.86008, 0.8992]
Firouz_fig_11_a_f2_y = [2.51302, 2.51302, 2.51302, 2.51302, 2.51302, 2.51637, 2.51804, 2.51804, 2.51971, 2.52473]
Firouz_fig_11_a_f2 = [Firouz_fig_11_a_f2_x,Firouz_fig_11_a_f2_y]

Firouz_fig_11_a_f3 = [[0.1,0.9],[2.66,2.66]]

Firouz_fig_11_b_f1_x=[0.100445434298440,0.150556792873051,0.2,0.250111358574610,0.299554565701559,0.350334075723830,0.401113585746102,0.451224944320712,0.501336302895322,0.550779510022271,0.600222717149220,0.651002227171492,0.701113585746102,0.751224944320712,0.801336302895322,0.850779510022271,0.900890868596882]
Firouz_fig_11_b_f1_y=[1.6985915492957, 1.6882629107981, 1.6798122065727, 1.6704225352112, 1.6600938967136, 1.6516431924882, 1.6460093896713, 1.6431924882629, 1.6403755868544, 1.6431924882629, 1.6450704225352, 1.6516431924882, 1.6600938967136, 1.6685446009389, 1.6788732394366, 1.6892018779342, 1.6957746478873]

##############################################
pl.figure(10)
plottab=eigen_freq_classic_tab_tet10_e_h55
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'sy',label='283147 nodes, 198495 TET10, compatible mesh')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_classic_e__h55=freq1
        nnodes_classic_e_h55=283147
    freq2=eigen_freq[2]
    pl.plot(param,freq1,'sy')
    pl.plot(param,freq2,'sy')

plottab=eigen_freq_xfem_tab_tet10_e_h15
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'ok',label='7206 nodes, 4281 TET10')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_e_h15=freq1
        nnodes_e_h15=7206
    freq2=eigen_freq[2]
    #freq3=eigen_freq[3]
    #freq4=eigen_freq[4]
    #freq5=eigen_freq[5]
    #freq6=eigen_freq[6]
    #freq7=eigen_freq[7]
    #freq8=eigen_freq[8]
    #freq9=eigen_freq[9]
    pl.plot(param,freq1,'ob')
    pl.plot(param,freq2,'og')
    #pl.plot(param,freq3,'og')
    #pl.plot(param,freq4,'om')
    #pl.plot(param,freq5,'or')
    #pl.plot(param,freq6,'oy')
    #pl.plot(param,freq7,'oc')
    #pl.plot(param,freq8,'ok')
    #pl.plot(param,freq9,'^k')

plottab=eigen_freq_xfem_tab_tet10_e_h25
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'^k',label='28896 nodes, 18861 TET10')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_e_h25=freq1
        nnodes_e_h25=28896
    freq2=eigen_freq[2]
    #freq3=eigen_freq[3]
    #freq4=eigen_freq[4]
    #freq5=eigen_freq[5]
    #freq6=eigen_freq[6]
    #freq7=eigen_freq[7]
    #freq8=eigen_freq[8]
    #freq9=eigen_freq[9]
    pl.plot(param,freq1,'^b')
    pl.plot(param,freq2,'^g')
    #pl.plot(param,freq3,'og')
    #pl.plot(param,freq4,'^m')
    #pl.plot(param,freq5,'^r')
    #pl.plot(param,freq6,'^y')
    #pl.plot(param,freq7,'^c')
    #pl.plot(param,freq8,'^k')
    #pl.plot(param,freq9,'^k')

plottab=eigen_freq_classic_tab_tet10_e_h25
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'sk',label='28896 nodes, 18861 TET10, classic')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_classic_e__h25=freq1
        nnodes_classic_e_h25=28896
    freq2=eigen_freq[2]
#    #freq3=eigen_freq[3]
#    #freq4=eigen_freq[4]
#    #freq5=eigen_freq[5]
#    #freq6=eigen_freq[6]
#    #freq7=eigen_freq[7]
#    #freq8=eigen_freq[8]
#    #freq9=eigen_freq[9]
#    pl.plot(param,freq1,'sb')
#    pl.plot(param,freq2,'sg')
#    #pl.plot(param,freq3,'og')
#    #pl.plot(param,freq4,'sm')
#    #pl.plot(param,freq5,'sr')
#    #pl.plot(param,freq6,'sy')
#    #pl.plot(param,freq7,'sc')
#    #pl.plot(param,freq8,'sk')
#    #pl.plot(param,freq9,'^k')


plottab=eigen_freq_classic_tab_tet10_e_h35
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'sk',label='72522 nodes, 49311 TET10, classic')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_classic_e__h35=freq1
        nnodes_classic_e_h35=72522
    freq2=eigen_freq[2]
    pl.plot(param,freq1,'sk')
    pl.plot(param,freq2,'sk')

plottab=eigen_freq_classic_tab_tet10_e_h45
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'sg',label='150712 nodes, 104952 TET10, classic')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_classic_e__h45=freq1
        nnodes_classic_e_h45=150712
    freq2=eigen_freq[2]
    pl.plot(param,freq1,'sg')
    pl.plot(param,freq2,'sg')


plottab=eigen_freq_xfem_tab_tet10_e_h35
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'>k',label='72522 nodes, 49311 TET10')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_e_h35=freq1
        nnodes_e_h35=72522
    freq2=eigen_freq[2]
    #freq3=eigen_freq[3]
    #freq4=eigen_freq[4]
    #freq5=eigen_freq[5]
    #freq6=eigen_freq[6]
    #freq7=eigen_freq[7]
    #freq8=eigen_freq[8]
    #freq9=eigen_freq[9]
    pl.plot(param,freq1,'>b')
    pl.plot(param,freq2,'>g')
    #pl.plot(param,freq4,'>m')
    #pl.plot(param,freq5,'>r')
    #pl.plot(param,freq6,'>y')
    #pl.plot(param,freq7,'>c')
    #pl.plot(param,freq8,'>k')

plottab=eigen_freq_xfem_tab_tet10_e_h45
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'*k',label='150712 nodes, 104952 TET10')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_e_h45=freq1
        nnodes_e_h45=150712
    freq2=eigen_freq[2]
    #freq3=eigen_freq[3]
    #freq4=eigen_freq[4]
    #freq5=eigen_freq[5]
    #freq6=eigen_freq[6]
    #freq7=eigen_freq[7]
    #freq8=eigen_freq[8]
    #freq9=eigen_freq[9]
    pl.plot(param,freq1,'*b')
    pl.plot(param,freq2,'*g')
    #pl.plot(param,freq4,'*m')
    #pl.plot(param,freq5,'*r')
    #pl.plot(param,freq6,'*y')
    #pl.plot(param,freq7,'*c')
    #pl.plot(param,freq8,'*k')

plottab=eigen_freq_xfem_tab_tet10_e_h55
nb_param_steps=len(plottab)-1
pl.plot(plottab[0][0],plottab[1][0],'sk',label='269355 nodes, 190430 TET10')
for k in range(nb_param_steps):
    param=plottab[0][k]
    eigen_freq=plottab[k+1]
    freq1=eigen_freq[1]
    if k==0:
        f_e_h55=freq1
        nnodes_e_h55=269355
    freq2=eigen_freq[2]
    pl.plot(param,freq1,'sb')
    pl.plot(param,freq2,'sg')


pl.plot(Firouz_fig_11_a_f1[0],np.array(Firouz_fig_11_a_f1[1])/(2*np.pi*(2/9.81)**(0.5)),'-ob',label='1st eigenmode')
#pl.plot(Firouz_fig_11_a_f2[0],np.array(Firouz_fig_11_a_f2[1])/(2*np.pi*(2/9.81)**(0.5)),'g',label='2sd eigenmode')
#pl.plot(Firouz_fig_11_a_f3[0],np.array(Firouz_fig_11_a_f3[1])/(2*np.pi*(2/9.81)**(0.5)),'m',label='3th eigenmode')

pl.xlim(0.35,0.65)
pl.ylim(0.4,0.7)

pl.xlabel('$e/h$')
pl.ylabel('Eigen frequency [Hz]')
pl.grid('on')
pl.legend(loc=3)

f_e_xfem=[f_e_h15,f_e_h25,f_e_h35,f_e_h45,f_e_h55]
nnodes_xfem=[nnodes_e_h15,nnodes_e_h25,nnodes_e_h35,nnodes_e_h45,nnodes_e_h55]

f_e_classic=[f_classic_e__h25,f_classic_e__h35,f_classic_e__h45,f_classic_e__h55]
nnodes_classic=[nnodes_classic_e_h25,nnodes_classic_e_h35,nnodes_classic_e_h45,nnodes_classic_e_h55]
            


print(f_e_xfem)
print(nnodes_xfem)
print([nnodes_classic_e_h55,f_classic_e__h55])

pl.figure(12)
pl.plot(nnodes_xfem,f_e_xfem,'ok',label='xfem mesh')
pl.plot(nnodes_classic_e_h55,f_classic_e__h55,'og',label='compatible mesh')
#pl.plot([1,1e6],[1.6/(2*np.pi*(2/9.81)**(0.5)),1.6/(2*np.pi*(2/9.81)**(0.5))],'-',label='ref. Firouz')
pl.plot([1,1e6],[0.578,0.578],'-',label='ref. Firouz')
pl.xscale('log')
pl.xlim(1000,500000)
pl.ylim(0.4,0.6)
pl.xlabel('Number of dof')
pl.ylabel('First eigen frequency [Hz]')
pl.grid('on')

pl.legend(loc=3)

#
###################################################################
#pl.figure(11)
#
#plottab=eigen_freq_xfem_tab_tet10_d_h15
#pl.plot(plottab[0][0],plottab[1][0],'ok',label='xx nodes, xx TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'ob')
#
#
#plottab=eigen_freq_xfem_tab_tet10_d_h25
#nb_param_steps=len(plottab)-1
#pl.plot(plottab[0][0],plottab[1][0],'^k',label='28896 nodes, 18861 TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'^b')
#
#plottab=eigen_freq_xfem_tab_tet10_d_h35
#nb_param_steps=len(plottab)-1
#pl.plot(plottab[0][0],plottab[1][0],'>k',label='xx nodes, xx TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'>b')
#
#plottab=eigen_freq_xfem_tab_tet10_d_h45
#nb_param_steps=len(plottab)-1
#pl.plot(plottab[0][0],plottab[1][0],'*k',label='xx nodes, xx TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'*b')
#
#
#plottab=eigen_freq_classic_tab_tet10_d_h15
#nb_param_steps=len(plottab)-1
#pl.plot(plottab[0][0],plottab[1][0],'or',label='xx nodes, xx TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'or')
#
#plottab=eigen_freq_classic_tab_tet10_d_h25
#nb_param_steps=len(plottab)-1
#pl.plot(plottab[0][0],plottab[1][0],'^r',label='xx nodes, xx TET10')
#for k in range(nb_param_steps):
#    param=plottab[0][k]
#    eigen_freq=plottab[k+1]
#    freq1=eigen_freq[1]
#    pl.plot(param,freq1,'^r')
#
#
#
#pl.plot(Firouz_fig_11_b_f1_x,np.array(Firouz_fig_11_b_f1_y)/(2*np.pi*(2/9.81)**(0.5)),'b',label='1st eigenmode')
#
#pl.xlim(0.45,0.65)
#pl.ylim(0.5,0.7)
#
#pl.xlabel('$d/h$')
#pl.ylabel('Eigen frequency [Hz]')
#pl.grid('on')
#pl.legend(loc=3)
pl.show()


