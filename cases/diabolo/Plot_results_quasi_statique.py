# -*- coding: utf-8 -*-
#################################################################################################################
#                                                POST-TRAITEMENT                                                #
#################################################################################################################

import scipy
import math
import numpy
import pickle
import pylab

font = {'size'   : 14}
pylab.rc('font', **font)



# Dessin ou non?

plot1 = 1 # traction 
plot2 = 0
plot3 = 0
plot4 = 1 # flexion
plot5 = 1 # torsion
plot6 = 0 # CSMA 2015 : compression + torsion + flexion

D = 10e-2  # diametre du diabolo


##
##
##############################



##f = open('U_compression','r')
##U2 = pickle.load(f)
##f.close()
##U2 = scipy.array(U2)
##
##file='lsdyna_compression'
##data = scipy.loadtxt(file)
##T2d = data[:,0]
##U2d = data[:,1]
##
##f = open('U_cisaillement','r')
##U3 = pickle.load(f)
##f.close()
##U3 = scipy.array(U3)
##
##file='lsdyna_cisaillement'
##data = scipy.loadtxt(file)
##T3d = data[:,0]
##U3d = data[:,1]

if plot1 == 1: # traction-compression

    f = open('U_traction.pkl','rb')
    Utraction,Ftraction = pickle.load(f)
    f.close()
    f = open('U_compression.pkl','rb')
    Ucompression,Fcompression = pickle.load(f)
    f.close()

    file='lsdyna_Uy_traction.txt'
    data = scipy.loadtxt(file)
    T1d = data[:,0]
    U1d = data[:,1]

    file='lsdyna_Uy_compression.txt'
    data = scipy.loadtxt(file)
    T2d = data[:,0]
    U2d = data[:,1]

    fig = pylab.figure(1)
    pylab.plot(Utraction,Ftraction,color="blue",label="SPARK")
    pylab.plot(Ucompression,Fcompression,color="blue")
    pylab.plot(U1d,T1d*7725.4248593736938*.5,color="red",label="LS-DYNA")
    pylab.plot(U2d,-T2d*7725.4248593736938,color="red")
    pylab.ylabel('Force de traction [N]')
    pylab.xlabel('Deplacement de la face superieure [m]')
    pylab.grid('on')
    pylab.legend(loc='lower right')
    pylab.title('Traction-compression')
    #pylab.ylim([min(signal)*1.1,max(signal)*1.1])
    #fig.savefig('plot_traction.png')


if plot4 == 1: # Flexion
    
    f = open('U_flexion.pkl','rb')
    Ux1f,Uy1f,Uz1f,Ux2f,Uy2f,Uz2f,couple = pickle.load(f)
    f.close()

    file='lsdyna_Ux1_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Ux1f = data[:,1]

    file='lsdyna_Uy1_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Uy1f = data[:,1]

    file='lsdyna_Uz1_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Uz1f = data[:,1]

    file='lsdyna_Ux2_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Ux2f = data[:,1]

    file='lsdyna_Uy2_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Uy2f = data[:,1]

    file='lsdyna_Uz2_flexion.txt'
    data = scipy.loadtxt(file)
    LS_Uz2f = data[:,1]
    Flsdyna = data[:,0]

    theta = []
    LS_theta = []

    for i in range(len(Ux1f)):
        theta.append(math.atan((Uy2f[i]-Uy1f[i])/(D+Ux2f[i]-Ux1f[i]))*180/math.pi)

    for i in range(len(LS_Uz2f)):
        LS_theta.append(math.atan((LS_Uy2f[i]-LS_Uy1f[i])/(D+LS_Ux2f[i]-LS_Ux1f[i]))*180/math.pi)

    couplemax=1e-3

    fig = pylab.figure(4)
    pylab.plot(theta,couple,color="blue",label="SPARK")
    pylab.plot(LS_theta,Flsdyna*couplemax,color="red",label="LS-DYNA")
    pylab.xlabel('ROTATION DE LA FACE SUPERIEURE [degre]')
    pylab.ylabel('Couple de flexion [N.m]')
    pylab.grid('on')
    pylab.title('Flexion')
    #pylab.legend(loc='lower right')
    #pylab.ylim([min(signal)*1.1,max(signal)*1.1])
    #fig.savefig('plot_flexion.png')

if plot5 == 1: # torsion

    f = open('U_torsion.pkl','rb')
    Ux1tor,Uy1tor,Uz1tor,Ux2tor,Uy2tor,Uz2tor,couple = pickle.load(f)
    f.close()

    file='lsdyna_Ux1_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Ux1tor = data[:,1]

    file='lsdyna_Uy1_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Uy1tor = data[:,1]

    file='lsdyna_Uz1_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Uz1tor = data[:,1]

    file='lsdyna_Ux2_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Ux2tor = data[:,1]

    file='lsdyna_Uy2_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Uy2tor = data[:,1]

    file='lsdyna_Uz2_torsion.txt'
    data = scipy.loadtxt(file)
    LS_Uz2tor = data[:,1]
    Flsdyna = data[:,0]

    theta = []
    LS_theta = []

    for i in range(len(Ux1tor)):
        theta.append(math.atan((Uz1tor[i]-Uz2tor[i])/(D-Ux1tor[i]+Ux2tor[i]))*180/math.pi)

    for i in range(len(LS_Uz2tor)):
        LS_theta.append(math.atan((LS_Uz1tor[i]-LS_Uz2tor[i])/(D-LS_Ux1tor[i]+LS_Ux2tor[i]))*180/math.pi)

    couplemax=1e-2

    fig = pylab.figure(5)
    pylab.plot(theta,couple,color="blue",label="SPARK")
    pylab.plot(LS_theta,Flsdyna*couplemax,color="red",label="LS-DYNA")
    pylab.xlabel('ROTATION DE LA FACE SUPERIEURE [degre]')
    pylab.ylabel('Couple de torsion [N.m]')
    pylab.grid('on')
    pylab.title('Torsion')

pylab.show()
