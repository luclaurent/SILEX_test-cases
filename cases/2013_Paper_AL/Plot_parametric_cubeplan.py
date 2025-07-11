from numpy import *
import string
import time
from scipy.sparse import *
from scipy import *
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve,use_solver,minres,eigen
from numpy.linalg import solve, norm
import os
from scipy.linalg import decomp
import pylab as pl
import pickle

f=open('cubeplan-classic_results.frf','r')
frf_classic_finemesh=pickle.load(f)
f.close()

nbcase=11

pl.figure(1)
frf=[]
for p in range(nbcase):
    f=open('cubeplan-xfem_param_results_'+str(p)+'.frf','r')
    frftmp=pickle.load(f)
    frf.append(frftmp)
    f.close()

prefsquare=20e-6*20e-6

x=0.5
position=[]
for p in range(nbcase):
    position.append(x)
    x=x+0.025

frf[10][1][407]=0.5*(frf[10][1][406]+frf[10][1][408])
frf[10][1][370]=0.5*(frf[10][1][369]+frf[10][1][371])
frf[10][1][446]=0.5*(frf[10][1][444]+frf[10][1][448])
frf[10][1][447]=0.5*(frf[10][1][445]+frf[10][1][449])
frf[0][1][27]=0.5*(frf[0][1][26]+frf[0][1][28])

for p in [10 ,6 ,2, 0]:
    lab='a='+str(position[p])
    pl.plot(frf[p][0],10*log10(frf[p][1]/prefsquare),label=lab, linewidth=2)

for i in range(len(frf[p][0])):
#for i in range(5):
    frfmax=frf[0][1][i]
    frfmin=frf[0][1][i]
    for p in range(nbcase):
#        print "frf[p][1][i]",frf[p][1][i]
        if ( frf[p][1][i]>frfmax ):
            frfmax=frf[p][1][i]
        if ( frf[p][1][i]<frfmin ):
            frfmin=frf[p][1][i]

#    print "frfmax",frfmax
#    print "frfmin",frfmin
#    pl.plot([frf[0][0][i],frf[0][0][i]],[log10(frfmin),log10(frfmax)],'k-')


#pl.plot(frf_classic_finemesh[0],10*log10(frf_classic_finemesh[1]/prefsquare),'k--',label='Reference, a=0.75', linewidth=2)
pl.axis([0.0, 600.0, 80, 150])

pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend()


# Analytical frequencies of the large cavity
celerity=340.0
lx=position[6]
ly=0.6
lz=0.4
n=5

freq=[]
I=[]
J=[]
K=[]
for k in range(n):
    for j in range(n):
        for i in range(n):
            freq.append(celerity*(sqrt((i/lx)**2+(j/ly)**2+(k/lz)**2))/2)
            I.append(i)
            J.append(j)
            K.append(k)

print "large cavity"
print sort(freq)



# Analytical frequencies of the small cavity
celerity=340.0
lx=1.0-position[6]
ly=0.6
lz=0.4
n=5

freq=[]
I=[]
J=[]
K=[]
for k in range(n):
    for j in range(n):
        for i in range(n):
            freq.append(celerity*(sqrt((i/lx)**2+(j/ly)**2+(k/lz)**2))/2)
            I.append(i)
            J.append(j)
            K.append(k)

print "small cavity"
print sort(freq)


# Analytical frequencies of a structure: SS-SS-SS-SS
E=70000.0e6
nu=0.27
rho=2700.0
h=4.0e-3
ly = 0.6
lz = 0.4
D=E*h**3/(12*(1-nu**2))

n=5

freq=[]
I=[]
J=[]

for j in range(n):
    for i in range(n):
        omega=sqrt( (D/(rho*h)))*( ((i+1)*pi/ly)**2+((j+1)*pi/lz)**2 )
        freq.append(omega/(2*pi))
        I.append(i+1)
        J.append(j+1)

print "structure"
print sort(freq)




pl.show()
stop

pl.figure(2)
for i in range(len(frf[p][0])):
#for i in range(5):
    frfmax=frf[0][1][i]
    frfmin=frf[0][1][i]
    for p in range(nbcase):
#        print "frf[p][1][i]",frf[p][1][i]
        if ( frf[p][1][i]>frfmax ):
            frfmax=frf[p][1][i]
        if ( frf[p][1][i]<frfmin ):
            frfmin=frf[p][1][i]

#    print "frfmax",frfmax
#    print "frfmin",frfmin

    pl.plot([frf[0][0][i],frf[0][0][i]],[log10(frfmin),log10(frfmax)],'k-', linewidth=2)

pl.axis([10.0, 600.0, -2, 5])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure')
pl.grid('on')
pl.legend()

pl.show()
