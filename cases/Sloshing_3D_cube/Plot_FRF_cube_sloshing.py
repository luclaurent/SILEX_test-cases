import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path

#
filename = Path(__file__).parent / 'cube_ballotement_tet4_results.frf'
f=open(filename,'rb')
frf_tet4=pickle.load(f)
f.close()


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
pl.plot(frf_tet4[0],10*np.log10(abs(frf_tet4[1])),'g-',label='tet4', linewidth=1)

for i in range(len(fmn)):
    pl.plot(fmn[i],min(10*np.log10(abs(frf_tet4[1]))),'ro')


pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)
pl.show()


