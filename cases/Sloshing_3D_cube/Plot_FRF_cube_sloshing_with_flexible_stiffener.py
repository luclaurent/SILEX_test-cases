import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')

# with rigid baffle
filename = Path(__file__).parent / 'cube_xfem_rigid_tet10_results.frf'
f=open(filename,'rb')
frf_tet10_xfem_rigid=pickle.load(f)
f.close()


# with rigid baffle
filename = Path(__file__).parent / 'cube_classic_rigid_tet10_results.frf'
f=open(filename,'rb')
frf_tet10_classic_rigid=pickle.load(f)
f.close()

# with flexible baffle
filename = Path(__file__).parent / 'cube_xfem_flexible_tet10_results.frf'
f=open(filename,'rb')
frf_tet10_xfem_flex=pickle.load(f)
f.close()

pl.figure(1)
pl.plot(frf_tet10_xfem_rigid[0],10*np.log10(abs(frf_tet10_xfem_rigid[1])),'g-',label='Xfem, rigide baffle / tet10', linewidth=1)
pl.plot(frf_tet10_classic_rigid[0],10*np.log10(abs(frf_tet10_classic_rigid[1])),'b-',label='classic, rigide baffle / tet10', linewidth=1)
pl.plot(frf_tet10_xfem_flex[0],10*np.log10(abs(frf_tet10_xfem_flex[1])),'m-',label='Xfem flexible baffle / tet10', linewidth=1)

pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)

pl.show()


