import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')

######################################
filename = Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h20_results.frf'
f=open(filename,'rb')
frf_tet4_xfem =pickle.load(f)
f.close()


pl.figure(1)
pl.plot(frf_tet4_xfem[0],10*np.log10(abs(frf_tet4_xfem[1])),'r-',label='XFEM, with stiffener / tet4', linewidth=1)


pl.xlabel('Frequency (Hz)')
pl.ylabel('Pressure [Pa] : Upper corner point 8')
pl.grid('on')
pl.legend(loc=4)
pl.show()
