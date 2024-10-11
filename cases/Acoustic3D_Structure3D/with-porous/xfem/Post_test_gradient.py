import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle

#import mumps

import sys
sys.path.append('../../../librairies')
import silex_lib_gmsh
import silex_lib_tri3_acou

from mpi4py import MPI
comm = MPI.COMM_WORLD

nproc=comm.Get_size()
rank = comm.Get_rank()

# mpirun -np 2 python Main_xfem.py

##############################################################
results_file='results/cavity9_with_porous_air_flexible_structure_CB_ang30_gradient'

f=open(results_file+'_KAF_MAF_lx3_2_000.pck','rb')
[dKFA_dtheta_2000,dMFA_dtheta_2000,dCSA_dtheta_2000,KFA_2000,MFA_2000,CSA_2000]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_lx3_2_001.pck','rb')
[dKFA_dtheta_2001,dMFA_dtheta_2001,dCSA_dtheta_2001,KFA_2001,MFA_2001,CSA_2001]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_lx3_2_002.pck','rb')
[dKFA_dtheta_2002,dMFA_dtheta_2002,dCSA_dtheta_2002,KFA_2002,MFA_2002,CSA_2002]=pickle.load(f)
f.close()

dKFA_dtheta_diff=(KFA_2002-KFA_2000)/(2.002-2.000)
dMFA_dtheta_diff=(MFA_2002-MFA_2000)/(2.002-2.000)
dCSA_dtheta_diff=(CSA_2002-CSA_2000)/(2.002-2.000)

#print(dKFA_dtheta_2001-dKFA_dtheta_diff)
#print(scipy.sparse.find(abs(dKFA_dtheta_2001-dKFA_dtheta_diff)>1e-3))

#print(dMFA_dtheta_2001-dMFA_dtheta_diff)
#print(scipy.sparse.find(abs(dMFA_dtheta_2001-dMFA_dtheta_diff)>1e-3))

print(dCSA_dtheta_2001-dCSA_dtheta_diff)
print(scipy.sparse.find(abs(dMFA_dtheta_2001-dMFA_dtheta_diff)>1e-3))
