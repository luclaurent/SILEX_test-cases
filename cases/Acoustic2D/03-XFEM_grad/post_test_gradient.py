import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle

#import mumps

import sys
sys.path.append('../../librairies')
import silex_lib_gmsh
import silex_lib_tri3_acou

from mpi4py import MPI
comm = MPI.COMM_WORLD

nproc=comm.Get_size()
rank = comm.Get_rank()

# mpirun -np 2 python Main_xfem.py

##############################################################
results_file='xfem_1'

f=open(results_file+'_KAF_MAF_6200.pck','rb')
[dKFA_dtheta_6200,dMFA_dtheta_6200,KFA_6200,MFA_6200]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_6201.pck','rb')
[dKFA_dtheta_6201,dMFA_dtheta_6201,KFA_6201,MFA_6201]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_6202.pck','rb')
[dKFA_dtheta_6202,dMFA_dtheta_6202,KFA_6202,MFA_6202]=pickle.load(f)
f.close()



dKFA_dtheta_diff=(KFA_6202-KFA_6200)/(0.6202-0.6200)
dMFA_dtheta_diff=(MFA_6202-MFA_6200)/(0.6202-0.6200)

comp_gradient_num_num_K=abs(dKFA_dtheta_6202-dKFA_dtheta_6200).max()/abs(dKFA_dtheta_6200).max()
comp_gradient_num_num_M=abs(dMFA_dtheta_6202-dMFA_dtheta_6200).max()/abs(dMFA_dtheta_6200).max()

comp_gradient_num_dif_K=abs(dKFA_dtheta_6201-dKFA_dtheta_diff).max()/abs(dKFA_dtheta_6201).max()
comp_gradient_num_dif_M=abs(dMFA_dtheta_6201-dMFA_dtheta_diff).max()/abs(dMFA_dtheta_6201).max()

print(dKFA_dtheta_6201-dKFA_dtheta_diff)
scipy.sparse.find(abs(dKFA_dtheta_6201-dKFA_dtheta_diff)>1e-3)
