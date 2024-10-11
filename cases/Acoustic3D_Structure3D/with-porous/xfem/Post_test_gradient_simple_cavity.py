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
results_file='results/simple_cavity'

f=open(results_file+'_KAF_MAF_a_660.pck','rb')
[dKFA_dtheta_660,dMFA_dtheta_660,dCSA_dtheta_660,KFA_660,MFA_660,CSA_660]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_a_661.pck','rb')
[dKFA_dtheta_661,dMFA_dtheta_661,dCSA_dtheta_661,KFA_661,MFA_661,CSA_661]=pickle.load(f)
f.close()

f=open(results_file+'_KAF_MAF_a_662.pck','rb')
[dKFA_dtheta_662,dMFA_dtheta_662,dCSA_dtheta_662,KFA_662,MFA_662,CSA_662]=pickle.load(f)
f.close()

dKFA_dtheta_diff=(KFA_662-KFA_660)/(0.662-0.660)
dMFA_dtheta_diff=(MFA_662-MFA_660)/(0.662-0.660)
dCSA_dtheta_diff=(CSA_662-CSA_660)/(0.662-0.660)

#print('----------------------------------------------------------\n')
#print(dMFA_dtheta_661.todense()[1])
#print('----------------------------------------------------------\n')
#print(dMFA_dtheta_diff.todense()[1])
#print('----------------------------------------------------------\n')
#print((dMFA_dtheta_661-dMFA_dtheta_diff).todense()[1])


#print(dKFA_dtheta_661-dKFA_dtheta_diff)
#print(scipy.sparse.find(abs(dKFA_dtheta_661-dKFA_dtheta_diff)>1e-3))
#print('----------------------------------------------------------\n')
#print(dKFA_dtheta_661.todense()[0])
#print('----------------------------------------------------------\n')
#print(dKFA_dtheta_diff.todense()[0])
#print('----------------------------------------------------------\n')
#print((dKFA_dtheta_661-dKFA_dtheta_diff).todense())

print('----------------------------------------------------------\n')
print(dCSA_dtheta_661.todense()[3])
print('----------------------------------------------------------\n')
print(dCSA_dtheta_diff.todense()[3])
print('----------------------------------------------------------\n')
print(dCSA_dtheta_661)
print('----------------------------------------------------------\n')
print(dCSA_dtheta_diff)
print('----------------------------------------------------------\n')
print(dCSA_dtheta_661-dCSA_dtheta_diff)
#print(scipy.sparse.find(abs(dMFA_dtheta_661-dMFA_dtheta_diff)>1e-3))
