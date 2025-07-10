import string
import time
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.sparse.construct

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')
import mumps
import gmsh
# import utils as u
# import utils_acoustics as ua

import silex_lib_compute_sloshing

from SILEXlib import silex_lib_acou_tet10 as libF
from SILEXlib import silex_lib_xfem_acou_tet10 as libF_xfem
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_levelset

#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_acou_tri6 as libFreeSurf
from SILEXlib import silex_lib_gmsh

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4 python3.4 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3.4 Main_toto.py
#
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################
dataPb = dict()

# parallepipedic cavity with plane structure
mesh_file_fluid=Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet10_LF'
mesh_file_stiffener=Path(__file__).parent / 'cube_xfem_sloshing_Stiffener_DKT_h'

results_file=Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet10_LF'

dataPb['freq_ini'] = 0.1
dataPb['freq_ref'] = 0.1
dataPb['freq_end'] = 2.0
dataPb['nb_freq_step'] = 100

# Imposed acceleration on tank and stiffener surfaces 
dataPb['U_dot_dot_imposed'] = np.array([1.0,0.0,0.0])
# dataPb['U_dot_dot_imposed'] = np.array([1.0,0.0,1.0])

# Flags
dataPb['flag_eigen_vectors'] = 0
dataPb['flag_FRF'] = 1
dataPb['flag_write_gmsh_results'] = 1

# fluid
dataFluid = dict()
#data['celerity'] = 343.0
dataFluid['rho'] = 1000.0
#data['fluid_damping'] = 1.0

# structure
#dataStructure = dict()
#dataStructure['rho'] = 2700.0
#dataStructure['thickness'] = 5.0e-3
#dataStructure['young'] = 72000.0e6
#dataStructure['nu'] = 0.3

# LOAD on structure
#dataPb['loaddof'] = 0
#dataPb['valload'] = 1

silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10_xfem(dataPb,dataFluid,mesh_file_fluid,mesh_file_stiffener,results_file)
silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4_xfem(dataPb,dataFluid,mesh_file_fluid,mesh_file_stiffener,results_file, convert_from_tet10=True)
