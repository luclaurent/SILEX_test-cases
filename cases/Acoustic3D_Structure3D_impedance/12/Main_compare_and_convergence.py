
###########################################################
# Libraries
###########################################################

import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pickle

import sys
sys.path.append('../../../librairies')
from pathlib import Path

from SILEXlib import silex_lib_xfem_acou_tet4
from SILEXlib import silex_lib_acou_tet4
from SILEXlib import silex_lib_gmsh
from SILEXlib import silex_lib_dkt
from SILEXlib import silex_lib_porous_tet4

#import pymumps as mumps
#import mumps

import silex_lib_compute_vibroacX

import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_fluid_flex_struc_CB_reduction_impedance_paroi_12.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

##############################################################
# Datas
##############################################################
dataPb = dict()

dataPb['mesh_file   '] = Path(__file__).parent /  'geom/cavity12_with_porous_air'
dataPb['results_file'] = Path(__file__).parent /  'results/cavity12_with_impedance_air_flexible_structure_CB_reduction_test'

dataPb['flag_write_gmsh_results']=1

dataPb['nb_mode_F ']= 100
dataPb['nb_mode_S ']= 20
dataPb['freq_ini  ']   = 10.0
dataPb['freq_end   ']  = 200.0
dataPb['nb_freq_step_per_proc']=100
dataPb['nproc'] = nproc


# air
dataPb['celerity']=343.0 # ok
dataPb['rho']=1.21 # ok

# shell structure
material_Struc=[]
material_Struc.append(75000.0e6) # E Young
material_Struc.append(0.33) # nu
material_Struc.append(5.0e-3) # thickness
material_Struc.append(2700.0) # rho

dataPb['material_Struc']=material_Struc

# structure damping
dataPb['modal_damping_S']=0.0

# impedance paroi : article Walid-JFD : CMAME 2008
dataPb['d_imp_paroi']= 50.0 # Pa.s/m
dataPb['k_imp_paroi']= 5.0e6 # Pa/m


silex_lib_compute_vibroacX.vibroac_Xfem_flex_struc_CB_reduction_impedance_paroi(dataPb)

