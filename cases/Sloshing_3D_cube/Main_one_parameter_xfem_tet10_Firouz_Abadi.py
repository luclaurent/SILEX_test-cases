# librairies
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
import pymumps as mumps

import gmsh
import utils as u
import utils_acoustics as ua

import silex_lib_compute_sloshing
import silex_lib_cube_tank_gmsh_geometry

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

# cube geom
lx = 2.0 # b
ly = 1.0 # a
lz = 1.0 # h

# Baffle position and geom
lx_baffle = 1.0 # d
lz_baffle = 0.40 # e
thickness_baffle = 0.002 # only for classic conforming mesh
lx_baffle_shift_down = 0.0
lx_baffle_shift_up = 0.0

#size of elements
h_fluid_elts =  lx/55

#mesh_file_fluid_tet4        =Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet4'
#results_file_tet4           =Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4_h30'
#mesh_file_tet10_classic     =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10'
#results_file_tet10_classic  =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h30'
#mesh_file_tet4_classic      =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4'
#results_file_tet4_classic   =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h30'

#silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,1,mesh_file_fluid_tet4)
#silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,2,mesh_file_tet10_classic)
#silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,1,mesh_file_tet4_classic)

dataPb['freq_ini'] = 0.3
#dataPb['freq_ref'] = 0.5
dataPb['freq_end'] = 1.2
dataPb['nb_freq_step'] = 3

# Imposed acceleration on tank and stiffener surfaces 
dataPb['U_dot_dot_imposed'] = np.array([1.0,0.0,0.0])

# Flags
dataPb['flag_eigen_vectors'] = 1
dataPb['flag_nb_eigen_modes'] = 3
dataPb['flag_FRF'] = 0
dataPb['flag_write_gmsh_results'] = 1

# fluid
dataFluid = dict()
#data['celerity'] = 343.0
dataFluid['rho'] = 1000.0
#data['fluid_damping'] = 1.0


# parallepipedic cavity with curved structure
mesh_file_fluid_tet10       =Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet10_Firouz'
mesh_file_stiffener         =Path(__file__).parent / 'cube_xfem_sloshing_Stiffener_DKT_Firouz'

results_file_tet10          =Path(__file__).parent / 'cube_xfem_tet10_Firouz_e_h55'
results_file_one_parameter  =Path(__file__).parent / 'cube_xfem_Firouz_tet10_e_h55'

#results_file_one_parameter  =Path(__file__).parent / 'cube_classic_Firouz_tet10_e_h55'
#mesh_file_tet10_classic     =Path(__file__).parent / 'cube_sloshing_with_stiffener_classic_tet10_Firouz_e'
#results_file_tet10_classic  =Path(__file__).parent / 'cube_sloshing_with_stiffener_classic_tet10_Firouz_e_h55'
#results_file_tet10          =results_file_tet10_classic

silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,2,mesh_file_fluid_tet10)


# Firouz, Fig 11 a :
#param_min=0.111
#param_max=0.895
#nb_param_steps=10
param_min=0.4001
param_max=0.6001
nb_param_steps=3

results_frf=[]
results_frf.append(np.linspace(param_min,param_max,nb_param_steps))
results_eigen_frequencies=[]
results_eigen_frequencies.append(np.linspace(param_min,param_max,nb_param_steps))

for param in np.linspace(param_min,param_max,nb_param_steps):
    lz_baffle = param # Firouz, Fig 11 a : lz_baffle = param
    #lx_baffle = 2*param # Firouz, Fig 11 b : lx_baffle = param
    print('------------')
    print('Parameter = ',param)
    print('------------')
 
    silex_lib_cube_tank_gmsh_geometry.Stiffener_DKT(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,h_fluid_elts*0.5,1,mesh_file_stiffener)
    silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10_xfem(dataPb,dataFluid,mesh_file_fluid_tet10,mesh_file_stiffener,results_file_tet10)

    #silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,2,mesh_file_tet10_classic)
    #silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10(dataPb,dataFluid,mesh_file_tet10_classic,results_file_tet10_classic)

    f=open(results_file_tet10.as_posix() +'_eigen_frequencies.pck','rb')
    freq_eigv=pickle.load(f)
    f.close()
    results_eigen_frequencies.append(freq_eigv)

    #f=open(results_file_tet10.as_posix() +'_results.frf','rb')
    #frf_tet10=pickle.load(f)
    #f.close()
    #results_frf.append(frf_tet10)

f=open(results_file_one_parameter.as_posix() +'.pkl','wb')
#pickle.dump([results_frf,results_eigen_frequencies], f)
pickle.dump([results_eigen_frequencies], f)
f.close()


