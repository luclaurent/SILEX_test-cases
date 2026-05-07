
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
import silex_lib_make_meshes

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

dataPb['mesh_file'] = Path(__file__).parent /  'geom/cavity12_automatic_air_cavity'
#dataPb['results_file'] = Path(__file__).parent /  'results/cavity12_with_impedance_air_flexible_structure'
results_file_name_base='results/cavity12_with_impedance_air_flexible_structure'


dataPb['flag_write_gmsh_results']=1

dataPb['nb_mode_F']= 100
dataPb['nb_mode_S']= 20
dataPb['freq_ini']   = 10.0
dataPb['freq_end']  = 50.0
dataPb['nb_freq_step_per_proc']=40
dataPb['nproc'] = nproc
dataPb['rank'] = rank

# air
dataPb['celerity']=343.0 # ok
dataPb['rho']=1.21 # ok

# shell structure
material_Struc=[]
# alu
# material_Struc.append(75000.0e6) # E Young
# material_Struc.append(0.33) # nu
# material_Struc.append(5.0e-3) # thickness
# material_Struc.append(2700.0) # rho
# WOOD
material_Struc.append(20000.0e6*1e10) # E Young
material_Struc.append(0.30) # nu
material_Struc.append(20.0e-3*1e10) # thickness
material_Struc.append(600.0*1e10) # rho


dataPb['material_Struc']=material_Struc

# structure damping
dataPb['modal_damping_S']=0.0

# impedance paroi : article Walid-JFD : CMAME 2008 : d=50.0 and k=5e6
#dataPb['d_imp_paroi']= 50.0 # Pa.s/m
#dataPb['k_imp_paroi']= 5.0e6 # Pa/m
dataPb['d_imp_paroi']= 0.0 # Pa.s/m
dataPb['k_imp_paroi']= 5.0e16 # Pa/m

##############################################################
# Make fluid mesh
##############################################################
# acoustic cavity geometry

# Parameters: acoustic cavity
#
#          Points at z=0
#                        
#     9 ______________10
#      |              | 
#      |              | ly2
#     4| . . . . . .  |_______________________3
#      |     lx2     11                       |
#      |                            h5        |
#      |                       35 _______37   | 
#      |                         |control|    | ly1    control vol, points at h5-lz5/2
# .....|.........................|  vol  | h5 | 
#      |                       31|_______|33  |
# ly5  |                             .        |
#      |                             .        |
# ....1|_____________________________.________|2
#      .                  lx1        .
#      .                             .
#      .             lx5             .
#
#          Points at z=lz5
# 
#    12 ______________13
#      |              | 
#      |              | ly2
#     8| . . . . . .  |_______________________7
#      |     lx2     14                       |
#      |                            h5        |
#      |                       36 _______38   | 
#      |                         |control|    | ly1 / control vol, points at h5+lz5/2
# .....|.........................|  vol  | h5 | 
#      |                       32|_______|34  |
# ly5  |                             .        |
#      |                             .        |
# ....5|_____________________________.________|6
#      .                  lx1        .
#      .                             .
#      .             lx5             .
#
# 

dataPb['lx1'] = 7.0 # large length of cavity
dataPb['ly1'] = 4.0 # small length of cavity
dataPb['lz1'] = 2.5 # height of cavity
dataPb['ly2'] = 1.0 #
dataPb['lx2'] = 3.5
dataPb['lx5'] = 6.0
dataPb['ly5'] = 1.0
dataPb['lz5'] = 1.5
dataPb['h5']  = 1.0
dataPb['h'] = dataPb['lx1']/20 #size of elements
dataPb['ElementOrder'] = 1


silex_lib_make_meshes.xfem_fluid_cavity(dataPb)

##############################################################
# Make structure mesh
##############################################################
#           
#                    structure side
#               21                         22
#                 __________________________   .  .  .  .  .
#               -             ^              \               r4 : radius
#              /              |                . 18  .  .  .
#            17|              |                |
#              |              |                |
#              |              |                |
#              |              |                |
#              |              | lz4 : height   |  lz4-r4
#              |              |                |
#              |              |                |
#              |              |                |
#       Z      |              |                |
#       ^      |              |                |
#       |    15|______________\/_______________|16
#       |                 l4 : length
#       |______>x_local    
#   
#   
#       
#    Structure position in cavity
#   
#               ^  /
#              /  / 
#            l4  / angle
#   Y        /  /. . . . . . . 
#   ^       /  / .   
#   |      \/ /  .           ly3
#   |            .             
#   |_____> X. . . . . . . . .
#   .            .
#   .     lx3    .
#

dataPb['lx3']   = 2.0
dataPb['ly3']   = 2.0
dataPb['l4']    = 2.5
dataPb['lz4']   = 2.0
dataPb['r4']    = 0.5
dataPb['deg']   = 30 # angle
dataPb['h_struc'] = dataPb['l4']/30 #size of elements
dataPb['ElementOrder_struc'] = 1

silex_lib_make_meshes.structure(dataPb)


#dataPb['results_file']=Path(__file__).parent /  (results_file_name_base + '_with_CB_and_struc_reduction')
#silex_lib_compute_vibroacX.vibroac_Xfem_flex_struc_CB_and_struc_reduction_impedance_paroi(dataPb)
#
#dataPb['results_file']=Path(__file__).parent /  (results_file_name_base + '_no_CB_with_struc_reduction')
#silex_lib_compute_vibroacX.vibroac_Xfem_flex_struc_reduction_impedance_paroi(dataPb)
#
#dataPb['results_file']=Path(__file__).parent /  (results_file_name_base + '_no_CB_no_struc_reduction')
#silex_lib_compute_vibroacX.vibroac_Xfem_flex_struc_impedance_paroi(dataPb)
#
#dataPb['results_file']=Path(__file__).parent /  (results_file_name_base + '_rigid_struc_no_CB_no_struc_reduction')
#silex_lib_compute_vibroacX.vibroac_Xfem_rigid_struc_impedance_paroi(dataPb)

dataPb['results_file']=Path(__file__).parent /  (results_file_name_base + '_rigid_struc_CB_reduction')
silex_lib_compute_vibroacX.vibroac_Xfem_rigid_struc_CB_reduction_impedance_paroi(dataPb)
