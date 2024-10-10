import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

#import os
import pylab as pl
import pickle

import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt_fortran
import silex_lib_gmsh

import mumps

##from mpi4py import MPI
##comm = MPI.COMM_WORLD
##nproc=comm.Get_size()
##rank = comm.Get_rank()
##
##class comm_mumps_one_proc:
##    rank = 0
##    def py2f(self):
##        return 0
##mycomm=comm_mumps_one_proc()

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

#mesh_file='geom/cavity2'
#results_file='results/cavity2_damping'

import ComputeCavity3
ANGLES=[0, scipy.pi/8, 2*scipy.pi/8, 3*scipy.pi/8, 4*scipy.pi/8, 5*scipy.pi/8, 6*scipy.pi/8, 7*scipy.pi/8, 8*scipy.pi/8]
#ANGLES=scipy.linspace(0.0, scipy.pi, num=20)
frf=[]
for angle in ANGLES:
    frf.append(ComputeCavity3.ComputeFRF([angle,1.5,2.0,3.0]))

f=open('results/Cavity3_frf_results_angle_study.frf','wb')
pickle.dump([frf, ANGLES], f)
f.close()
