import string
import time
#import scipy
#import scipy.sparse
#import scipy.sparse.linalg

#import os
#import pylab as pl
import pickle

#import sys
#sys.path.append('../../../librairies')


#import silex_lib_xfem_acou_tet4
#import silex_lib_dkt_fortran
#import silex_lib_gmsh
#import silex_lib_porous_tet4_fortran

#import mumps

#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#nproc=comm.Get_size()
#rank = comm.Get_rank()

#class comm_mumps_one_proc:
#    rank = 0
#    def py2f(self):
#        return 0
#mycomm=comm_mumps_one_proc()

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

import ComputePorousCavity5
#ANGLES=[0,30,60,90,120,150,180,210,240,270,300,330]
ANGLES=scipy.linspace(0.0, 270.0, num=4)
frf=[]
for angle in ANGLES:
    frf.append(ComputePorousCavity5.ComputeFRF([angle]))

f=open('results/CavityPorous5_frf_results_angle_study.frf','wb')
pickle.dump([frf, ANGLES], f)
f.close()
