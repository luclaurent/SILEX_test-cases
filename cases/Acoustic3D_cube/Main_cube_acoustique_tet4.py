import string
import time
import numpy
import scipy
import scipy.sparse
import scipy.sparse.linalg

#import os
import pylab as pl
import pickle

import sys
sys.path.append('../../librairies')

import silex_lib_acou_tet4
import silex_lib_gmsh

import mumps

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

# parallepipedic cavity with plane structure
mesh_file='cube_acou_tet4'
results_file='cube_acou_tet4'


freq_ini     = 150.0
freq_end     = 1000.0
nb_freq_step = 1000

#deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=340.0
rho=1.2
#fluid_damping=(1.0+0.002j)
fluid_damping=1.0

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
fluid_elements_S6,IdNodesS6 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,6) # external surface

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

fluid_nnodes5 = IdNodes5.shape[0]
fluid_nnodes6 = IdNodesS6.shape[0]

print ("Number of fluid nodes:",fluid_nnodes)
print ("Number of fluid elements:",fluid_nelem)

#silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
#silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf,JJf,Vffk,Vffm=silex_lib_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=list(range(fluid_ndof))

##################################################################
# Construct the whole system
##################################################################

K=KFF[SolvedDofF,:][:,SolvedDofF]
M=MFF[SolvedDofF,:][:,SolvedDofF]

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )

F0=numpy.zeros((fluid_ndof))
F0[1-1]=1.0

##############################################################
# FRF computation
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("frequency=",freq)

        F=scipy.array(omega**2*F0 , dtype='d')
        #F=scipy.array(omega**2*F0 , dtype='c16')
        #F[SolvedDofF]=-(KFF[SolvedDofF,:][:,1-1]-(omega**2)*MFF[SolvedDofF,:][:,1-1])*(scipy.zeros((len([1])))+1.0)


        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , scipy.array(F.todense() , dtype=complex) )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )

        sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='d') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        

        #press = scipy.zeros((fluid_ndof),dtype=float)
        press = scipy.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        frf.append(silex_lib_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,press))

        if rank==0:
            press_save.append(press.real)

    frfsave=[scipy.array(frequencies),scipy.array(frf)]

    if rank==0:
        silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

    print ("Time at the end of the FRF:",time.ctime())

    Allfrfsave=[frequencies,frf]
    f=open(results_file+'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()

