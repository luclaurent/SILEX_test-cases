#######################################################################
#                                                                     #
#                      TRACTION SUR LE DIABOLO                        #
#                                                                     #
#######################################################################

#######################################################################
#      Import libraries
#######################################################################

# Librairies generales

import pylab
import time
import scipy
import scipy.linalg
import scipy.sparse
import copy
import mumps
import pickle

# Librairies pour le calcul parallele

#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#nproc = comm.Get_size()
#rank = comm.Get_rank()

#class comm_mumps_one_proc:
#    rank = 0
#   def py2f(self):
#        return 0

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 10  python3 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py


# Fortran librairies 

import sys
sys.path.append('../../librairies')

import silex_lib_hex8_fortran as silex_lib_elt
import silex_lib_gmsh
import silex_lib_extra

#################################################################################################################
print('Calcul en traction statique du diabolo SUPER simplifie : poutre 12 ddls')
#################################################################################################################
#                                                 MATERIALS                                                     #
#################################################################################################################

flag_write_fields=0

#################################################################################################################
#                                                     MESH                                                      #
#################################################################################################################


# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='diabolo'

# Output result file: define the name of the result file 
ResultsFileName='Results_simplifie_traction'

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',3)

# read surfaces where to impose boundary conditions 
elemS1,IdnodS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3) # Surface sup plaque du haut
elemS2,IdnodS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,4) # Surface inf plaque du bas

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf_S1',nodes,elemS1,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf_S2',nodes,elemS2,3)

# get number of nodes, dof and elements from the mesh
nnodes   = nodes.shape[0]
ndof     = 12

#################################################################################
#               BUILD THE R MATRIX                                              #
#################################################################################

R1 = silex_lib_extra.rigidify_surface(IdnodS1,nodes,[0.0,0.104,0.0])
dofS1 = scipy.hstack([(IdnodS1-1)*3,(IdnodS1-1)*3+1,(IdnodS1-1)*3+2])

R2 = silex_lib_extra.rigidify_surface(IdnodS2,nodes,[0.0,0.0,0.0])
dofS2 = scipy.hstack([(IdnodS2-1)*3,(IdnodS2-1)*3+1,(IdnodS2-1)*3+2])

R = scipy.sparse.construct.bmat( [ [ R1 ,  R2[:,list(range(ndof,ndof+6,1))] ] ] )

#################################################################################################################
#                                            BOUNDARY CONDITIONS                                                #
#################################################################################################################

# Boundary conditions
IdNodesFixed_x=2
IdNodesFixed_y=2
IdNodesFixed_z=2

#################################################################################################################
#                                              LOADING PART                                                     #
#################################################################################################################

# load calculation
F = scipy.zeros(12)
# FACE 1 superieure:
F[0]=100.0 # Newton / x
#F[1]=100.0 # Newton / y
#F[2]=100.0 # Newton / z
#F[3]=1.0 # Newton.metre / MX 
#F[4]=1.0 # Newton.metre / MY 
#F[5]=1.0 # Newton.metre / MZ 

#################################################################################################################
#                                              INITIALIZATION                                                   #
#################################################################################################################

# Global initialization
Q          = scipy.zeros(ndof)
Qprime     = scipy.zeros(ndof+6+6)

#################################################################################################################
#                                                 EXPERT PART                                                   #
#################################################################################################################

# define fixed dof
#Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
#SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

#SolvedDofsPrime = scipy.hstack([SolvedDofs,list(range(ndof,ndof+6))])

#SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,dofS1)

SolvedDofsPrime = scipy.setdiff1d(list(range(ndof+6+6)),dofS1)
SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,dofS2)

Fixed_DofsPrime = list(range(ndof+6,ndof+6+6))

SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,Fixed_DofsPrime)

K1_i,K1_j,K1_v = silex_lib_elt.ktan_dd(nodes,elemV1,Q,flag1,param1)
K1 = scipy.sparse.csc_matrix((K1_v,(K1_i,K1_j)),shape=(ndof,ndof))
K2_i,K2_j,K2_v = silex_lib_elt.ktan_dd(nodes,elemV2,Q,flag2,param2)
K2 = scipy.sparse.csc_matrix((K2_v,(K2_i,K2_j)),shape=(ndof,ndof))           
K = K1 + K2                                                                      

Kprime=scipy.sparse.csc_matrix(R.T*K*R)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
Qprime[SolvedDofsPrime] = mumps.spsolve(Kprime[SolvedDofsPrime,:][:,SolvedDofsPrime],Fprime[SolvedDofsPrime])
toc = time.clock()
print("time to solve the problem:",toc-tic)

Q=R*Qprime
#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,3))
disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
disp[range(nnodes),2]=Q[list(range(2,ndof,3))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',3,'displacement']
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,5,fields_to_write)

toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")




















