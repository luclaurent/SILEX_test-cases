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
print('Calcul en traction statique du diabolo simplifie')
#################################################################################################################
#                                                 MATERIALS                                                     #
#################################################################################################################


# (Units : Pa, Kg/m^3, s)

##########################################################
# Material 1 : Aluminium plate over/under the rubber part
E1      = 70.0e9                              
nu1     = 0.3
K1      = E1/(3.0*(1.0-2.0*nu1))
Lambda1 = (E1*nu1)/((1.0+nu1)*(1.0-2.0*nu1))
mu1     = E1/(2.0*(1.0+nu1))
rho1    = 2700.0

flag1  = 2                                      # flag for Neo Hooke hyperelastic model
param1 = [mu1,Lambda1,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

##########################################################
# Material 2 : Rubber part of the model
c1 = 0.1634e6            #
c2 = -1.198e3            # Yeoh parameters
c3 = 3.781e1             #

# for quasi-incompressible rubber nu is in [0.495,0.5]
# mu = E/2(1+nu) -> mu~= E/3
# mu = 2c1 -> E = 6c1

nu2     = 0.45
E2      = 2*c1*2*(1+nu2)
K2      = E2/(3.0*(1.0-2.0*nu2))
Lambda2 = (E2*nu2)/((1.0+nu2)*(1.0-2.0*nu2))
mu2     = E2/(2.0*(1.0+nu2))
rho2    = 1000.0

flag2  = 5                                   # flag for Yeoh hyperelastic model
param2 = [mu2,Lambda2,c1,c2,c3,0.0,0.0,0.0]  # Material parameter in a vector

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
elemV1,IdnodV1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',5,1) # volume V1 -> les plaques sur et sous la liaison
elemV2,IdnodV2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',5,2) # volume V2 -> la partie en caoutchouc de la liaison

# read surfaces where to impose boundary conditions 
elemS1,IdnodS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3) # Surface sup plaque du haut
elemS2,IdnodS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,4) # Surface inf plaque du bas

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol_V1',nodes,elemV1,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol_V2',nodes,elemV2,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf_S1',nodes,elemS1,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf_S2',nodes,elemS2,3)

# get number of nodes, dof and elements from the mesh
nnodes   = nodes.shape[0]
ndof     = nodes.shape[0]*nodes.shape[1]
nelemV1  = elemV1.shape[0]
nelemV2  = elemV2.shape[0]
nelem    = nelemV1 + nelemV2
nGauss   = 8

print ("Number of nodes:",nnodes)
print ("Number of elements:",nelemV1+nelemV2)
print (" ")

# node to be saved, and its direction
nvisu = 13
visu_dir = 1   # x = 0 | y = 1 | z = 2

elements = scipy.vstack([elemV1,elemV2])

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
#IdNodesFixed_x=IdnodS2
#IdNodesFixed_y=IdnodS2
#IdNodesFixed_z=IdnodS2

#################################################################################################################
#                                              LOADING PART                                                     #
#################################################################################################################

# load calculation
Fprime = scipy.zeros(ndof+6+6)
# FACE 1 superieure:
Fprime[ndof+0]=100.0 # Newton / x
#Fprime[ndof+1]=100.0 # Newton / y
#Fprime[ndof+2]=100.0 # Newton / z
#Fprime[ndof+3]=1.0 # Newton.metre / MX 
#Fprime[ndof+4]=1.0 # Newton.metre / MY 
#Fprime[ndof+5]=1.0 # Newton.metre / MZ 

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




















