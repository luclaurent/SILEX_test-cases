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

#################################################################################################################
print('Calcul en traction quasi-statique du diabolo')
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


#################################################################################################################
#                                                     MESH                                                      #
#################################################################################################################


# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='diabolo'

# Output result file: define the name of the result file 
ResultsFileName='Results_traction'

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',3)
elemV1,IdnodV1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',5,1) # volume V1 -> les plaques sur et sous la liaison
elemV2,IdnodV2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',5,2) # volume V2 -> la partie en caoutchouc de la liaison

# read surfaces where to impose boundary conditions 
elemS1,IdnodS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3) # Surface sup plaque du haut
elemS2,IdnodS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,4) # Surface inf plaque du bas

# write the surface mesh in a gmsh-format file to verify if its correct
#gmsh_lib.WriteResults(ResultsFileName+'_vol_V1',nodes,elemV1,1)
#gmsh_lib.WriteResults(ResultsFileName+'_vol_V2',nodes,elemV2,2)
#gmsh_lib.WriteResults(ResultsFileName+'_surf_S1',nodes,elemS1,3)
#gmsh_lib.WriteResults(ResultsFileName+'_surf_S2',nodes,elemS2,4)

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

#################################################################################################################
#                                            BOUNDARY CONDITIONS                                                #
#################################################################################################################


# Dof fixed in the x direction
Fixed_Dofs_x = scipy.hstack([(IdnodS2-1)*3])

# Dof fixed in the x direction
Fixed_Dofs_y = scipy.hstack([(IdnodS2-1)*3+1])

# Dof fixed in the x direction
Fixed_Dofs_z = scipy.hstack([(IdnodS2-1)*3+2])

# Free dof
SolvedDofs = scipy.setdiff1d(list(range(ndof)),Fixed_Dofs_x)
SolvedDofs = scipy.setdiff1d(SolvedDofs,Fixed_Dofs_y)
SolvedDofs = scipy.setdiff1d(SolvedDofs,Fixed_Dofs_z)


#################################################################################################################
#                                              LOADING PART                                                     #
#################################################################################################################

# Quasistatic parameters

n = 10                          # Number of increment in the quasistatic part
load = 1.0e6                     # Quasistatic traction loading (Pa)
direction = scipy.array([0,1,0]) # force in direction +y

# Computation of the scale factor
scale  = scipy.linspace(0,1,n)

# load calculation
Force = silex_lib_elt.forceonsurface(nodes,elemS1,load,direction)
NormFext = scipy.linalg.norm(Force)

#################################################################################################################
#                                              INITIALIZATION                                                   #
#################################################################################################################


# Global initialization
Q          = scipy.zeros(ndof)
QQ         = scipy.zeros(ndof)
QQQ        = scipy.zeros(ndof)
niter      = scipy.zeros(n)
Fext       = scipy.zeros(ndof)
Fint1      = scipy.zeros(ndof)
Fint2      = scipy.zeros(ndof)
disp       = scipy.zeros((nnodes,3))
disp_save  = []
load       = scipy.zeros((nnodes,3))
load_save  = []
sigma_save = []
Usave      = []
Fsave      = []



#################################################################################################################
#                                                 EXPERT PART                                                   #
#################################################################################################################

# NEWTON-RAPHSON PARAMETERS
#==========================

# LOOP OVER QUASISTATIC INCREMENT
#================================

print('Increment : ',n)
print('')

tic = time.clock()

for t in range(n):

    if t%int(n/10.0) == 0:
    
        print('increment : ',t)
        print('max iteration : ',max(niter))

    # Force scaling
    Fext = Force*scale[t]

    load[range(nnodes),0]=Fext[list(range(0,ndof,3))]
    load[range(nnodes),1]=Fext[list(range(1,ndof,3))]
    load[range(nnodes),2]=Fext[list(range(2,ndof,3))]
    load_save.append(load.copy())

    # Internal force prediction
    Fint1 = silex_lib_elt.fint_dd(nodes,elemV1,Q,flag1,param1)
    Fint2 = silex_lib_elt.fint_dd(nodes,elemV2,Q,flag2,param2)
    Fint = Fint1 + Fint2

    # Residual
    R = Fext - Fint 
    rr = R[SolvedDofs]

    # Test for residual nullity
    while (scipy.linalg.norm(R[SolvedDofs]) > (1e-3)*NormFext):  
       
        niter[t] = niter[t]+1        
          
        # Divergence detection
        if ( niter[t] >= 10 ) :
            print('Too many iterations ----> Divergence!')
            stop

        # Stifness matrix 
        K1_i,K1_j,K1_v = silex_lib_elt.ktan_dd(nodes,elemV1,Q,flag1,param1)
        K1 = scipy.sparse.csc_matrix((K1_v,(K1_i,K1_j)),shape=(ndof,ndof))
        K2_i,K2_j,K2_v = silex_lib_elt.ktan_dd(nodes,elemV2,Q,flag2,param2)
        K2 = scipy.sparse.csc_matrix((K2_v,(K2_i,K2_j)),shape=(ndof,ndof))           
        K = K1 + K2                                                                      
        kk = K[SolvedDofs,:][:,SolvedDofs]

        # Correction of displacement
        dQ = scipy.zeros(ndof)
        dQ[SolvedDofs] = mumps.spsolve(kk,rr)
        Q = Q + dQ

        # Internal force correction
        Fint1 = silex_lib_elt.fint_dd(nodes,elemV1,Q,flag1,param1)
        Fint2 = silex_lib_elt.fint_dd(nodes,elemV2,Q,flag2,param2)
        Fint = Fint1 + Fint2

        # Residual
        R = Fext - Fint
        rr = R[SolvedDofs]
    
    # Stored quantities

    Usave.append(Q[(nvisu-1)*3+visu_dir])
    #Fsave.append(Fint[(nvisu-1)*3+visu_dir])
    Fsave.append(sum(Fext))

    disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))]
    disp_save.append(disp.copy())

toc = time.clock()
tps = toc-tic

f=open('U_traction.pkl','wb')
pickle.dump([Usave,Fsave], f)
f.close()

print('')
print('End of Expert part')
print('Time for resolution: ',tps)
print('')

#################################################################################################################
#                                                POST-TRAITEMENT                                                #
#################################################################################################################

# Write GMSH files
field_to_write=[[disp_save,'nodal',3,'displacement'],[load_save,'nodal',3,'Force']]
silex_lib_gmsh.WriteResults2(ResultsFileName,nodes,elements,5,field_to_write)
































