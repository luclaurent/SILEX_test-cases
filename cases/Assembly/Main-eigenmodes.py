#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

import silex_lib_hex8_fortran as silex_lib_elt
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'un assemblage - Analyse Modale")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='assemblage'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_assemblage-eigenmodes'

# choose the element type
eltype=5

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# volume des plaque en acier dessus et dessous
elementV10,IdnodeV10=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# volume en caoutchouc
elementV11,IdnodeV11=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,11)

# volume structure au dessus
elementV12,IdnodeV12=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,12)

# face du bas: pied 1
elementS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,1)

# face du bas: pieds 2, 3 , 4
elementS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3)


# face du haut
elementS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,2)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1',nodes,elementS1,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2',nodes,elementS2,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3',nodes,elementS3,3)
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol10',nodes,elementV10,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol11',nodes,elementV11,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'_vol12',nodes,elementV12,5)

elements=scipy.vstack([elementV10,elementV11,elementV12])

silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,5)


# Define material

# (Units : Pa, Kg/m^3, s)

##########################################################
# Material 1 : Aluminium plate over/under the rubber part : V10
E1      = 70.0e9                              
nu1     = 0.3
K1      = E1/(3.0*(1.0-2.0*nu1))
Lambda1 = (E1*nu1)/((1.0+nu1)*(1.0-2.0*nu1))
mu1     = E1/(2.0*(1.0+nu1))
rho1    = 2700.0

flag1  = 2                                      # flag for Neo Hooke hyperelastic model
param1 = [mu1,Lambda1,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

##########################################################
# Material 2 : Rubber part of the model : V11
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

##########################################################
# Material 3 : structure above : V12
E3      = 70.0e9                              
nu3     = 0.3
K3      = E3/(3.0*(1.0-2.0*nu3))
Lambda3 = (E3*nu3)/((1.0+nu3)*(1.0-2.0*nu3))
mu3     = E3/(2.0*(1.0+nu3))
rho3    = 2700.0

flag3  = 2                                      # flag for Neo Hooke hyperelastic model
param3 = [mu3,Lambda3,0.0,0.0,0.0,0.0,0.0,0.0]  # Material parameter in a vector

############################################################
## NUMBER OF EIGENMODES

nb_modes = 100

# Boundary conditions
##IdNodesFixed_x=scipy.hstack([IdnodeS1,IdnodeS3])
##IdNodesFixed_y=scipy.hstack([IdnodeS1,IdnodeS3])
##IdNodesFixed_z=scipy.hstack([IdnodeS1,IdnodeS3])
IdNodesFixed_x=scipy.hstack([IdnodeS1,IdnodeS3])
IdNodesFixed_y=scipy.hstack([IdnodeS1,IdnodeS3])
IdNodesFixed_z=scipy.hstack([IdnodeS1,IdnodeS3])

toc = time.clock()
print("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*3
nelem  = elements.shape[0]
print("Total : ")
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)
print("")
print("links : ")
print("Number of elements:",elementV10.shape[0]+elementV11.shape[0])
print("")
print("structure : ")
print("Number of elements:",elementV12.shape[0])
print("")


# define fixed dof
Fixed_Dofs = scipy.hstack([
    (scipy.array(IdNodesFixed_x)-1)*3,
    (scipy.array(IdNodesFixed_y)-1)*3+1,
    (scipy.array(IdNodesFixed_z)-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()

K1_i,K1_j,K1_v = silex_lib_elt.ktan_dd(nodes,elementV10,Q,flag1,param1)
K1 = scipy.sparse.csc_matrix((K1_v,(K1_i,K1_j)),shape=(ndof,ndof))
K2_i,K2_j,K2_v = silex_lib_elt.ktan_dd(nodes,elementV11,Q,flag2,param2)
K2 = scipy.sparse.csc_matrix((K2_v,(K2_i,K2_j)),shape=(ndof,ndof))           
K3_i,K3_j,K3_v = silex_lib_elt.ktan_dd(nodes,elementV12,Q,flag3,param3)
K3 = scipy.sparse.csc_matrix((K3_v,(K3_i,K3_j)),shape=(ndof,ndof))
K = K1 + K2 + K3

Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementV10,[1,1,rho1])
M1=scipy.sparse.csc_matrix( (Vm1,(Ik1,Jk1)), shape=(ndof,ndof) )
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementV11,[1,1,rho2])
M2=scipy.sparse.csc_matrix( (Vm2,(Ik2,Jk2)), shape=(ndof,ndof) )
Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementV12,[1,1,rho3])
M3=scipy.sparse.csc_matrix( (Vm3,(Ik3,Jk3)), shape=(ndof,ndof) )
M = M1 + M2 + M3

toc = time.clock()
print("time to compute the stiffness and mass matrix :",toc-tic)

#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],nb_modes,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))

eigen_vector_S_list=[]
for i in range(eigen_values_S.shape[0]):
    Q=scipy.zeros(ndof)
    Q[SolvedDofs]=eigen_vectors_S[:,i]
    disp=scipy.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))]
    eigen_vector_S_list.append(disp)

toc = time.clock()
print ("structure eigen frequencies : ",freq_eigv_S)
print ("time for computing the structure modes:",toc-tic)
#############################################################################
#         Write results to gmsh format
#############################################################################

silex_lib_gmsh.WriteResults2(ResultsFileName+'_structure_modes',nodes,elements,5,[[eigen_vector_S_list,'nodal',3,'modes']])

print("----- END -----")



