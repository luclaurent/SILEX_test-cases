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
print("SILEX CODE - calcul d'une plaque trouee avec des hexaedres a 8 noeuds")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='diabolo'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_diabolo'

# choose the element type
eltype=5

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elementsV1,IdnodesV1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)
elementsV2,IdnodesV2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)

elements=scipy.vstack([elementsV1,elementsV2])

# read surfaces where to impose boundary conditions
# haut
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3)
# bas
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,4)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,8)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,8)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,8)

# Define material
Young1  = 70.0e9 # alu
nu1     = 0.3
rho1    = 2700.0

c1 = 0.1634e6            # rubber
c2 = -1.198e3            # Yeoh parameters
c3 = 3.781e1             #
nu2     = 0.45
Young2  = 2*c1*2*(1+nu2)
rho2    = 1000.0

# Boundary conditions

# define fixed dof
IdNodesFixed_x=IdnodeS4
IdNodesFixed_y=IdnodeS4
IdNodesFixed_z=IdnodeS4

toc = time.clock()
print("time for the reading data part:",toc-tic)

tic0 = time.clock()

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*ndim
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)


#############################################################################
#      compute stiffness matrix
#############################################################################
tic = time.clock()


Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementsV1,[Young1,nu1,rho1])
K1=scipy.sparse.csc_matrix( (Vk1,(Ik1,Jk1)), shape=(ndof,ndof) )
M1=scipy.sparse.csc_matrix( (Vm1,(Ik1,Jk1)), shape=(ndof,ndof) )
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementsV2,[Young2,nu2,rho2])
K2=scipy.sparse.csc_matrix( (Vk2,(Ik2,Jk2)), shape=(ndof,ndof) )
M2=scipy.sparse.csc_matrix( (Vm2,(Ik2,Jk2)), shape=(ndof,ndof) )

K=K1+K2
M=M1+M2

toc = time.clock()
print("time to compute the stiffness matrix:",toc-tic)

#############################################################################
#       Solve the problem
#############################################################################
tic = time.clock()

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],50,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

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
tic = time.clock()

silex_lib_gmsh.WriteResults2(ResultsFileName+'_structure_modes',nodes,elements,eltype,[[eigen_vector_S_list,'nodal',3,'modes']])



toc = time.clock()
print ("time to write results:",toc-tic)
print ("----- END -----")




