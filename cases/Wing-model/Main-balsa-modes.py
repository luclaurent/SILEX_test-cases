#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

import silex_lib_dkt_fortran as silex_lib_elt
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une aile d'avion")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='aile-maquette-balsa'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_aile-maquette-balsa'

# choose the element type
eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# longeron bord de fuite
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# longeron central
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)

# longeron bord d'attaque
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)

# nervures sauf emplanture
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# nervure emplanture
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,5)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf4',nodes,elementsS4,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf5',nodes,elementsS5,2)

elements=scipy.vstack([elementsS1,elementsS2,elementsS3,elementsS4,elementsS5])

silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,2)

# Define material
##Young     = 5140e6 #Pa
##nu        = 0.3
##rho       = 140.0 # kg/m3
Young     = 70000e6 #Pa
nu        = 0.3
rho       = 2800.0 # kg/m3

##thickness1 = 2e-3
##thickness2 = 4e-3
##thickness3 = 2e-3
##thickness4 = 2e-3
##thickness5 = 2e-3

thickness1 = 1e-3
thickness2 = 1e-3
thickness3 = 1e-3
thickness4 = 1e-3
thickness5 = 1e-3


# Boundary conditions
IdNodesFixed_x=IdnodeS5
IdNodesFixed_y=IdnodeS5
IdNodesFixed_z=IdnodeS5
IdNodesFixed_rotx=IdnodeS5
IdNodesFixed_roty=IdnodeS5
IdNodesFixed_rotz=IdnodeS5


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
ndof   = nnodes*6
nelem  = elementsS1.shape[0]+elementsS2.shape[0]+elementsS3.shape[0]+elementsS4.shape[0]+elementsS5.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

#      CLEAN MESH
Id_nodes_used=scipy.unique(elements)
Id_nodes_nonused=scipy.setdiff1d(range(1,nnodes),Id_nodes_used)
IdNodesFixed_x=scipy.hstack([IdNodesFixed_x,Id_nodes_nonused])
IdNodesFixed_y=scipy.hstack([IdNodesFixed_y,Id_nodes_nonused])
IdNodesFixed_z=scipy.hstack([IdNodesFixed_z,Id_nodes_nonused])
IdNodesFixed_rotx=scipy.hstack([IdNodesFixed_rotx,Id_nodes_nonused])
IdNodesFixed_roty=scipy.hstack([IdNodesFixed_roty,Id_nodes_nonused])
IdNodesFixed_rotz=scipy.hstack([IdNodesFixed_rotz,Id_nodes_nonused])

# define fixed dof
Fixed_Dofs = scipy.hstack([
    (scipy.array(IdNodesFixed_x)-1)*6,
    (scipy.array(IdNodesFixed_y)-1)*6+1,
    (scipy.array(IdNodesFixed_z)-1)*6+2,
    (scipy.array(IdNodesFixed_rotx)-1)*6+3,
    (scipy.array(IdNodesFixed_roty)-1)*6+4,
    (scipy.array(IdNodesFixed_rotz)-1)*6+5])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
#print (silex_lib_elt.stiffnessmatrix.__doc__)
Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementsS1,[Young,nu,thickness1,rho])
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementsS2,[Young,nu,thickness2,rho])
Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementsS3,[Young,nu,thickness3,rho])
Ik4,Jk4,Vk4,Vm4=silex_lib_elt.stiffnessmatrix(nodes,elementsS4,[Young,nu,thickness4,rho])
Ik5,Jk5,Vk5,Vm5=silex_lib_elt.stiffnessmatrix(nodes,elementsS5,[Young,nu,thickness5,rho])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K1=scipy.sparse.csc_matrix( (Vk1,(Ik1,Jk1)), shape=(ndof,ndof) ,dtype=float)
K2=scipy.sparse.csc_matrix( (Vk2,(Ik2,Jk2)), shape=(ndof,ndof) ,dtype=float)
K3=scipy.sparse.csc_matrix( (Vk3,(Ik3,Jk3)), shape=(ndof,ndof) ,dtype=float)
K4=scipy.sparse.csc_matrix( (Vk4,(Ik4,Jk4)), shape=(ndof,ndof) ,dtype=float)
K5=scipy.sparse.csc_matrix( (Vk5,(Ik5,Jk5)), shape=(ndof,ndof) ,dtype=float)

K=K1+K2+K3+K4+K5

M1=scipy.sparse.csc_matrix( (Vm1,(Ik1,Jk1)), shape=(ndof,ndof) ,dtype=float)
M2=scipy.sparse.csc_matrix( (Vm2,(Ik2,Jk2)), shape=(ndof,ndof) ,dtype=float)
M3=scipy.sparse.csc_matrix( (Vm3,(Ik3,Jk3)), shape=(ndof,ndof) ,dtype=float)
M4=scipy.sparse.csc_matrix( (Vm4,(Ik4,Jk4)), shape=(ndof,ndof) ,dtype=float)
M5=scipy.sparse.csc_matrix( (Vm5,(Ik5,Jk5)), shape=(ndof,ndof) ,dtype=float)

M=M1+M2+M3+M4+M5


#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],10,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))

eigen_vector_S_list=[]
for i in range(eigen_values_S.shape[0]):
    Q=scipy.zeros(ndof)
    Q[SolvedDofs]=eigen_vectors_S[:,i]
    disp=scipy.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,6))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,6))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,6))]
    eigen_vector_S_list.append(disp)

toc = time.clock()
print ("structure eigen frequencies : ",freq_eigv_S)
print ("time for computing the structure modes:",toc-tic)
#############################################################################
#         Write results to gmsh format
#############################################################################

silex_lib_gmsh.WriteResults2(ResultsFileName+'_structure_modes',nodes,elements,2,[[eigen_vector_S_list,'nodal',3,'modes']])

print("----- END -----")



