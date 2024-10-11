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
print("SILEX CODE - calcul d'un cubesat")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='cubesat1'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_cubesat1'

# choose the element type
eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# surface de la structure du cubesat
elementsS10,IdnodeS10=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# surfaces superieures des plots
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf10',nodes,elementsS10,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf4',nodes,elementsS4,2)

# Pour la visu
elementsS20,IdnodeS20=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,20)
silex_lib_gmsh.WriteResults(ResultsFileName+'_pour_la_visu',nodes,elementsS20,2)

#elements=scipy.vstack([elementsS10,elementsS1,elementsS2,elementsS3,elementsS4])
elements=elementsS10

silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,2)

# Define material: alu
Young     = 70000e6 #Pa
nu        = 0.3
rho       = 2800.0 # kg/m3

thickness = 2e-3

# Boundary conditions
IdNodesFixed_x=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])
IdNodesFixed_y=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])
IdNodesFixed_z=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])
IdNodesFixed_rotx=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])
IdNodesFixed_roty=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])
IdNodesFixed_rotz=scipy.hstack([IdnodeS1,IdnodeS2,IdnodeS3,IdnodeS4])


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
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

#      CLEAN MESH
Id_nodes_used=scipy.unique(elements)
Id_nodes_nonused=scipy.setdiff1d(range(1,nnodes),Id_nodes_used)
IdNodesFixed_x=scipy.unique(scipy.hstack([IdNodesFixed_x,Id_nodes_nonused]))
IdNodesFixed_y=scipy.unique(scipy.hstack([IdNodesFixed_y,Id_nodes_nonused]))
IdNodesFixed_z=scipy.unique(scipy.hstack([IdNodesFixed_z,Id_nodes_nonused]))
IdNodesFixed_rotx=scipy.unique(scipy.hstack([IdNodesFixed_rotx,Id_nodes_nonused]))
IdNodesFixed_roty=scipy.unique(scipy.hstack([IdNodesFixed_roty,Id_nodes_nonused]))
IdNodesFixed_rotz=scipy.unique(scipy.hstack([IdNodesFixed_rotz,Id_nodes_nonused]))

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
Ik,Jk,Vk,Vm=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness,rho])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

M=scipy.sparse.csc_matrix( (Vm,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)


#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],10,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')
#eigen_values_S,eigen_vectors_S= scipy.linalg.eig(K[SolvedDofs,:][:,SolvedDofs],M[SolvedDofs,:][:,SolvedDofs])
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



