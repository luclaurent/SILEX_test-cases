#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

import silex_lib_hex8 as silex_lib_elt
import silex_lib_gmsh
#############################################################################
print("SILEX CODE - calcul d'une plaque trouee avec des hexaedres a 8 noeuds")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='plaque-trouee-hex8'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_plaque-trouee-hex8'

# choose the element type
eltype=5

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,5)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',3,4)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,8)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,8)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,8)

# Define material
Young  = 200000.0
nu     = 0.3

# Boundary conditions

# define fixed dof
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS2
IdNodesFixed_z=IdnodeS3

F = silex_lib_elt.forceonsurface(nodes,elementsS4,10.0,[0.0,1.0,0.0])

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


Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

toc = time.clock()
print("time to compute the stiffness matrix:",toc-tic)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.clock()
print("time to solve the problem:",toc-tic)

#############################################################################
#       compute stress in elements
#############################################################################

tic = time.clock()
SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)
toc = time.clock()
print("time to compute stresses:",toc-tic)
print("The global error is:",ErrorGlobal)


#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,ndim))
disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
disp[range(nnodes),2]=Q[list(range(2,ndof,3))]

# external load written on 3 columns:
load=scipy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,3))]
load[range(nnodes),1]=F[list(range(1,ndof,3))]
load[range(nnodes),2]=F[list(range(2,ndof,3))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[range(nelem),[6]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[6]],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaElem[range(nelem),[0]],'elemental',1,'Sigma 11'],
                      [SigmaElem[range(nelem),[1]],'elemental',1,'Sigma 22'],
                      [SigmaElem[range(nelem),[2]],'elemental',1,'Sigma 33'],
                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma 23'],
                      [SigmaElem[range(nelem),[4]],'elemental',1,'Sigma 13'],
                      [SigmaElem[range(nelem),[5]],'elemental',1,'Sigma 12'],
                      [SigmaElem[range(nelem),[6]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[0]],'nodal',1,'Smooth Sigma 11'],
                      [SigmaNodes[range(nnodes),[1]],'nodal',1,'Smooth Sigma 22'],
                      [SigmaNodes[range(nnodes),[2]],'nodal',1,'Smooth Sigma 33'],
                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Smooth Sigma 23'],
                      [SigmaNodes[range(nnodes),[4]],'nodal',1,'Smooth Sigma 13'],
                      [SigmaNodes[range(nnodes),[5]],'nodal',1,'Smooth Sigma 12'],
                      [SigmaNodes[range(nnodes),[6]],'nodal',1,'Smooth Sigma V.M.'],
                      [EpsilonElem[range(nelem),[0]],'elemental',1,'Epsilon 11'],
                      [EpsilonElem[range(nelem),[1]],'elemental',1,'Epsilon 22'],
                      [EpsilonElem[range(nelem),[2]],'elemental',1,'Epsilon 33'],
                      [EpsilonElem[range(nelem),[3]]/2.0,'elemental',1,'Epsilon 23'],
                      [EpsilonElem[range(nelem),[4]]/2.0,'elemental',1,'Epsilon 13'],
                      [EpsilonElem[range(nelem),[5]]/2.0,'elemental',1,'Epsilon 12'],
                      [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Smooth Epsilon 11'],
                      [EpsilonNodes[range(nnodes),[1]],'nodal',1,'Smooth Epsilon 22'],
                      [EpsilonNodes[range(nnodes),[2]],'nodal',1,'Smooth Epsilon 33'],
                      [EpsilonNodes[range(nnodes),[3]]/2.0,'nodal',1,'Smooth Epsilon 23'],
                      [EpsilonNodes[range(nnodes),[4]]/2.0,'nodal',1,'Smooth Epsilon 13'],
                      [EpsilonNodes[range(nnodes),[5]]/2.0,'nodal',1,'Smooth Epsilon 12'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

toc = time.clock()
print ("time to write results:",toc-tic)
print ("----- END -----")




