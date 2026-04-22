#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
#import mumps

import sys
#sys.path.append('../../librairies')
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib-new')
#import mumps
import gmsh

#import silex_lib_tri3_python as silex_lib_elt
from SILEXlib import silex_lib_tri3_fortran as silex_lib_elt
from SILEXlib import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une plaque trouee avec des tri3")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.process_time()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='lame'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_lame-tri3_2026'

# choose the element type
eltype=2

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)


# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,1)

# Define material
Young  = 200000.0
nu     = 0.3
thickness = 1.0

# Boundary conditions

# define fixed dof
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS1

F=silex_lib_elt.forceonline(nodes,elementsS2,[0.0,10.0,0.0,10.0],[0.0,200.0,100.0,200.0])

toc = time.process_time()
print("time for the reading data part: {}".format(toc-tic))

tic0 = time.process_time()
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
Fixed_Dofs = np.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

# define free dof
SolvedDofs = np.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=np.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic = time.process_time()

Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness])

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

toc = time.process_time()
print("time to compute the stiffness matrix: {}".format(toc-tic))

#############################################################################
#       Solve the problem
#############################################################################

tic = time.process_time()
Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.process_time()
print("time to solve the problem: {}".format(toc-tic))

#############################################################################
#       compute stress, smooth stress, strain and error
#############################################################################
tic = time.process_time()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness],Q)

toc = time.process_time()
print("time to compute stresses: {}".format(toc-tic))
print("The global error is:",ErrorGlobal)


#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.process_time()

# displacement written on 2 columns:
disp=np.zeros((nnodes,2))
disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
disp[range(nnodes),1]=Q[list(range(1,ndof,2))]

load=np.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,2))]
load[range(nnodes),1]=F[list(range(1,ndof,2))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'Displacement'],
                  [load,'nodal',2,'Force'],
                  [SigmaElem[range(nelem),[0]],'elemental',1,'Sigma xx'],
                  [SigmaElem[range(nelem),[1]],'elemental',1,'Sigma yy'],
                  [SigmaElem[range(nelem),[2]],'elemental',1,'Sigma xy'],
                  [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                  [SigmaNodes[range(nnodes),[0]],'nodal',1,'Sigma xx Smooth'],
                  [SigmaNodes[range(nnodes),[1]],'nodal',1,'Sigma yy Smooth'],
                  [SigmaNodes[range(nnodes),[2]],'nodal',1,'Sigma xy Smooth'],
                  [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                  [EpsilonElem[range(nelem),[0]],'elemental',1,'Epsilon xx'],
                  [EpsilonElem[range(nelem),[1]],'elemental',1,'Epsilon yy'],
                  [EpsilonElem[range(nelem),[2]]/2.0,'elemental',1,'Epsilon xy'],
                  [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Epsilon xx Smooth'],
                  [EpsilonNodes[range(nnodes),[1]],'nodal',1,'Epsilon yy Smooth'],
                  [EpsilonNodes[range(nnodes),[2]]/2.0,'nodal',1,'Epsilon xy Smooth'],
                  [ErrorElem,'elemental',1,'Error']
                  ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)


toc = time.process_time()
print("time to write results: {}".format(toc-tic))
print("----- END -----")
