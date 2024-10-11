#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

import silex_lib_tet4_fortran as silex_lib_elt
#import silex_lib_tet4_python as silex_lib_elt
import silex_lib_gmsh
#from mpi4py import MPI
import mumps




#  ^ y
#  |   3
#  |   O---
#  |     \ \_
#  |    _|  _O 2
#  |   |   /
#  |  4|  O 1
#  |   | /
#  |   |/            z
#  |_______________>



#############################################################################
print ("SILEX CODE - calcul d'une piece col de cygne - A380")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='Gooseneck'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_gooseneck_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,5)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,4)

# Hole 1:
# 43061.90 3813.64 -2516.34
# 43062.96 3809.73 -2535.93
# 43139.79 3813.66 -2512.13
# 43140.85 3809.75 -2531.71
l1=43139.79-43061.90
d1=2516.34-2512.13
S1=l1*scipy.pi*d1

# Hole 2:
# 43045.43 3852.62 -2380.57
# 43120.69 3849.48 -2370.46
# 43049.11 3846.70 -2409.75
# 43124.37 3843.56 -2399.64
l2=43124.37-43049.11
d2=2409.75-2399.64
S2=l2*scipy.pi*d2

# Hole 2:
# 43100.52 4096.18 -2579.67
# 43140.13 4094.53 -2574.35
# 43103.59 4091.33 -2604.00
# 43143.20 4089.67 -
l3=43140.13-43100.52
d3=2604.00-2598.68
S3=l3*scipy.pi*d3

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)

# Define material
Young  = 72000.0
nu     = 0.3

# Boundary conditions
IdNodesFixed_x=IdnodeS4
IdNodesFixed_y=IdnodeS4
IdNodesFixed_z=IdnodeS4

# compute external forces from pressure : HOLE 1
Load1x = 500 # N
Load1y = 500 # N
Load1z = 500 # N
press1 = scipy.sqrt(Load1x**2+Load1y**2+Load1z**2)/S1 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction1 = [Load1x,Load1y,Load1z]
F1 = silex_lib_elt.forceonsurface(nodes,elementsS1,press1,direction1)

# compute external forces from pressure : HOLE 2
Load2x = 500 # N
Load2y = 500 # N
Load2z = 500 # N
press2 = scipy.sqrt(Load2x**2+Load2y**2+Load2z**2)/S2 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction2 = [Load2x,Load2y,Load2z]
F2 = silex_lib_elt.forceonsurface(nodes,elementsS2,press2,direction2)

# compute external forces from pressure : HOLE 3
Load3x = 500 # N
Load3y = 500 # N
Load3z = 500 # N
press3 = scipy.sqrt(Load3x**2+Load3y**2+Load3z**2)/S3 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction3 = [Load3x,Load3y,Load3z]
F3 = silex_lib_elt.forceonsurface(nodes,elementsS3,press3,direction3)
F=F1+F2+F3

toc = time.clock()
print ("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*ndim
nelem  = elements.shape[0]
print ("Number of nodes:",nnodes)
print ("Number of elements:",nelem)

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
#print silex_lib_elt.stiffnessmatrix.__doc__
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])
toc = time.clock()

K=scipy.sparse.csr_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)
print ("time to compute the stiffness matrix / FORTRAN:",toc-tic)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs], use_umfpack=True)
Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.clock()
print ("time to solve the problem:",toc-tic)

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.clock()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)

toc = time.clock()
print ("time to compute stress and error:",toc-tic)
print ("The global error is:",ErrorGlobal)
print ("Total time for the computational part:",toc-tic0)

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
                      [load,'nodal',ndim,'Force'],
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



