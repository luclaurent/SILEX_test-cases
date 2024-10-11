#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

#import cProfile

import silex_lib_tet10 as silex_lib_elt
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une helice avec des tet10")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='propeller-tet10'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_propeller_tet10'

# choose the element type
eltype=11

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,3)
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',9,5)
# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,11)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,9)

# valeurs pvc trouvee ici: http://www.canplast.ch/pvc_02.html
# Define material
#Young  = 3600 #MPa
#nu     = 0.4

# hetre
Young = 14300
nu = 0.4

# Boundary conditions
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS2
IdNodesFixed_z=IdnodeS1

# compute external forces from pressure
# helice grosso modo du 21'' de diametre (269mm de rayon maxi)
# donc kv de 300 et voltage entre 4s et 6s
#press=0.5*1.2*((300*5*3.7*2*scipy.pi/60.0)*250.0e-3)**2*1e-6  # 1/2 * rho * v^2 avec v = omega R et omega = rpm * 2pi / 60 et rpm = kv * voltage
# calcul pour une traction de 20kg divise par la surface 
#press=(20.0*10.0)/(500.0*35.0)

press = 0.01 #MPa
F = silex_lib_elt.forceonsurface(nodes,elementsS3,press,[0.0,0.0,0.0])

toc = time.clock()
print("time for the user part:",toc-tic)

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
tic0 = time.clock()
tic = time.clock()
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.clock()
print("time to solve the problem:",toc-tic)

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.clock()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)

toc = time.clock()
print("time to compute stress and error:",toc-tic)
print("The global error is:",ErrorGlobal)
print("Total time for the computational part:",toc-tic0)

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

