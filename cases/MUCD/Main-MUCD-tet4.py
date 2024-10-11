#############################################################################
#      Import libraries
#############################################################################
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import sys
sys.path.append('../../librairies')

import silex_lib_gmsh
import silex_lib_tet4 as silex_lib_elt

#############################################################################
print 'SILEX CODE - calcul du MUCD avec des tet4'
#############################################################################

#############################################################################
#      Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='MUCD-tet4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_MUCD-tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',3)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# read surfaces where to impose boundary conditions
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)


# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)

# aluminium
Young = 75000
nu = 0.3

# Boundary conditions
IdNodesFixed_x=IdnodeS2
IdNodesFixed_y=IdnodeS2
IdNodesFixed_z=IdnodeS2

# compute external forces from pressure
# helice grosso modo du 21'' de diametre (269mm de rayon maxi)
# donc kv de 300 et voltage entre 4s et 6s
#press=0.5*1.2*((300*5*3.7*2*scipy.pi/60.0)*250.0e-3)**2*1e-6  # 1/2 * rho * v^2 avec v = omega R et omega = rpm * 2pi / 60 et rpm = kv * voltage
# calcul pour une traction de 20kg divise par la surface 
#press=(20.0*10.0)/(500.0*35.0)

press = 10 #MPa
F = silex_lib_elt.forceonsurface(nodes,elementsS3,press,[0.0,0.0,0.0])


toc = time.clock()
print "time for the user part:",toc-tic

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*ndim
nelem  = elements.shape[0]
print "Number of nodes:",nnodes
print "Number of elements:",nelem

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
print "time to compute the stiffness matrix / FORTRAN:",toc-tic

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
Q[scipy.ix_(SolvedDofs)] = scipy.sparse.linalg.spsolve(K[scipy.ix_(SolvedDofs,SolvedDofs)],F[scipy.ix_(SolvedDofs)])
toc = time.clock()
print "time to solve the problem:",toc-tic

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.clock()

SigmaElem,SigmaNodes,Epsilon,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)

toc = time.clock()
print "time to compute stress and error:",toc-tic
print "The global error is:",ErrorGlobal
print "Total time for the computational part:",toc-tic0

#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,ndim))
disp[range(nnodes),0]=Q[range(0,ndof,3)]
disp[range(nnodes),1]=Q[range(1,ndof,3)]
disp[range(nnodes),2]=Q[range(2,ndof,3)]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[scipy.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[scipy.ix_(range(nelem),[0])],'elemental',1,'Sigma 11'],
                      [SigmaElem[scipy.ix_(range(nelem),[1])],'elemental',1,'Sigma 22'],
                      [SigmaElem[scipy.ix_(range(nelem),[2])],'elemental',1,'Sigma 33'],
                      [SigmaElem[scipy.ix_(range(nelem),[3])],'elemental',1,'Sigma 12'],
                      [SigmaElem[scipy.ix_(range(nelem),[4])],'elemental',1,'Sigma 23'],
                      [SigmaElem[scipy.ix_(range(nelem),[5])],'elemental',1,'Sigma 13'],
                      [SigmaElem[scipy.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[0])],'nodal',1,'Smooth Sigma 11'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[1])],'nodal',1,'Smooth Sigma 22'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[2])],'nodal',1,'Smooth Sigma 33'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[3])],'nodal',1,'Smooth Sigma 12'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[4])],'nodal',1,'Smooth Sigma 23'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[5])],'nodal',1,'Smooth Sigma 13'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)


toc = time.clock()
print "time to write results:",toc-tic
print "total time:",toc-tic0
print "----- END -----"


