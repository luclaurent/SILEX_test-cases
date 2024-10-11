#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg

import sys
sys.path.append('../../librairies')

import silex_lib_tet4 as silex_lib_elt
import silex_lib_gmsh
import pickle


f=open('Results_parametric','r')
[VM,case]=pickle.load(f)
f.close()

VM=scipy.array(VM)
I=scipy.argmax(VM)

print "VM max = ",max(VM)
Fcombi_worst=case[I]

###################################################################
tic = time.clock()
# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='long_support'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_long_support_tet4_WorstCase'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,4)
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,5)
elementsS6,IdnodeS6=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,6)
elementsS7,IdnodeS7=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,7)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf5',nodes,elementsS5,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf6',nodes,elementsS6,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf7',nodes,elementsS7,2)

# Define material
Young  = 74000.0
nu     = 0.33

# Boundary conditions
IdNodesFixed_x=IdnodeS7
IdNodesFixed_y=IdnodeS7
IdNodesFixed_z=IdnodeS6

# compute external forces from pressure
press = 1.0/(10.0*2.0*scipy.pi*7.0) #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [1.0,0.0,0.0]
F1x = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2x = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3x = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4x = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5x = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)
direction = [0.0,1.0,0.0]
F1y = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2y = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3y = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4y = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5y = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)
direction = [0.0,0.0,1.0]
F1z = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2z = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3z = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4z = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5z = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)

#F=1000.0*F1x-2000.0*F2y+3000.0*F5z
Fworst=  Fcombi_worst[0]*F1x+Fcombi_worst[1]*F1y+Fcombi_worst[2]*F1z+Fcombi_worst[3]*F2x+Fcombi_worst[4]*F2y+Fcombi_worst[5]*F2z+Fcombi_worst[6]*F3x+Fcombi_worst[7]*F3y+Fcombi_worst[8]*F3z+Fcombi_worst[9]*F4x+Fcombi_worst[10]*F4y+Fcombi_worst[11]*F4z+Fcombi_worst[12]*F5x+Fcombi_worst[13]*F5y+Fcombi_worst[14]*F5z

F=Fworst

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
#print silex_lib_elt.stiffnessmatrix.__doc__
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

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)

toc = time.clock()
print "time to compute stres and error:",toc-tic
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

# external load written on 3 columns:
load=scipy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[range(0,ndof,3)]
load[range(nnodes),1]=F[range(1,ndof,3)]
load[range(nnodes),2]=F[range(2,ndof,3)]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[scipy.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      [load,'nodal',ndim,'Force'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaElem[scipy.ix_(range(nelem),[0])],'elemental',1,'Sigma 11'],
                      [SigmaElem[scipy.ix_(range(nelem),[1])],'elemental',1,'Sigma 22'],
                      [SigmaElem[scipy.ix_(range(nelem),[2])],'elemental',1,'Sigma 33'],
                      [SigmaElem[scipy.ix_(range(nelem),[3])],'elemental',1,'Sigma 23'],
                      [SigmaElem[scipy.ix_(range(nelem),[4])],'elemental',1,'Sigma 13'],
                      [SigmaElem[scipy.ix_(range(nelem),[5])],'elemental',1,'Sigma 12'],
                      [SigmaElem[scipy.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[0])],'nodal',1,'Smooth Sigma 11'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[1])],'nodal',1,'Smooth Sigma 22'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[2])],'nodal',1,'Smooth Sigma 33'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[3])],'nodal',1,'Smooth Sigma 23'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[4])],'nodal',1,'Smooth Sigma 13'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[5])],'nodal',1,'Smooth Sigma 12'],
                      [SigmaNodes[scipy.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [EpsilonElem[scipy.ix_(range(nelem),[0])],'elemental',1,'Epsilon 11'],
                      [EpsilonElem[scipy.ix_(range(nelem),[1])],'elemental',1,'Epsilon 22'],
                      [EpsilonElem[scipy.ix_(range(nelem),[2])],'elemental',1,'Epsilon 33'],
                      [EpsilonElem[scipy.ix_(range(nelem),[3])]/2.0,'elemental',1,'Epsilon 23'],
                      [EpsilonElem[scipy.ix_(range(nelem),[4])]/2.0,'elemental',1,'Epsilon 13'],
                      [EpsilonElem[scipy.ix_(range(nelem),[5])]/2.0,'elemental',1,'Epsilon 12'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[0])],'nodal',1,'Smooth Epsilon 11'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[1])],'nodal',1,'Smooth Epsilon 22'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[2])],'nodal',1,'Smooth Epsilon 33'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[3])]/2.0,'nodal',1,'Smooth Epsilon 23'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[4])]/2.0,'nodal',1,'Smooth Epsilon 13'],
                      [EpsilonNodes[scipy.ix_(range(nnodes),[5])]/2.0,'nodal',1,'Smooth Epsilon 12'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

toc = time.clock()
print "time to write results:",toc-tic
print "----- END -----"



