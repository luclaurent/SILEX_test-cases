#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg

import sys
sys.path.append('../../librairies')

import cProfile

import silex_lib_tet4 as silex_lib_elt
#import silex_lib_tet4_python as silex_lib_tet4
import silex_lib_gmsh

# MPI PARALLEL LIB
from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank  = comm.Get_rank()

# to run it:
# mpirun -np 3 python Main-para-support-tet4.py

#############################################################################
print 'SILEX CODE - calcul d un support avec des tet4 sur 3 processeurs'
#############################################################################

#logo='<img src="logo-silex.png">'
#print logo

tic = time.process_time()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='support-para-2'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_para_support_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

elementsV1,IdnodesV1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum1',nodes,elementsV1,4)

elementsV2,IdnodesV2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum2',nodes,elementsV2,4)

# Interface
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,5)
silex_lib_gmsh.WriteResults(ResultsFileName+'Interface',nodes,elementsS5,2)

# read surfaces where to impose boundary conditions
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,4)

IdnodesV1 = np.setdiff1d(IdnodesV1,IdnodeS5)
IdnodesV2 = np.setdiff1d(IdnodesV2,IdnodeS5)


# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,2)

# Define material
Young  = 74000.0
nu     = 0.33

# Boundary conditions
IdNodesFixed_x=IdnodeS4
IdNodesFixed_y=IdnodeS4
IdNodesFixed_z=IdnodeS2

# compute external forces from pressure
press = 5000.0/(10.0*2.0*np.pi*7.0) #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0,1.0,0.0]
F = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)

toc = time.process_time()
print 'Proc ',rank,' time for the user part:',toc-tic

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*ndim
nelem  = elementsV1.shape[0]+elementsV2.shape[0]
nelem1  = elementsV1.shape[0]
nelem2  = elementsV2.shape[0]
print "Number of nodes:",nnodes
print "Number of elements:",nelem
print "Number of elements in V1:",nelem1
print "Number of elements in V2:",nelem2

elements=scipy.vstack([elementsV1,elementsV2])

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = np.setdiff1d(range(ndof),Fixed_Dofs)

Dof1=scipy.hstack([(IdnodesV1-1)*3,(IdnodesV1-1)*3+1,(IdnodesV1-1)*3+2])
Dof2=scipy.hstack([(IdnodesV2-1)*3,(IdnodesV2-1)*3+1,(IdnodesV2-1)*3+2])
Dof3=scipy.hstack([(IdnodeS5-1)*3,(IdnodeS5-1)*3+1,(IdnodeS5-1)*3+2])

SolvedDofs1 = np.setdiff1d(Dof1,Fixed_Dofs)
SolvedDofs2 = np.setdiff1d(Dof2,Fixed_Dofs)
SolvedDofs3 = np.setdiff1d(Dof3,Fixed_Dofs)

# initialize displacement vector
Q=np.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.process_time()
tic = time.process_time()
#print silex_lib_elt.stiffnessmatrix.__doc__
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])
toc = time.process_time()

print "time to compute the stiffness matrix :",toc-tic

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.process_time()

MySolve = scipy.sparse.linalg.factorized(K[np.ix_(SolvedDofs1,SolvedDofs1)]) # Makes LU decomposition
U1sl    = MySolve( F[np.ix_(SolvedDofs1)] )
K11inv_K13=np.zeros((len(SolvedDofs1),len(SolvedDofs3)))
for i in range(len(SolvedDofs3)):
    One_dof                               = [SolvedDofs3[i]]
    #K13_i_column                          = np.zeros(len(SolvedDofs1))
    #K13_i_column[range(len(SolvedDofs1))] =
    K13_i_column                          = K[np.ix_(SolvedDofs1),One_dof]
    #K11inv_K13[range(len(SolvedDofs1)),[i]]= MySolve( K13_i_column )
    tmp = np.zeros(len(SolvedDofs1))
    tmp[np.ix_(range(len(SolvedDofs1)))]=K13_i_column.todense()
    xi = MySolve( tmp )
    K11inv_K13[range(len(SolvedDofs1)),[i]]=xi

MySolve = scipy.sparse.linalg.factorized(K[np.ix_(SolvedDofs2,SolvedDofs2)]) # Makes LU decomposition
U2sl    = MySolve( F[np.ix_(SolvedDofs2)] )
K22inv_K23=np.zeros((len(SolvedDofs2),len(SolvedDofs3)))
for i in range(len(SolvedDofs3)):
    One_dof                               = [SolvedDofs3[i]]
    K23_i_column                          = np.zeros(len(SolvedDofs2))
    K23_i_column[range(len(SolvedDofs2))] = K[np.ix_(SolvedDofs2),One_dof].todense()
    K22inv_K23[range(len(SolvedDofs2)),[i]]= MySolve( K23_i_column )

K11inv_K13=scipy.sparse.csc_matrix(K11inv_K13)
K22inv_K23=scipy.sparse.csc_matrix(K22inv_K23)

#Schur = K[np.ix_(SolvedDofs3,SolvedDofs3)]-K[np.ix_(SolvedDofs3,SolvedDofs1)]*K11inv_K13-K[np.ix_(SolvedDofs3,SolvedDofs2)]*K22inv_K23
Schur = K[np.ix_(SolvedDofs3,SolvedDofs3)]-K[np.ix_(SolvedDofs3,SolvedDofs1)]*K11inv_K13-K[np.ix_(SolvedDofs3,SolvedDofs2)]*K22inv_K23
Schur=np.array(Schur.todense())
#F=scipy.sparse.csc_matrix(F)
#U1sl=scipy.sparse.csc_matrix(U1sl)
#U2sl=scipy.sparse.csc_matrix(U2sl)
tmp=F[np.ix_(SolvedDofs3)]-K[np.ix_(SolvedDofs3,SolvedDofs1)]*U1sl-K[np.ix_(SolvedDofs3,SolvedDofs2)]*U2sl
U3=scipy.linalg.solve(Schur,tmp)
print tmp
#print U3
#comm.send(U3, dest=1, tag=11)
#comm.send(U3, dest=2, tag=11)
#U1=U1sl-scipy.dot(K11inv_K13,U3)
#U2=U2sl-scipy.dot(K22inv_K23,U3)
U1=U1sl-scipy.dot(np.array(K11inv_K13.todense()),U3)
U2=U2sl-scipy.dot(np.array(K22inv_K23.todense()),U3)
Q[np.ix_(SolvedDofs1)] = U1
Q[np.ix_(SolvedDofs2)] = U2
Q[np.ix_(SolvedDofs3)] = U3
   
#Q[np.ix_(SolvedDofs)] = scipy.sparse.linalg.spsolve(K[np.ix_(SolvedDofs,SolvedDofs)],F[np.ix_(SolvedDofs)])

toc = time.process_time()
print "time to solve the problem:",toc-tic

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)
toc = time.process_time()
print "time to compute stres and error:",toc-tic
print "The global error is:",ErrorGlobal
print "Total time for the computational part:",toc-tic0

#############################################################################
#         Write results to gmsh format
#############################################################################

tic = time.process_time()

# displacement written on 3 columns:
disp=np.zeros((nnodes,ndim))
disp[range(nnodes),0]=Q[range(0,ndof,3)]
disp[range(nnodes),1]=Q[range(1,ndof,3)]
disp[range(nnodes),2]=Q[range(2,ndof,3)]

# external load written on 3 columns:
load=np.zeros((nnodes,ndim))
load[range(nnodes),0]=F[range(0,ndof,3)]
load[range(nnodes),1]=F[range(1,ndof,3)]
load[range(nnodes),2]=F[range(2,ndof,3)]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[np.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[np.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaElem[np.ix_(range(nelem),[0])],'elemental',1,'Sigma 11'],
                      [SigmaElem[np.ix_(range(nelem),[1])],'elemental',1,'Sigma 22'],
                      [SigmaElem[np.ix_(range(nelem),[2])],'elemental',1,'Sigma 33'],
                      [SigmaElem[np.ix_(range(nelem),[3])],'elemental',1,'Sigma 23'],
                      [SigmaElem[np.ix_(range(nelem),[4])],'elemental',1,'Sigma 13'],
                      [SigmaElem[np.ix_(range(nelem),[5])],'elemental',1,'Sigma 12'],
                      [SigmaElem[np.ix_(range(nelem),[6])],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[np.ix_(range(nnodes),[0])],'nodal',1,'Smooth Sigma 11'],
                      [SigmaNodes[np.ix_(range(nnodes),[1])],'nodal',1,'Smooth Sigma 22'],
                      [SigmaNodes[np.ix_(range(nnodes),[2])],'nodal',1,'Smooth Sigma 33'],
                      [SigmaNodes[np.ix_(range(nnodes),[3])],'nodal',1,'Smooth Sigma 23'],
                      [SigmaNodes[np.ix_(range(nnodes),[4])],'nodal',1,'Smooth Sigma 13'],
                      [SigmaNodes[np.ix_(range(nnodes),[5])],'nodal',1,'Smooth Sigma 12'],
                      [SigmaNodes[np.ix_(range(nnodes),[6])],'nodal',1,'Smooth Sigma V.M.'],
                      [EpsilonElem[np.ix_(range(nelem),[0])],'elemental',1,'Epsilon 11'],
                      [EpsilonElem[np.ix_(range(nelem),[1])],'elemental',1,'Epsilon 22'],
                      [EpsilonElem[np.ix_(range(nelem),[2])],'elemental',1,'Epsilon 33'],
                      [EpsilonElem[np.ix_(range(nelem),[3])]/2.0,'elemental',1,'Epsilon 23'],
                      [EpsilonElem[np.ix_(range(nelem),[4])]/2.0,'elemental',1,'Epsilon 13'],
                      [EpsilonElem[np.ix_(range(nelem),[5])]/2.0,'elemental',1,'Epsilon 12'],
                      [EpsilonNodes[np.ix_(range(nnodes),[0])],'nodal',1,'Smooth Epsilon 11'],
                      [EpsilonNodes[np.ix_(range(nnodes),[1])],'nodal',1,'Smooth Epsilon 22'],
                      [EpsilonNodes[np.ix_(range(nnodes),[2])],'nodal',1,'Smooth Epsilon 33'],
                      [EpsilonNodes[np.ix_(range(nnodes),[3])]/2.0,'nodal',1,'Smooth Epsilon 23'],
                      [EpsilonNodes[np.ix_(range(nnodes),[4])]/2.0,'nodal',1,'Smooth Epsilon 13'],
                      [EpsilonNodes[np.ix_(range(nnodes),[5])]/2.0,'nodal',1,'Smooth Epsilon 12'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

toc = time.process_time()
print "time to write results:",toc-tic
print "----- END -----"



