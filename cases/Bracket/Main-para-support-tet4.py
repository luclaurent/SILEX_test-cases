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
# mpirun -np 3 python3.4 Main-para-support-tet4.py

import mumps

#############################################################################
print ("SILEX CODE - calcul d'un support avec des tet4 sur 3 processeurs")
#############################################################################

#logo='<img src="logo-silex.png">'
#print logo

tic = time.clock()
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

IdnodesV1 = scipy.setdiff1d(IdnodesV1,IdnodeS5)
IdnodesV2 = scipy.setdiff1d(IdnodesV2,IdnodeS5)

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
press = 5000.0/(10.0*2.0*scipy.pi*7.0) #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0,1.0,0.0]
F = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)

toc = time.clock()
if rank==0:
    print( 'Proc ',rank,' time for the user part:',toc-tic)

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
if rank==0:
    print( "Number of nodes:",nnodes)
    print( "Number of elements:",nelem)
    print( "Number of elements in V1:",nelem1)
    print( "Number of elements in V2:",nelem2)

elements=scipy.vstack([elementsV1,elementsV2])

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

Dof1=scipy.hstack([(IdnodesV1-1)*3,(IdnodesV1-1)*3+1,(IdnodesV1-1)*3+2])
Dof2=scipy.hstack([(IdnodesV2-1)*3,(IdnodesV2-1)*3+1,(IdnodesV2-1)*3+2])
Dof3=scipy.hstack([(IdnodeS5-1)*3,(IdnodeS5-1)*3+1,(IdnodeS5-1)*3+2])

SolvedDofs1 = scipy.setdiff1d(Dof1,Fixed_Dofs)
SolvedDofs2 = scipy.setdiff1d(Dof2,Fixed_Dofs)
SolvedDofs3 = scipy.setdiff1d(Dof3,Fixed_Dofs)

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
if rank==0:
    print("time to compute the stiffness matrix :",toc-tic)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Solve the problem
#############################################################################

if rank==0:
    tic = time.clock()

if rank==1:
    print("start k1 factorized")
    MySolve = scipy.sparse.linalg.factorized(K[SolvedDofs1,:][:,SolvedDofs1]) # Makes LU decomposition
    print("k1 factorized done")
    U1sl    = MySolve( F[SolvedDofs1] )
    K11inv_K13=scipy.zeros((len(SolvedDofs1),len(SolvedDofs3)))
    for i in range(len(SolvedDofs3)):
        One_dof                               = [SolvedDofs3[i]]
        K13_i_column                          = scipy.zeros(len(SolvedDofs1))
        K13_i_column[list(range(len(SolvedDofs1)))] = K[SolvedDofs1,:][:,One_dof].todense()
        K11inv_K13[range(len(SolvedDofs1)),[i]]= MySolve( K13_i_column )
    comm.send([K11inv_K13,U1sl], dest=0, tag=11)
    print("k1 done")

if rank==2:
    print("start k2 factorized")
    MySolve = scipy.sparse.linalg.factorized(K[SolvedDofs2,:][:,SolvedDofs2]) # Makes LU decomposition
    print("k2 factorized done")
    U2sl    = MySolve( F[SolvedDofs2] )
    K22inv_K23=scipy.zeros((len(SolvedDofs2),len(SolvedDofs3)))
    for i in range(len(SolvedDofs3)):
        One_dof                               = [SolvedDofs3[i]]
        K23_i_column                          = scipy.zeros(len(SolvedDofs2))
        K23_i_column[list(range(len(SolvedDofs2)))] = K[SolvedDofs2,:][:,One_dof].todense()
        K22inv_K23[range(len(SolvedDofs2)),[i]]= MySolve( K23_i_column )
    comm.send([K22inv_K23,U2sl], dest=0, tag=11)
    print("k2 done")
    
if rank==0:
    data1 = comm.recv(source=1, tag=11)
    K11inv_K13 = data1[0]
    U1sl = data1[1]
    data2 = comm.recv(source=2, tag=11)
    K22inv_K23 = data2[0]
    U2sl = data2[1]
    Schur = K[SolvedDofs3,:][:,SolvedDofs3]-K[SolvedDofs3,:][:,SolvedDofs1]*K11inv_K13-K[SolvedDofs3,:][:,SolvedDofs2]*K22inv_K23
    Schur=scipy.array(Schur)
    tmp=F[SolvedDofs3]-K[SolvedDofs3,:][:,SolvedDofs1]*U1sl-K[SolvedDofs3,:][:,SolvedDofs2]*U2sl
    print("start solve U3 on proc 0")
    U3=scipy.linalg.solve(Schur,tmp)
    print("end solve U3 on 0")
    #print U3
    #comm.send(U3, dest=1, tag=11)
    #comm.send(U3, dest=2, tag=11)
    U1=U1sl-scipy.dot(K11inv_K13,U3)
    U2=U2sl-scipy.dot(K22inv_K23,U3)
    Q[SolvedDofs1] = U1
    Q[SolvedDofs2] = U2
    Q[SolvedDofs3] = U3
   
#Q[scipy.ix_(SolvedDofs)] = scipy.sparse.linalg.spsolve(K[scipy.ix_(SolvedDofs,SolvedDofs)],F[scipy.ix_(SolvedDofs)])
#if rank==0:
#    toc = time.clock()
#    print "time to solve the problem:",toc-tic

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
if rank==0:
    SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)
    toc = time.clock()
    print("time to compute stres and error:",toc-tic)
    print("The global error is:",ErrorGlobal)
    print("Total time for the computational part:",toc-tic0)

#############################################################################
#         Write results to gmsh format
#############################################################################
if rank==0:
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



