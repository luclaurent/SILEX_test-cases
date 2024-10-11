#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

#import silex_lib_tet4_fortran as silex_lib_elt
import silex_lib_tet4_python as silex_lib_elt
import silex_lib_gmsh
import silex_lib_extra_python as silex_lib_extra

#############################################################################
print("SILEX CODE - calcul d'une fourche avec des tet4")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='fourche-tet4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_fourche_rigid_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)

# Define material
Young  = 75000.0
nu     = 0.27

# Boundary conditions
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS1
IdNodesFixed_z=IdnodeS1

# compute external forces from pressure

# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1

#direction2 = [-2000.0+3000,-2000.0,0.0]
#direction3 = [-2000.0-3000,-2000.0,0.0]
#press2 = scipy.sqrt(direction2[0]**2+direction2[1]**2+direction2[2]**2)/(10.0*scipy.pi*10.0)
#press3 = scipy.sqrt(direction3[0]**2+direction3[1]**2+direction3[2]**2)/(10.0*scipy.pi*10.0)
#F2 = silex_lib_elt.forceonsurface(nodes,elementsS2,press2,direction2)
#F3 = silex_lib_elt.forceonsurface(nodes,elementsS3,press3,direction3)
#F=F2+F3

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

#################################################################################
#               BUILD THE R MATRIX                                              #
#################################################################################

R2 = silex_lib_extra.rigidify_surface(IdnodeS2,nodes,[52.0,39.0,36.0])
dofS2 = scipy.hstack([(IdnodeS2-1)*3,(IdnodeS2-1)*3+1,(IdnodeS2-1)*3+2])

R3 = silex_lib_extra.rigidify_surface(IdnodeS3,nodes,[52.0,39.0,-36.0])
dofS3 = scipy.hstack([(IdnodeS3-1)*3,(IdnodeS3-1)*3+1,(IdnodeS3-1)*3+2])

sparse_ones = scipy.sparse.csc_matrix( (list(scipy.ones(ndof)),(list(range(ndof)),list(range(ndof)))), shape=(ndof,ndof) )

R = scipy.sparse.construct.bmat( [ [ sparse_ones
                                     ,R2[:,list(range(ndof,ndof+6,1))]
                                     ,R3[:,list(range(ndof,ndof+6,1))]
                                     ]
                                   ]
                                 )
#################################################################################################################
#                                              LOADING PART                                                     #
#################################################################################################################

# load calculation
Fprime = scipy.zeros(ndof+6+6)
# FACE 2 :
Fprime[ndof+0]=-2000.0+3000 # Newton / x
Fprime[ndof+1]=-2000.0 # Newton / y
Fprime[ndof+2]=0.0 # Newton / z
Fprime[ndof+3]=0.0 # Newton.metre / MX 
Fprime[ndof+4]=0.0 # Newton.metre / MY
Fprime[ndof+5]=0.0 # Newton.metre / MZ 

# FACE 3 :
Fprime[ndof+6]=-2000.0-3000 # Newton / x
Fprime[ndof+7]=-2000.0 # Newton / y
Fprime[ndof+8]=0.0 # Newton / z
Fprime[ndof+9]=0.0 # Newton.metre / MX 
Fprime[ndof+10]=0.0 # Newton.metre / MY 
Fprime[ndof+11]=0.0 # Newton.metre / MZ 
#################################################################################################################
#                                              INITIALIZATION                                                   #
#################################################################################################################

# Global initialization
Q          = scipy.zeros(ndof)
Qprime     = scipy.zeros(ndof+6+6)

SolvedDofsPrime = list(range(ndof+6+6))
SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,Fixed_Dofs)
SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,dofS2)
SolvedDofsPrime = scipy.setdiff1d(SolvedDofsPrime,dofS3)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
#print silex_lib_elt.stiffnessmatrix.__doc__
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

Kprime=scipy.sparse.csc_matrix(R.T*K*R)

#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

Qprime[SolvedDofsPrime] = mumps.spsolve(Kprime[SolvedDofsPrime,:][:,SolvedDofsPrime],Fprime[SolvedDofsPrime])

toc = time.clock()
print("time to solve the problem:",toc-tic)

Q=R*Qprime


#############################################################################
#       compute smooth stress and error in elements
#############################################################################
tic = time.clock()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)

toc = time.clock()
print("time to compute stres and error:",toc-tic)
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
##load=scipy.zeros((nnodes,ndim))
##load[range(nnodes),0]=F[list(range(0,ndof,3))]
##load[range(nnodes),1]=F[list(range(1,ndof,3))]
##load[range(nnodes),2]=F[list(range(2,ndof,3))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[range(nelem),[6]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[6]],'nodal',1,'Smooth Sigma V.M.'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      #[load,'nodal',ndim,'Force'],
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



