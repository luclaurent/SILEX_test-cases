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
import silex_lib_extra_python

#############################################################################
print ("SILEX CODE - calcul d un piston avec des tet4")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='nervure-tet4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_nervure_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=3

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

nodes_tmp=nodes.copy()
nodes[:,1]=nodes[:,2].copy()
nodes[:,2]=nodes_tmp[:,1].copy()

t0=time.clock()
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)
t1=time.clock()

print ("time for reading elements:",t1-t0)


# read surfaces where to impose boundary conditions
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2) # charge
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3) # trou bord attaque
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,4) # trou barod fuite
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,5) # toutes les surfaces pour la sortie
elementsS6,IdnodeS6=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,6) # surface de "devant"

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,2)

# Define material
Young  = 200000.0
nu     = 0.3

# Boundary conditions
IdNodesFixed_x=scipy.hstack([IdnodeS3,IdnodeS4])
IdNodesFixed_y=[]
IdNodesFixed_z=scipy.hstack([IdnodeS6])
R3=silex_lib_extra_python.turn_dof3D(IdnodeS3,nodes,[(159.78+156.36)/2.0,3.02])
R4=silex_lib_extra_python.turn_dof3D(IdnodeS4,nodes,[(336.90+333.48)/2.0,3.02])

# compute external forces from pressure
press = -1 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0,1.0,0.0]
F = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)


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
#Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1,(IdNodesFixed_z-1)*3+2])
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_z-1)*3+2])

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

K=scipy.sparse.csc_matrix(R4.T*R3.T*K*R3*R4)
#############################################################################
#       Solve the problem
#############################################################################

tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs], use_umfpack=True)
Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
toc = time.clock()
print ("time to solve the problem:",toc-tic)
Q=R3*R4*Q
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
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],[load,'nodal',ndim,'Force'],
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
if flag_write_fields==3:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaNodes[range(nnodes),[6]],'nodal',1,'Smooth_Sigma_VM_on_+H/2'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)
#elementsshift = scipy.vstack([[elementsS6[:,1]],[elementsS6[:,0]],[elementsS6[:,2]]]).T
#elementsx3d = scipy.vstack([elementsS6, elementsshift])
silex_lib_gmsh.WriteResults(ResultsFileName+'Surf_Model',nodes,elementsS5,2,fields_to_write)

toc = time.clock()
print ("time to write results:",toc-tic)
print ("----- END -----")



