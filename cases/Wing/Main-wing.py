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

import silex_lib_dkt as silex_lib_elt
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une aile d'avion")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='wing'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_wing_dkt'

# choose the element type
eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# nervures
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# peau inferieure longeron
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)

# ame longeron
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)

# premiere nervure
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# peau superieure longeron
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,5)

# bord d'attaque
elementsS6,IdnodeS6=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,6)

# bord de fuite
elementsS7,IdnodeS7=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,7)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1_nervures',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2_peau_inferieure_longeron',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3_ame_longeron',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf4_nervure_emplanture',nodes,elementsS4,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf5_peau_superieure_longeron',nodes,elementsS5,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf6_bord_attaque',nodes,elementsS6,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf7_bord_fuite',nodes,elementsS7,2)

elements=scipy.vstack([elementsS1,elementsS2,elementsS3,elementsS4,elementsS5,elementsS6,elementsS7])

# Define material
Young     = 200000.0
nu        = 0.3
thickness = 2.0

Young1     = Young
nu1        = nu
thickness1 = thickness

Young2     = Young
nu2        = nu
thickness2 = thickness

Young3     = Young
nu3        = nu
thickness3 = thickness

Young4     = Young
nu4        = nu
thickness4 = thickness

Young5     = Young
nu5        = nu
thickness5 = thickness

Young6     = Young
nu6        = nu
thickness6 = thickness

Young7     = Young
nu7        = nu
thickness7 = thickness

# Boundary conditions
IdNodesFixed_x=IdnodeS4
IdNodesFixed_y=IdnodeS4
IdNodesFixed_z=IdnodeS4
IdNodesFixed_rotx=IdnodeS4
IdNodesFixed_roty=IdnodeS4
IdNodesFixed_rotz=IdnodeS4

# compute external forces from pressure
press = 0.01 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [0.0,1.0,0.0]
F = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)


toc = time.clock()
print("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*6
nelem  = elementsS1.shape[0]+elementsS2.shape[0]+elementsS3.shape[0]+elementsS4.shape[0]+elementsS5.shape[0]+elementsS6.shape[0]+elementsS7.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

#      CLEAN MESH
Id_nodes_used=scipy.unique(elements)
Id_nodes_nonused=scipy.setdiff1d(range(1,nnodes),Id_nodes_used)
IdNodesFixed_x=scipy.hstack([IdNodesFixed_x,Id_nodes_nonused])
IdNodesFixed_y=scipy.hstack([IdNodesFixed_y,Id_nodes_nonused])
IdNodesFixed_z=scipy.hstack([IdNodesFixed_z,Id_nodes_nonused])
IdNodesFixed_rotx=scipy.hstack([IdNodesFixed_rotx,Id_nodes_nonused])
IdNodesFixed_roty=scipy.hstack([IdNodesFixed_roty,Id_nodes_nonused])
IdNodesFixed_rotz=scipy.hstack([IdNodesFixed_rotz,Id_nodes_nonused])

# define fixed dof
Fixed_Dofs = scipy.hstack([
    (scipy.array(IdNodesFixed_x)-1)*6,
    (scipy.array(IdNodesFixed_y)-1)*6+1,
    (scipy.array(IdNodesFixed_z)-1)*6+2,
    (scipy.array(IdNodesFixed_rotx)-1)*6+3,
    (scipy.array(IdNodesFixed_roty)-1)*6+4,
    (scipy.array(IdNodesFixed_rotz)-1)*6+5])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
#print (silex_lib_elt.stiffnessmatrix.__doc__)
Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementsS1,[Young1,nu1,thickness1,1000.0])
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementsS2,[Young2,nu2,thickness2,1000.0])
Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementsS3,[Young3,nu3,thickness3,1000.0])
Ik4,Jk4,Vk4,Vm4=silex_lib_elt.stiffnessmatrix(nodes,elementsS4,[Young4,nu4,thickness4,1000.0])
Ik5,Jk5,Vk5,Vm5=silex_lib_elt.stiffnessmatrix(nodes,elementsS5,[Young5,nu5,thickness5,1000.0])
Ik6,Jk6,Vk6,Vm6=silex_lib_elt.stiffnessmatrix(nodes,elementsS6,[Young6,nu6,thickness6,1000.0])
Ik7,Jk7,Vk7,Vm7=silex_lib_elt.stiffnessmatrix(nodes,elementsS7,[Young7,nu7,thickness7,1000.0])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K1=scipy.sparse.csc_matrix( (Vk1,(Ik1,Jk1)), shape=(ndof,ndof) ,dtype=float)
K2=scipy.sparse.csc_matrix( (Vk2,(Ik2,Jk2)), shape=(ndof,ndof) ,dtype=float)
K3=scipy.sparse.csc_matrix( (Vk3,(Ik3,Jk3)), shape=(ndof,ndof) ,dtype=float)
K4=scipy.sparse.csc_matrix( (Vk4,(Ik4,Jk4)), shape=(ndof,ndof) ,dtype=float)
K5=scipy.sparse.csc_matrix( (Vk5,(Ik5,Jk5)), shape=(ndof,ndof) ,dtype=float)
K6=scipy.sparse.csc_matrix( (Vk6,(Ik6,Jk6)), shape=(ndof,ndof) ,dtype=float)
K7=scipy.sparse.csc_matrix( (Vk7,(Ik7,Jk7)), shape=(ndof,ndof) ,dtype=float)

K=K1+K2+K3+K4+K5+K6+K7

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

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness,7000.0],Q)

toc = time.clock()
print("time to compute stres and error:",toc-tic)
#print "The global error is:",ErrorGlobal
print("Total time for the computational part:",toc-tic0)

#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 3 columns:
disp=scipy.zeros((nnodes,ndim))
disp[range(nnodes),0]=Q[list(range(0,ndof,6))]
disp[range(nnodes),1]=Q[list(range(1,ndof,6))]
disp[range(nnodes),2]=Q[list(range(2,ndof,6))]

# external load written on 3 columns:
load=scipy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,6))]
load[range(nnodes),1]=F[list(range(1,ndof,6))]
load[range(nnodes),2]=F[list(range(2,ndof,6))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaNodes[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                      [SigmaNodes[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2']
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [load,'nodal',ndim,'Force'],
                      [SigmaElem[:,0],'elemental',1,'Sigma xx on +H/2'],
                      [SigmaElem[:,1],'elemental',1,'Sigma yy on +H/2'],
                      [SigmaElem[:,2],'elemental',1,'Sigma xy on +H/2'],
                      [SigmaElem[:,3],'elemental',1,'Sigma 1 on +H/2'],
                      [SigmaElem[:,4],'elemental',1,'Sigma 2 on +H/2'],
                      [SigmaElem[:,5],'elemental',1,'Alpha0 on +H/2'],
                      [SigmaElem[:,6],'elemental',1,'Sigma V.M. on +H/2'],
                      [SigmaElem[:,7],'elemental',1,'Sigma xx on -H/2'],
                      [SigmaElem[:,8],'elemental',1,'Sigma yy on -H/2'],
                      [SigmaElem[:,9],'elemental',1,'Sigma xy on -H/2'],
                      [SigmaElem[:,10],'elemental',1,'Sigma 1 on -H/2'],
                      [SigmaElem[:,11],'elemental',1,'Sigma 2 on -H/2'],
                      [SigmaElem[:,12],'elemental',1,'Alpha0 on -H/2'],
                      [SigmaElem[:,13],'elemental',1,'Sigma V.M. on -H/2'],
                      [SigmaNodes[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                      [SigmaNodes[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                      [EpsilonElem[:,0],'elemental',1,'Epsilon xx on +H/2'],
                      [EpsilonElem[:,1],'elemental',1,'Epsilon yy on +H/2'],
                      [EpsilonElem[:,2],'elemental',1,'Epsilon xy on +H/2'],
                      [EpsilonElem[:,3],'elemental',1,'Epsilon xx on -H/2'],
                      [EpsilonElem[:,4],'elemental',1,'Epsilon yy on -H/2'],
                      [EpsilonElem[:,5],'elemental',1,'Epsilon xy on -H/2']
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")



