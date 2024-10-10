#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

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
MeshFileName='aile-maquette-balsa'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_aile-maquette-balsa'

# choose the element type
eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# longeron bord de fuite
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# longeron central
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,2)

# longeron bord d'attaque
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)

# nervures sauf emplanture
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# nervure emplanture
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,5)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf1',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf3',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf4',nodes,elementsS4,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_surf5',nodes,elementsS5,2)

elements=scipy.vstack([elementsS1,elementsS2,elementsS3,elementsS4,elementsS5])

silex_lib_gmsh.WriteResults(ResultsFileName+'_complet',nodes,elements,2)

# Define material
Young     = 4300e6 #Pa
nu        = 0.3

thickness1 = 2.0e-3
thickness2 = 4.0e-3
thickness3 = 2.0e-3
thickness4 = 2.0e-3
thickness5 = 2.0e-3

# Boundary conditions
IdNodesFixed_x=IdnodeS5
IdNodesFixed_y=IdnodeS5
IdNodesFixed_z=IdnodeS5
IdNodesFixed_rotx=IdnodeS5
IdNodesFixed_roty=IdnodeS5
IdNodesFixed_rotz=IdnodeS5

# compute external forces from pressure
#press = 0.01 #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
#direction = [0.0,1.0,0.0]
#F = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)

F=scipy.zeros(6*nodes.shape[0])
#F[(581-1)*6+2]=-20e-3*9.81 # force en A
F[(602-1)*6+2]=-20e-3*9.81 # force en B
#F[(566-1)*6+2]=-20e-3*9.81 # force en C


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
nelem  = elementsS1.shape[0]+elementsS2.shape[0]+elementsS3.shape[0]+elementsS4.shape[0]+elementsS5.shape[0]
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
Ik1,Jk1,Vk1,Vm1=silex_lib_elt.stiffnessmatrix(nodes,elementsS1,[Young,nu,thickness1,140.0])
Ik2,Jk2,Vk2,Vm2=silex_lib_elt.stiffnessmatrix(nodes,elementsS2,[Young,nu,thickness2,140.0])
Ik3,Jk3,Vk3,Vm3=silex_lib_elt.stiffnessmatrix(nodes,elementsS3,[Young,nu,thickness3,140.0])
Ik4,Jk4,Vk4,Vm4=silex_lib_elt.stiffnessmatrix(nodes,elementsS4,[Young,nu,thickness4,140.0])
Ik5,Jk5,Vk5,Vm5=silex_lib_elt.stiffnessmatrix(nodes,elementsS5,[Young,nu,thickness5,140.0])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K1=scipy.sparse.csc_matrix( (Vk1,(Ik1,Jk1)), shape=(ndof,ndof) ,dtype=float)
K2=scipy.sparse.csc_matrix( (Vk2,(Ik2,Jk2)), shape=(ndof,ndof) ,dtype=float)
K3=scipy.sparse.csc_matrix( (Vk3,(Ik3,Jk3)), shape=(ndof,ndof) ,dtype=float)
K4=scipy.sparse.csc_matrix( (Vk4,(Ik4,Jk4)), shape=(ndof,ndof) ,dtype=float)
K5=scipy.sparse.csc_matrix( (Vk5,(Ik5,Jk5)), shape=(ndof,ndof) ,dtype=float)

K=K1+K2+K3+K4+K5

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

output1=silex_lib_elt.compute_stress_strain_error(nodes,scipy.vstack([elementsS1,elementsS3,elementsS4,elementsS5]),[Young,nu,thickness1,140.0],Q)
output2=silex_lib_elt.compute_stress_strain_error(nodes,elementsS2,[Young,nu,thickness2,140.0],Q)

SigmaElemS1    = output1[0]
SigmaNodesS1   = output1[1]
EpsilonElemS1  = output1[2]
EpsilonNodesS1 = output1[3]
ErrorElemS1    = output1[4]
ErrorGlobalS1  = output1[5]
DirPrinPlu1S1  = output1[6]
DirPrinPlu2S1  = output1[7]
DirPrinMin1S1  = output1[8]
DirPrinMin2S1  = output1[9]

SigmaElemS2    = output2[0]
SigmaNodesS2   = output2[1]
EpsilonElemS2  = output2[2]
EpsilonNodesS2 = output2[3]
ErrorElemS2    = output2[4]
ErrorGlobalS2  = output2[5]
DirPrinPlu1S2  = output2[6]
DirPrinPlu2S2  = output2[7]
DirPrinMin1S2  = output2[8]
DirPrinMin2S2  = output2[9]

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

SigmaNodesS1[IdnodeS2-1,0]=SigmaNodesS2[IdnodeS2-1,0]
SigmaNodesS1[IdnodeS2-1,1]=SigmaNodesS2[IdnodeS2-1,1]
DirPrinPlu1S1=scipy.vstack([DirPrinPlu1S1,DirPrinPlu1S2])
DirPrinPlu2S1=scipy.vstack([DirPrinPlu2S1,DirPrinPlu2S2])
DirPrinMin1S1=scipy.vstack([DirPrinMin1S1,DirPrinMin1S2])
DirPrinMin2S1=scipy.vstack([DirPrinMin2S1,DirPrinMin2S2])
SigmaElemS1=scipy.vstack([SigmaElemS1,SigmaElemS2])

if flag_write_fields==0:
    fields_to_writeS1=[[disp,'nodal',ndim,'displacement'],
                       [load,'nodal',ndim,'Force'],
                       [SigmaNodesS1[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                       [SigmaNodesS1[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                       [DirPrinPlu1S1,'elemental',ndim,'Dir. Princ.1 +H/2'],
                       [DirPrinPlu2S1,'elemental',ndim,'Dir. Princ.2 +H/2'],
                       [DirPrinMin1S1,'elemental',ndim,'Dir. Princ.1 -H/2'],
                       [DirPrinMin2S1,'elemental',ndim,'Dir. Princ.2 -H/2'],
                       [SigmaElemS1[:,3],'elemental',1,'Sigma 1 on +H/2'],
                       [SigmaElemS1[:,4],'elemental',1,'Sigma 2 on +H/2'],
                       [SigmaElemS1[:,10],'elemental',1,'Sigma 1 on -H/2'],
                       [SigmaElemS1[:,11],'elemental',1,'Sigma 2 on -H/2'],
                      ]

    fields_to_writeS2=[[disp,'nodal',ndim,'displacement'],
                       [load,'nodal',ndim,'Force'],
                       [SigmaNodesS2[:,0],'nodal',1,'Smooth Sigma V.M. on +H/2'],
                       [SigmaNodesS2[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                       [DirPrinPlu1S2,'elemental',ndim,'Dir. Princ.1 +H/2'],
                       [DirPrinPlu2S2,'elemental',ndim,'Dir. Princ.2 +H/2'],
                       [DirPrinMin1S2,'elemental',ndim,'Dir. Princ.1 -H/2'],
                       [DirPrinMin2S2,'elemental',ndim,'Dir. Princ.2 -H/2'],
                       [SigmaElemS2[:,3],'elemental',1,'Sigma 1 on +H/2'],
                       [SigmaElemS2[:,4],'elemental',1,'Sigma 2 on +H/2'],
                       [SigmaElemS2[:,10],'elemental',1,'Sigma 1 on -H/2'],
                       [SigmaElemS2[:,11],'elemental',1,'Sigma 2 on -H/2'],
                      ]
    

silex_lib_gmsh.WriteResults(ResultsFileName+'_epaisseur_2',nodes,scipy.vstack([elementsS1,elementsS3,elementsS4,elementsS5,elementsS2]),2,fields_to_writeS1)
silex_lib_gmsh.WriteResults(ResultsFileName+'_epaisseur_4',nodes,elementsS2,2,fields_to_writeS2)

print('Point A : deplacement MAXI = ',disp[581-1,:]*4300/3000)
print('Point B : deplacement MAXI = ',disp[602-1,:]*4300/3000)
print('Point C : deplacement MAXI = ',disp[566-1,:]*4300/3000)
print('Point A : deplacement MINI = ',disp[581-1,:])
print('Point B : deplacement MINI = ',disp[602-1,:])
print('Point C : deplacement MINI = ',disp[566-1,:])

toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")



