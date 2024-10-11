#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import sys
sys.path.append('../../librairies')

import silex_lib_tri3_python as silex_lib_elt
#import silex_lib_tri3 as silex_lib_elt
import silex_lib_gmsh

#############################################################################
print("SILEX CODE - calcul d'une plaque trouee avec des tri3")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='motif-tri3'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_motif-super-tri3'

# choose the element type
eltype=2

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# read surfaces where to impose boundary conditions
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,3)


# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_Surface_mesh',nodes,elements,2)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S1_mesh',nodes,elementsS1,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S2_mesh',nodes,elementsS2,1)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_S3_mesh',nodes,elementsS3,1)

# Define material
Young  = 70000.0
nu     = 0.27
thickness = 5.0

# Boundary conditions

# define fixed dof
IdNodesFixed_x=IdnodeS3
IdNodesFixed_y=IdnodeS3

F=silex_lib_elt.forceonline(nodes,elementsS2,[0.0,4500.0/100.0,0.0,4500.0/100.0],[100.0,0.0,100.0,100.0])

toc = time.clock()
print("time for the reading data part:",toc-tic)

tic0 = time.clock()
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
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic = time.clock()

Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness])

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

toc = time.clock()
print("time to compute the stiffness matrix:",toc-tic)

InterfaceIdNodes=scipy.hstack([IdnodeS3,IdnodeS2])

InterfaceDofs=scipy.hstack([(InterfaceIdNodes-1)*2,(InterfaceIdNodes-1)*2+1])
InternalDofs=scipy.setdiff1d(range(ndof),InterfaceDofs)

#############################################################################
#       compute "super" matrix condensed on interface
#############################################################################
InvK_internal_internal=scipy.sparse.linalg.inv(K[InternalDofs,:][:,InternalDofs])
Ksuper = K[InterfaceDofs,:][:,InterfaceDofs] - K[InterfaceDofs,:][:,InternalDofs]*InvK_internal_internal*K[InternalDofs,:][:,InterfaceDofs]

#############################################################################
#       compute the problem for 1 square+hole
#############################################################################

ndofsuper = len(InterfaceDofs)

Fixed_super_Dofs  = scipy.intersect1d(Fixed_Dofs,InterfaceDofs)
Solved_super_Dofs = scipy.intersect1d(SolvedDofs,InterfaceDofs)

#Qsuper = scipy.zeros(ndofsuper)
#Fsuper = F[InterfaceDofs]
InterfaceDofs=list(InterfaceDofs)
Solved_super_Dofs=list(Solved_super_Dofs)

Solved_super_local_Dofs=[]
for i in range(len(Solved_super_Dofs)):
    Solved_super_local_Dofs.append(InterfaceDofs.index(Solved_super_Dofs[i]))

Q[Solved_super_Dofs] = mumps.spsolve(Ksuper[Solved_super_local_Dofs,:][:,Solved_super_local_Dofs],F[Solved_super_Dofs])

#Q[InterfaceDofs] = Qsuper
InvK_internal_interface=InvK_internal_internal*K[InternalDofs,:][:,InterfaceDofs]
Q[InternalDofs] = -InvK_internal_interface*Q[InterfaceDofs]



#############################################################################
#       Solve the problem
#############################################################################

#tic = time.clock()
#Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
#toc = time.clock()
#print("time to solve the problem:",toc-tic)

#############################################################################
#       compute stress, smooth stress, strain and error
#############################################################################
tic = time.clock()

SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness],Q)

toc = time.clock()
print("time to compute stresses:",toc-tic)
print("The global error is:",ErrorGlobal)


#############################################################################
#         Write results to gmsh format
#############################################################################
tic = time.clock()

# displacement written on 2 columns:
disp=scipy.zeros((nnodes,2))
disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
disp[range(nnodes),1]=Q[list(range(1,ndof,2))]

load=scipy.zeros((nnodes,ndim))
load[range(nnodes),0]=F[list(range(0,ndof,2))]
load[range(nnodes),1]=F[list(range(1,ndof,2))]

if flag_write_fields==0:
    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                      [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                      [ErrorElem,'elemental',1,'error'],
                      ]

if flag_write_fields==1:
    fields_to_write=[ [disp,'nodal',ndim,'Displacement'],
                  [load,'nodal',2,'Force'],
                  [SigmaElem[range(nelem),[0]],'elemental',1,'Sigma xx'],
                  [SigmaElem[range(nelem),[1]],'elemental',1,'Sigma yy'],
                  [SigmaElem[range(nelem),[2]],'elemental',1,'Sigma xy'],
                  [SigmaElem[range(nelem),[3]],'elemental',1,'Sigma V.M.'],
                  [SigmaNodes[range(nnodes),[0]],'nodal',1,'Sigma xx Smooth'],
                  [SigmaNodes[range(nnodes),[1]],'nodal',1,'Sigma yy Smooth'],
                  [SigmaNodes[range(nnodes),[2]],'nodal',1,'Sigma xy Smooth'],
                  [SigmaNodes[range(nnodes),[3]],'nodal',1,'Sigma V.M. Smooth'],
                  [EpsilonElem[range(nelem),[0]],'elemental',1,'Epsilon xx'],
                  [EpsilonElem[range(nelem),[1]],'elemental',1,'Epsilon yy'],
                  [EpsilonElem[range(nelem),[2]]/2.0,'elemental',1,'Epsilon xy'],
                  [EpsilonNodes[range(nnodes),[0]],'nodal',1,'Epsilon xx Smooth'],
                  [EpsilonNodes[range(nnodes),[1]],'nodal',1,'Epsilon yy Smooth'],
                  [EpsilonNodes[range(nnodes),[2]]/2.0,'nodal',1,'Epsilon xy Smooth'],
                  [ErrorElem,'elemental',1,'Error']
                  ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName+'_1_hole',nodes,elements,eltype,fields_to_write)

#############################################################################
#       compute the problem for 3 square+hole
#############################################################################

# Make the whole new mesh:
nodes_3holes = scipy.vstack([nodes,nodes+[100.0,0.0],nodes+[200.0,0.0]])
elements_3holes = scipy.vstack([elements,elements+nnodes,elements+2*nnodes])

if len(IdnodeS2)!=len(IdnodeS3):
    print("Be carefull...................................................")

elt1 = list(range(0,len(InterfaceIdNodes),1))
elt2 = list(range(int(0.5*len(InterfaceIdNodes)),int(1.5*len(InterfaceIdNodes)),1))
elt3 = list(range(int(1*len(InterfaceIdNodes)),int(2*len(InterfaceIdNodes)),1))

InterfaceIdNodes_1_left  = scipy.array(list(range(0,int(0.5*len(InterfaceIdNodes)),1)))
#InterfaceIdNodes_1_right = IdnodeS2
#InterfaceIdNodes_2_left  = IdnodeS3+nnodes
#InterfaceIdNodes_2_right = IdnodeS2+nnodes
#InterfaceIdNodes_3_left  = IdnodeS3+2*nnodes
InterfaceIdNodes_3_right = scipy.array(list(range(int(1.5*len(InterfaceIdNodes)),int(2*len(InterfaceIdNodes)),1)))

InterfaceIdDoflocal_1  = scipy.hstack([scipy.array(elt1)*2,scipy.array(elt1)*2+1])
InterfaceIdDoflocal_2  = scipy.hstack([scipy.array(elt2)*2,scipy.array(elt2)*2+1])
InterfaceIdDoflocal_3  = scipy.hstack([scipy.array(elt3)*2,scipy.array(elt3)*2+1])


elements_super_3holes=scipy.zeros((3,len(elt1)))
elements_super_3holes[0,:]=elt1
elements_super_3holes[1,:]=elt2
elements_super_3holes[2,:]=elt3

nbelem_super_3holes=3

nnodes_super_3holes=len( scipy.unique(scipy.array(elements_super_3holes)) )
ndof_super_3holes=nnodes_super_3holes*2

# define fixed dof
Fixed_super_3holes_Dofs = scipy.hstack([(InterfaceIdNodes_1_left)*2,(InterfaceIdNodes_1_left)*2+1])
# define free dof
Solved_super_3holes_Dofs = scipy.setdiff1d(range(ndof_super_3holes),Fixed_super_3holes_Dofs)

Ik = scipy.zeros(4*elements_super_3holes.shape[1]**2*nbelem_super_3holes,dtype=int)
Jk = scipy.zeros(4*elements_super_3holes.shape[1]**2*nbelem_super_3holes,dtype=int)
Vk = scipy.zeros(4*elements_super_3holes.shape[1]**2*nbelem_super_3holes,dtype=float)

p=0

for e in range(nbelem_super_3holes):
    idnodes=elements_super_3holes[e]
    dofx = (idnodes)*2
    dofy = (idnodes)*2+1
    dofelem = scipy.hstack([dofx,dofy])

    for i in range(len(dofelem)):
        for j in range(len(dofelem)):
            Ik[p]=dofelem[i]
            Jk[p]=dofelem[j]
            Vk[p]=Ksuper[i,j]
            p=p+1

K_super_3holes=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_super_3holes,ndof_super_3holes) )

Q_super_3holes=scipy.zeros(ndof_super_3holes)
F_3holes=scipy.zeros(ndof_super_3holes)
F_3holes[scipy.hstack([InterfaceIdNodes_3_right*2,InterfaceIdNodes_3_right*2+1])]=F[scipy.hstack([(IdnodeS2-1)*2,(IdnodeS2-1)*2+1])]

Q_super_3holes[Solved_super_3holes_Dofs] = mumps.spsolve(K_super_3holes[Solved_super_3holes_Dofs,:][:,Solved_super_3holes_Dofs],F_3holes[Solved_super_3holes_Dofs])
#Q_super_3holes[Solved_super_3holes_Dofs] = scipy.sparse.linalg.spsolve(K_super_3holes[Solved_super_3holes_Dofs,:][:,Solved_super_3holes_Dofs],F_3holes[Solved_super_3holes_Dofs])

# compute into super elements

Qall=scipy.zeros(nodes_3holes.shape[0]*2)

InternalDofs_1=InternalDofs
InternalDofs_2=InternalDofs+ndof
InternalDofs_3=InternalDofs+2*ndof

InterfaceDofs_1=scipy.array(InterfaceDofs)
InterfaceDofs_2=scipy.array(InterfaceDofs)+ndof
InterfaceDofs_3=scipy.array(InterfaceDofs)+2*ndof

# Elt 1:
Qall[InternalDofs_1] = -InvK_internal_interface*Q_super_3holes[InterfaceIdDoflocal_1]
Qall[InternalDofs_2] = -InvK_internal_interface*Q_super_3holes[InterfaceIdDoflocal_2]
Qall[InternalDofs_3] = -InvK_internal_interface*Q_super_3holes[InterfaceIdDoflocal_3]

Qall[InterfaceDofs_1] = Q_super_3holes[InterfaceIdDoflocal_1]
Qall[InterfaceDofs_2] = Q_super_3holes[InterfaceIdDoflocal_2]
Qall[InterfaceDofs_3] = Q_super_3holes[InterfaceIdDoflocal_3]

# displacement written on 2 columns:
dispAll=scipy.zeros((nodes_3holes.shape[0],2))
dispAll[range(nodes_3holes.shape[0]),0]=Qall[list(range(0,3*ndof,2))]
dispAll[range(nodes_3holes.shape[0]),1]=Qall[list(range(1,3*ndof,2))]

loadAll=scipy.zeros((nodes_3holes.shape[0],2))
#loadAll[range(nodes_3holes.shape[0]),0]=F[list(range(0,ndof,2))]
loadAll[IdnodeS2+2*nnodes-1,1]=F[(IdnodeS2-1)*2+1]

fields_to_write=[ [dispAll,'nodal',ndim,'displacement'], [loadAll,'nodal',ndim,'load']]

silex_lib_gmsh.WriteResults('Results_super_elements_3_holes',nodes_3holes,elements_3holes,eltype,fields_to_write)

toc = time.clock()
print("time to write results:",toc-tic)
print("----- END -----")
