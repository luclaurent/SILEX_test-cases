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
MeshFileName='plaque-dkt'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_plaque_dkt'

# choose the element type
eltype=2

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)

# encastrement (fixe en x)
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,1)

# force
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)

# surface
elements,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)

# appui sur y
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,4)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'_1',nodes,elementsS1,1)
silex_lib_gmsh.WriteResults(ResultsFileName+'_2',nodes,elementsS2,1)
silex_lib_gmsh.WriteResults(ResultsFileName+'_3',nodes,elements,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'_4',nodes,elementsS4,1)

# Define material
Young     = 200000.0
nu        = 0.3
thickness = 5.0

# Boundary conditions
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS4
IdNodesFixed_z=IdnodeS1
IdNodesFixed_rotx=IdnodeS1
IdNodesFixed_roty=IdnodeS1
IdNodesFixed_rotz=IdnodeS1

forcez=0.0
forcex=10000.0
F=silex_lib_elt.forceonline(nodes,elementsS2,[forcex/100.0,0.0,forcez/100.0,forcex/100.0,0.0,forcez/100.0],[1000.0,0.0,0.0,1000.0,100.0,0.0])

toc = time.clock()
print("time for the user part:",toc-tic)

#############################################################################
#      EXPERT PART
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*6
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

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
Ik,Jk,Vk,Vm=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness,7000.0])
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

output=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu,thickness,7000.0],Q)

SigmaElem    = output[0]
SigmaNodes   = output[1]
EpsilonElem  = output[2]
EpsilonNodes = output[3]
ErrorElem    = output[4]
ErrorGlobal  = output[5]
DirPrinPlu1  = output[6]
DirPrinPlu2  = output[7]
DirPrinMin1  = output[8]
DirPrinMin2  = output[9]

toc = time.clock()
print("time to compute stres and error:",toc-tic)
#print("The global error is:",ErrorGlobal)
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
                      [SigmaNodes[:,1],'nodal',1,'Smooth Sigma V.M. on -H/2'],
                      [DirPrinPlu1,'elemental',ndim,'Dir. Princ.1 +H/2'],
                      [DirPrinPlu2,'elemental',ndim,'Dir. Princ.2 +H/2'],
                      [DirPrinMin1,'elemental',ndim,'Dir. Princ.1 -H/2'],
                      [DirPrinMin2,'elemental',ndim,'Dir. Princ.2 -H/2'],
                      [SigmaElem[:,3],'elemental',1,'Sigma 1 on +H/2'],
                      [SigmaElem[:,4],'elemental',1,'Sigma 2 on +H/2'],
                      [SigmaElem[:,10],'elemental',1,'Sigma 1 on -H/2'],
                      [SigmaElem[:,11],'elemental',1,'Sigma 2 on -H/2'],
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
                      [EpsilonElem[:,5],'elemental',1,'Epsilon xy on -H/2'],
                      [DirPrinPlu1,'elemental',ndim,'Dir. Princ.1 +H/2'],
                      [DirPrinPlu2,'elemental',ndim,'Dir. Princ.2 +H/2'],
                      [DirPrinMin1,'elemental',ndim,'Dir. Princ.1 -H/2'],
                      [DirPrinMin2,'elemental',ndim,'Dir. Princ.2 -H/2'],
                      ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

toc = time.clock()
print("time to write results: ",toc-tic)
print("Analytic displacement along x: ", forcex*1000.0/Young/(100.0*thickness)  )
print("Analytic displacement along z: ", forcez*1000.0**3/3.0/Young/(100.0*thickness**3/12.0)  )
print("Analytic stress xx for traction : ", forcex/100.0/thickness  )
print("Analytic stress xx for flexion  : ", forcez*1000.0/(100.0*thickness**3/12.0)*thickness/2.0  )
print("Analytic stress xx for traction+flexion on +H/2 : ", forcex/100.0/thickness-forcez*1000.0/(100.0*thickness**3/12.0)*thickness/2.0  )
print("Analytic stress xx for traction+flexion on -H/2 : ", forcex/100.0/thickness+forcez*1000.0/(100.0*thickness**3/12.0)*thickness/2.0  )
print("----- END -----")



