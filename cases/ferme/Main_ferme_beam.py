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

import silex_lib_beam_python as silex_lib_elt

#############################################################################
print("SILEX CODE - calcul d'une ferme de charpente - poutre Bernoulli")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='ferme'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_ferme_beam'

# choose the element type
eltype=1

# choose geometry dimension
ndim=2

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

#silex_lib_gmsh.WriteResults('toto',nodes,elements,1)

# Define material
Young  = 210000e6 #Pa

# cross section
Section = (80.0*40.0-2.0*36.0**2)*1e-6 # area of the section
Inertia = (80.0*40.0**3.0/12.0-2.0*(36.0**4.0/12.0+2.0**2.0*36.0**2.0))*1e-12

# Boundary conditions
IdNodesFixed_x=scipy.array([1,5],dtype=int)
IdNodesFixed_y=scipy.array([1,5],dtype=int)
# Load on x direction : [ [ node , force_x ]  ,  [ node , force_x ] , .....  ]
LoadX=[[1,0.0],
       [2,0.0],
       [3,0.0],
       [4,0.0],
       [5,0.0],
       [6,0.0],
       [7,0.0],
       [8,0.0]]

# Load on y direction : [ [ node , force_y ]  ,  [ node , force_y ] , .....  ]
LoadY=[[1,-3403.75/2.0],
       [2,0.0],
       [3,0.0],
       [4,0.0],
       [5,-3403.75/2.0],
       [6,-3403.75],
       [7,-3403.75],
       [8,-3403.75]]
       
# Load on theta direction : [ [ node , force_x ]  ,  [ node , force_x ] , .....  ]
#LoadTheta=[[1,0.0],
#       [2,0.0],
#       [3,0.0],
#       [4,0.0],
#       [5,0.0],
#       [6,0.0],
#       [7,0.0],
#       [8,0.0]]
#############################################################################
#      expert part: define load and prescribed displacements
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*3
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*3,(IdNodesFixed_y-1)*3+1])
# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

# initialize force vector
F=scipy.zeros((ndof))

for i in range(len(LoadX)):
    F[(LoadX[i][0]-1)*3]=LoadX[i][1]

for i in range(len(LoadY)):
    F[(LoadY[i][0]-1)*3+1]=LoadY[i][1]
    
#for i in range(len(LoadTheta)):
#    F[(LoadTheta[i][0]-1)*3+2]=LoadTheta[i][1]
    
#############################################################################
#      compute stiffness matrix
#############################################################################
Ik,Jk,Vk=silex_lib_elt.StiffnessMatrix(nodes,elements,[Young,Section,Inertia])
K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

#############################################################################
#       Solve the problem
#############################################################################

Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

#############################################################################
#       compute normal force, stress, buckling limit load
#############################################################################
#NormalForce,Sigma,Fcr=silex_lib_elt.compute_normal_force_stress_buckling(nodes,elements,[Young,Section,Inertia],Q)


#############################################################################
# write results in a gmsh file
#############################################################################

# displacement written on 2 columns:
disp=scipy.zeros((nnodes,2))
disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
# external forces written on 2 columns:
load=scipy.zeros((nnodes,2))
load[range(nnodes),0]=F[list(range(0,ndof,3))]
load[range(nnodes),1]=F[list(range(1,ndof,3))]
# get normal forces for compressive elements only, becomes positive
#CompressiveElts=abs(NormalForce)*(scipy.sign(NormalForce)-1)/(-2.0)

# compute ratio between normal force and buckling limit load for compressive elements only
#Ratio=CompressiveElts/Fcr

# syntax to define fields to write :
#fields_to_write=[ [variable_name1,'nodal' or'elemental' ,number of values per node,'name 1'],
#                  [variable_name2,'nodal' or'elemental' ,number of values per node,'name 2'],
#                  ...
#                  ]
#
# to write just the mesh with no fields :
#      fields_to_write=[]
# or even easier:
#      gmsh_lib.WriteResults('name.msh',nodes,elements)
                

# file out the fields list to write in the results gmsh-format file
fields_to_write=[ [disp,'nodal',2,'displacement'],
                  [load,'nodal',2,'forces']
                  ]

# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)
toc=time.clock()
print("Time for total computation",toc-tic)
