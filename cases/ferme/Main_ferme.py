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

#import silex_lib_truss as silex_lib_elt
import silex_lib_truss_python as silex_lib_elt

#############################################################################
print("SILEX CODE - calcul d'une ferme de charpente")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='ferme'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_ferme'

# choose the element type
eltype=1

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

# OR define by hands:
#nodes=scipy.array([[0.0  ,	   0.0	],
#                    [7.0/4.0  , 0.0	],
#                    [7.0/2.0   , 0.0	],
#                    [3.0*7.0/4  , 0.0],
#                    [7.0   , 0.0	],
#                    [7.0/4.0   , 0.75],
#                    [7.0/2.0    , 1.5],
#                    [3.0*7.0/4   , 0.75]])

#elements=scipy.array([[     1 ,	 2],
#                          [     2 ,	 3],
#                          [     3 ,	 4],
#                          [     4 ,	 5],
#                          [     1 ,	 6],
#                          [     6 ,	 7],
#                          [     7 ,	 8],
#                          [     8 ,	 5],
#                          [     6 ,	 2],
#                          [     8 ,	 4],
#                          [     6 ,	 3],
#                          [     7 ,	 3],
#                          [     3 ,	 8]])

#silex_lib_gmsh.WriteResults('toto',nodes,elements,1)

# Define material
Young  = 210000e6

# define geometry of the cross section
#      
#         Cross section of the rods
#      
#                  a
#        ___________________________ ___________________________
#       |                          ||                          |
#       |  h                       ||                          |
#       |______________________    ||    ______________________|
#                             |    ||    |
#                             |    ||    |
#                             |    ||    |
#                             |    || b  |
#                             |    ||    |
#                             |    ||    |
#                             |    ||    |
#                             |    ||    |
#                             |    ||    |
#                             |  h ||    |
#                             |____||____|

a = 40e-3  # m
b = a      # m
h = 4e-3   # m

# compute cross section
S1      = a*h
S2      = (b-h)*h
Section = 2*(S1+S2) # area of the section
# Inertia
yg1     = b-h/2
yg2     = (b-h)/2
yg      = ( yg1*S1 + yg2*S2 )/(S1+S2)
I1      = h*(a**3)/12  + S1*((yg1-yg)**2)
I2      = (b-h)*(h**3)/12 + S2*((yg2-yg)**2)
Inertia = 2*(I1+I2)

# Boundary conditions
IdNodesFixed_x=scipy.array([1,5],dtype=int)
IdNodesFixed_y=scipy.array([1,5],dtype=int)

# glass
glassweight=15.0 # Kg/m^2
glassarea=(3.8*2)*3.0 # m^2

# Define load case:
snow_press   = 450     #  N/m^2 = Pa - Paris

# weight of glass
weight = -glassweight*9.81*glassarea
# snow on a wire
snow = -snow_press*glassarea

LoadOnOneNode=2.0*((weight+snow)/4.0)/2.0

# If the user wants to have only the mesh for gmsh, uncomment next line
silex_lib_gmsh.WriteResults('maillage_seul',nodes,elements,eltype)

# Load on x direction
LoadX=[]

# Load on y direction
LoadY=[[6,LoadOnOneNode],
       [7,LoadOnOneNode],
       [8,LoadOnOneNode],
       ]

#############################################################################
#      expert part: define load and prescribed displacements
#############################################################################
#      initialisations
#############################################################################
# get number of nodes, dof and elements from the mesh
nnodes = nodes.shape[0]
ndof   = nnodes*2
nelem  = elements.shape[0]
print("Number of nodes:",nnodes)
print("Number of elements:",nelem)

# define fixed dof
Fixed_Dofs = scipy.hstack([(IdNodesFixed_x-1)*2,(IdNodesFixed_y-1)*2+1])

# define free dof
SolvedDofs = scipy.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=scipy.zeros(ndof)

# initialize force vector
F=scipy.zeros((ndof))

for i in range(len(LoadX)):
    F[(LoadX[i][0]-1)*2]=LoadX[i][1]

for i in range(len(LoadY)):
    F[(LoadY[i][0]-1)*2+1]=LoadY[i][1]

#############################################################################
#      compute stiffness matrix
#############################################################################

Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,Section])
K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

#############################################################################
#       Solve the problem
#############################################################################

Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

#############################################################################
#       compute normal force, stress, buckling limit load
#############################################################################
NormalForce,Sigma,Fcr=silex_lib_elt.compute_normal_force_stress_buckling(nodes,elements,[Young,Section,Inertia],Q)


#############################################################################
# write results in a gmsh file
#############################################################################

# displacement written on 2 columns:
disp=scipy.zeros((nnodes,2))
disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
disp[range(nnodes),1]=Q[list(range(1,ndof,2))]

# external forces written on 2 columns:
load=scipy.zeros((nnodes,2))
load[range(nnodes),0]=F[list(range(0,ndof,2))]
load[range(nnodes),1]=F[list(range(1,ndof,2))]

# get normal forces for compressive elements only, becomes positive
CompressiveElts=abs(NormalForce)*(scipy.sign(NormalForce)-1)/(-2.0)

# compute ratio between normal force and buckling limit load for compressive elements only
Ratio=CompressiveElts/Fcr

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
                  [load,'nodal',2,'forces'],
                  [Sigma,'elemental',1,'Normal stress'],
                  [NormalForce,'elemental',1,'Normal Force'],
                  [CompressiveElts,'elemental',1,'Compressive elements'],
                  [Fcr,'elemental',1,'Buckling limit load'],
                  [Ratio,'elemental',1,'Ratio'],
                  ]

# write the mesh and the results in a gmsh-format file
#truss_lib.WriteResults2Gmsh(ResultsFileName,nodes,elements,fields_to_write)
# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)

