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

import silex_lib_truss_proj7 as silex_lib_elt
#import silex_lib_optim_fortran as silex_lib_optim

#############################################################################
print("SILEX CODE - calcul d'une ferme de charpente")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='truss'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_truss'

# choose the element type
eltype=1

# choose geometry dimension
ndim=2

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,1)

#silex_lib_gmsh.WriteResults('toto',nodes,elements,1)

# Define material
Young  = 2e11

# define geometry of the cross section
#      
#         Cross section of the rods
Section = 0.05**2 # area of the section
# Inertia
Inertia = 0.05**4/12

# Boundary conditions
IdNodesFixed_x=scipy.array([1,7,13,19],dtype=int)
IdNodesFixed_y=scipy.array([1,7,13,19],dtype=int)

# If the user wants to have only the mesh for gmsh, uncomment next line
#silex_lib_gmsh.WriteResults('maillage_seul',nodes,elements,eltype)

# Load on x direction
LoadX=[]

# Load on y direction
LoadY=[[24,-100]
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

Lfrac = 0.21
xe =scipy.ones(nelem)*Lfrac
penal = 3.0
YoungMin = 2e9

dC = scipy.zeros(nelem)


##############################################################################
#       Length study
##############################################################################
Lelem,sumL = silex_lib_elt.getlength(nodes,elements)


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
    

##for i in range(nelem):
#idnodes=elements[25,:]
#X=nodes[[idnodes-1],0][0]
#Y=nodes[[idnodes-1],1][0]
#ke = silex_lib_elt.ElementalStiffness(X,Y,Youngelem,Section)
    
#############################################################################
#       boucle d'itération
#############################################################################
loop = 0
loop_list=[]
XE_to_plot_list=[]
change_to_list=[]
change=1.0
while change>0.01 and loop<300:
    loop += 1
    loop_list.append(loop)
    
    # update elementql young modulus    
    Youngelem  = YoungMin+xe**penal*(Young-YoungMin)
    #Youngelem = xe**penal*Young
#############################################################################
#      compute stiffness matrix
#############################################################################
    Ik,Jk,Vk=silex_lib_elt.stiffnessmatrixoptim(nodes,elements,[Youngelem,Section])
    K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )
#############################################################################
#       Solve the problem
#############################################################################
    Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])
#############################################################################
#       Optimization
#############################################################################
    #dC = 2*Q.T*F-(penal*xe**(penal-1))*(Q.T*K*Q)
#############################################################################
#       compute normal force, stress, buckling limit load
#############################################################################
    NormalForce,Sigma,Fcr=silex_lib_elt.compute_normal_force_stress_buckling(nodes,elements,[Young,Section,Inertia],Q)
    StrEner=silex_lib_elt.getelementalstrainenergy(nodes,elements,[Youngelem,Section],Q)
    dC=(-penal*xe**(penal-1)*(Young-YoungMin))*StrEner
    #dC=-penal*xe**(penal-1)*Young
    xeold=xe
    xe=silex_lib_optim.oldoc88(xe,dC,Lelem,Lfrac,1)    
    #xe=silex_lib_elt.OptimalityCriteria(xe,dC,Lelem,Lfrac,nelem)
    XE_to_plot_list.append(xe.copy())
    change=scipy.linalg.norm(xe-xeold,scipy.inf)
    change_to_list.append(change)
    print('Loop :',loop, 'change :',change)

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
                    [xe,'elemental',1,'xe']
                  ]

# write the mesh and the results in a gmsh-format file
#truss_lib.WriteResults2Gmsh(ResultsFileName,nodes,elements,fields_to_write)
# write the mesh and the results in a gmsh-format file
silex_lib_gmsh.WriteResults(ResultsFileName,nodes,elements,eltype,fields_to_write)
silex_lib_gmsh.WriteResults2(ResultsFileName+'_Density_'+str(Lfrac)+'Lfrac',nodes,elements,eltype,[[XE_to_plot_list,'elemental',1,'xe optim.']])

