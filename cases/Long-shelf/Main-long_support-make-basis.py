#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg

import sys
sys.path.append('../../librairies')

import silex_lib_tet4_fortran as silex_lib_elt
import silex_lib_gmsh
import pickle

#############################################################################
print('SILEX CODE - calcul d un support avec des tet4')
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='long_support'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_long_support_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,10)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)
elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,4)
elementsS5,IdnodeS5=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,5)
elementsS6,IdnodeS6=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,6)
elementsS7,IdnodeS7=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,7)

# write the surface mesh in a gmsh-format file to verify if its correct
silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,4)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf4',nodes,elementsS4,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf5',nodes,elementsS5,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf6',nodes,elementsS6,2)
silex_lib_gmsh.WriteResults(ResultsFileName+'surf7',nodes,elementsS7,2)

# Define material
Young  = 74000.0
nu     = 0.33

# Boundary conditions
IdNodesFixed_x=IdnodeS7
IdNodesFixed_y=IdnodeS7
IdNodesFixed_z=IdnodeS6

# compute external forces from pressure
press = 1.0/(10.0*2.0*scipy.pi*7.0) #MPa
# give the direction of the surfacic load:
#          if [0.0,0.0,0.0] then the local normal to the surface is used
#          otherwise, the direction is normalized to 1
direction = [1.0,0.0,0.0]
F1x = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2x = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3x = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4x = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5x = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)
direction = [0.0,1.0,0.0]
F1y = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2y = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3y = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4y = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5y = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)
direction = [0.0,0.0,1.0]
F1z = silex_lib_elt.forceonsurface(nodes,elementsS1,press,direction)
F2z = silex_lib_elt.forceonsurface(nodes,elementsS2,press,direction)
F3z = silex_lib_elt.forceonsurface(nodes,elementsS3,press,direction)
F4z = silex_lib_elt.forceonsurface(nodes,elementsS4,press,direction)
F5z = silex_lib_elt.forceonsurface(nodes,elementsS5,press,direction)

Fbasis=[F1x,F1y,F1z,F2x,F2y,F2z,F3x,F3y,F3z,F4x,F4y,F4z,F5x,F5y,F5z]

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
MySolve = scipy.sparse.linalg.factorized(K[SolvedDofs,:][:,SolvedDofs]) # Makes LU decomposition

#############################################################################
#       Solve the problem : find the 15 vector basis
#############################################################################
print("COMPUTE BASIS Q")
Qbasis=[]
i=1
for F in Fbasis:
    tic = time.clock()
    print('basis ',i)
    Q=scipy.zeros(ndof)
    Q[SolvedDofs] = MySolve( F[SolvedDofs] )
    Qbasis.append(Q)
    toc = time.clock()
    print("time to solve one basis problem:",toc-tic)
    i=i+1

print("---------------------")

#############################################################################
#       compute smooth stress and error in elements
#############################################################################
print("COMPUTE STRESS")

Sigmabasis=[]
i=1
for Q in Qbasis:
    tic = time.clock()
    SigmaNodes=scipy.zeros((nnodes))
    SigmaElem,SigmaNodes,EpsilonElem,EpsilonNodes,ErrorElem,ErrorGlobal=silex_lib_elt.compute_stress_strain_error(nodes,elements,[Young,nu],Q)
    Sigmabasis.append(SigmaNodes)
    toc = time.clock()
    print('basis ',i)
    print("time to compute stres and error:",toc-tic)
    print("The global error is:",ErrorGlobal)
    i=i+1

print("---------------------")

#############################################################################
#         Write results to gmsh format
#############################################################################
print("WRITE RESULTS")
for i in range(len(Qbasis)):
    tic = time.clock()
    print("basis ",i+1)
    # displacement written on 3 columns:
    Q=Qbasis[i]
    disp=scipy.zeros((nnodes,ndim))
    disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))]

    # external load written on 3 columns:
    F=Fbasis[i]
    load=scipy.zeros((nnodes,ndim))
    load[range(nnodes),0]=F[list(range(0,ndof,3))]
    load[range(nnodes),1]=F[list(range(1,ndof,3))]
    load[range(nnodes),2]=F[list(range(2,ndof,3))]

    SigmaNodes=Sigmabasis[i]

    fields_to_write=[ [disp,'nodal',ndim,'displacement'],
                      [SigmaNodes[range(nnodes),[6]],'nodal',1,'Smooth Sigma V.M.'],
                      [load,'nodal',ndim,'Force'],
                      ]

    # write the mesh and the results in a gmsh-format file
    silex_lib_gmsh.WriteResults(ResultsFileName+'_'+str(i+1),nodes,elements,eltype,fields_to_write)

    toc = time.clock()
    print("time to write results:",toc-tic)

print("WRITE BASIS")
f=open(ResultsFileName+'_basis','wb')
pickle.dump([Qbasis,Sigmabasis], f)
f.close()

print("----- END -----")


