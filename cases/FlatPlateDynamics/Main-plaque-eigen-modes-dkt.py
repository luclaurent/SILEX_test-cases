#############################################################################
#      Import libraries
#############################################################################
import string
import time
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from pathlib import Path

#import cProfile

from SILEXlib import silex_lib_dkt as silex_lib_elt
from SILEXlib import silex_lib_gmsh

#############################################################################
print("SILEX CODE - flat plate - Eigen modes")
#############################################################################

tic = time.perf_counter()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
#MeshFileName='plaque-dkt'
MeshFileName='flat_plate_dkt'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_flat_plate_dkt_eigenmodes'
ResultsFullPath = Path(__file__).parent / ResultsFileName

# choose the element type
eltype=2


# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=0

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(str(Path(__file__).parent /  (MeshFileName+'.msh')),ndim)

# surface
#elements,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,3)
elements,IdnodeS=silex_lib_gmsh.ReadGmshElements(str(Path(__file__).parent /  (MeshFileName+'.msh')),eltype,1)

# bords
#elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,1)
elementsL2,IdnodeL2=silex_lib_gmsh.ReadGmshElements(str(Path(__file__).parent /  (MeshFileName+'.msh')),1,2)

# force
#elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,2)
#elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(str(Path(__file__).parent /  (MeshFileName+'.msh')),1,2)

# appui sur y
#elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',1,4)
#elementsS4,IdnodeS4=silex_lib_gmsh.ReadGmshElements(str(Path(__file__).parent /  (MeshFileName+'.msh')),1,4)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'_1',nodes,elementsS1,1)
silex_lib_gmsh.WriteResults(str(Path(__file__).parent /  (ResultsFileName+'_S')),nodes,elements,eltype)
#silex_lib_gmsh.WriteResults(ResultsFileName+'_2',nodes,elementsS2,1)
silex_lib_gmsh.WriteResults(str(Path(__file__).parent /  (ResultsFileName+'_L2')),nodes,elementsL2,1)


# Define material
Young     = 200000e6
nu        = 0.3
thickness = 0.003
density   = 7000

# Boundary conditions
IdNodesFixed_x=[1,4]
IdNodesFixed_y=[1]
#IdNodesFixed_x=IdnodeL2
#IdNodesFixed_y=IdnodeL2
IdNodesFixed_z=IdnodeL2
IdNodesFixed_rotx=[]
IdNodesFixed_roty=[]
IdNodesFixed_rotz=IdnodeL2

# number of eigen modes
nbmodes=20

# external load
fx=0.0
fy=0.0
fz=-1000.0
#F=silex_lib_elt.forceonline(nodes,elementsS2,[forcex/100.0,0.0,forcez/100.0,forcex/100.0,0.0,forcez/100.0],[1000.0,0.0,0.0,1000.0,100.0,0.0])

toc = time.process_time()
print("time for the user part: {}".format(toc-tic))
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

# Define external load
F=np.zeros(ndof)   
F[(10-1)*6]=fx
F[(10-1)*6+1]=fy
F[(10-1)*6+2]=fz

# define fixed dof
Fixed_Dofs = np.hstack([
    (np.array(IdNodesFixed_x)-1)*6,
    (np.array(IdNodesFixed_y)-1)*6+1,
    (np.array(IdNodesFixed_z)-1)*6+2,
    (np.array(IdNodesFixed_rotx)-1)*6+3,
    (np.array(IdNodesFixed_roty)-1)*6+4,
    (np.array(IdNodesFixed_rotz)-1)*6+5])

# define free dof
SolvedDofs = np.setdiff1d(range(ndof),Fixed_Dofs)

# initialize displacement vector
Q=np.zeros(ndof)

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.process_time()
tic = time.process_time()
#print (silex_lib_elt.stiffnessmatrix.__doc__)
Ik,Jk,Vk,Vm=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu,thickness,density])
toc = time.process_time()
print("time to compute the stiffness matrix / FORTRAN: {}".format(toc-tic))

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)
M=scipy.sparse.csc_matrix( (Vm,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Solve the problem - EIGEN MODES
#############################################################################

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],nbmodes,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

freq_eigv_S=list(np.sqrt(eigen_values_S)/(2*np.pi))

eigen_vector_S_list=[]
for i in range(eigen_values_S.shape[0]):
    Q=np.zeros(ndof)
    Q[SolvedDofs]=eigen_vectors_S[:,i]
    disp=np.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,6))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,6))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,6))]
    eigen_vector_S_list.append(disp)

toc = time.process_time()
print ("structure eigen frequencies : ",freq_eigv_S)
print ("time for computing the structure modes: {}".format(toc-tic))

#############################################################################
#         Write results to gmsh format
#############################################################################
silex_lib_gmsh.WriteResults2(str(ResultsFullPath),nodes,elements,eltype,[[eigen_vector_S_list,'nodal',3,'modes']])

toc = time.process_time()
print("time to write results:  {}".format(toc-tic))
print("----- END -----")



