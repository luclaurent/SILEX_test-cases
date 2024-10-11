#############################################################################
#      Import libraries
#############################################################################
import string; import time; import scipy
import scipy.sparse
import scipy.sparse.linalg
import mumps

import pickle

import sys
sys.path.append('../../librairies')

#import cProfile

#import silex_lib_tet4_python as silex_lib_elt
import silex_lib_tet4_fortran as silex_lib_elt
import silex_lib_gmsh

import mumps

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
    
mycomm=comm_mumps_one_proc()

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4  python3 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#

#############################################################################
print("SILEX CODE - calcul d'une helice avec des tet4")
#############################################################################

tic = time.clock()
#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################

# Input mesh: define the name of the mesh file (*.msh)
MeshFileName='propeller-tet4'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_propeller_tet4'

# choose the element type
eltype=4

# choose geometry dimension
ndim=3

# choose the results in the results file
flag_write_fields=1

# read the mesh from gmsh
nodes=silex_lib_gmsh.ReadGmshNodes(MeshFileName+'.msh',ndim)
elements,Idnodes=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',eltype,4)

# read surfaces where to impose boundary conditions
elementsS1,IdnodeS1=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,1)
elementsS2,IdnodeS2=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,2)
elementsS3,IdnodeS3=silex_lib_gmsh.ReadGmshElements(MeshFileName+'.msh',2,3)

# write the surface mesh in a gmsh-format file to verify if its correct
#silex_lib_gmsh.WriteResults(ResultsFileName+'Volum',nodes,elements,11)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf1',nodes,elementsS1,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf2',nodes,elementsS2,9)
#silex_lib_gmsh.WriteResults(ResultsFileName+'surf3',nodes,elementsS3,9)

# valeurs pvc trouvee ici: http://www.canplast.ch/pvc_02.html
# Define material
#Young  = 3600 #MPa
#nu     = 0.4

# hetre
Young = 14300.0 #MPa
nu = 0.4
rho = 680e-12 # tonne / mm^3 # 680 kg.m-3 = 680.10^-12 tonne.mm-3 

# Boundary conditions
IdNodesFixed_x=IdnodeS1
IdNodesFixed_y=IdnodeS2
IdNodesFixed_z=IdnodeS1

# compute external forces from pressure
# helice grosso modo du 21'' de diametre (269mm de rayon maxi)
# donc kv de 300 et voltage entre 4s et 6s
#press=0.5*1.2*((300*5*3.7*2*scipy.pi/60.0)*250.0e-3)**2*1e-6  # 1/2 * rho * v^2 avec v = omega R et omega = rpm * 2pi / 60 et rpm = kv * voltage
# calcul pour une traction de 20kg divise par la surface 
#press=(20.0*10.0)/(500.0*35.0)

press = 0.01 #MPa
#F = silex_lib_elt.forceonsurface(nodes,elementsS3,press,[0.0,0.0,0.0])

frequencies=scipy.linspace(1,2000,500)

nbmodes=20

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
F=scipy.zeros(ndof)
F[(1-1)*3+1] = 1.0

#############################################################################
#      compute stiffness matrix
#############################################################################
tic0 = time.clock()
tic = time.clock()
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(nodes,elements,[Young,nu])
toc = time.clock()
print("time to compute the stiffness matrix / FORTRAN:",toc-tic)

K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

tic = time.clock()
Ik,Jk,Vk=silex_lib_elt.massmatrix(nodes,elements,rho)
toc = time.clock()
print("time to compute the mass matrix / FORTRAN:",toc-tic)

M=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) ,dtype=float)

#############################################################################
#       Eigen value problem
#############################################################################

#Q[SolvedDofs] = mumps.spsolve(K[SolvedDofs,:][:,SolvedDofs],F[SolvedDofs])

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],nbmodes,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))

eigen_vector_S_list=[]
for i in range(eigen_values_S.shape[0]):
    Q=scipy.zeros(ndof)
    Q[SolvedDofs]=eigen_vectors_S[:,i]
    disp=scipy.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))]
    eigen_vector_S_list.append(disp)

toc = time.clock()
print ("structure eigen frequencies : ",freq_eigv_S)
print ("time for computing the structure modes:",toc-tic)


#############################################################################
#         Write results to gmsh format  : EIGEN-MODES
#############################################################################

silex_lib_gmsh.WriteResults2(ResultsFileName+'_structure_modes',nodes,elements,eltype,[[eigen_vector_S_list,'nodal',3,'modes']])

##############################################################
# FRF computation 
##############################################################

frf=[]

disp_save=[]

for i in range(len(frequencies)):

    freq = frequencies[i]
    omega=2*scipy.pi*freq

    print ("frequency=",freq)

    #sol = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs], F[SolvedDofs].T)
    #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype='d') , scipy.array(F.todense() , dtype='d'), comm=comm_mumps_one_proc()).T
    Q[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype='d') , scipy.array(F[SolvedDofs],dtype='d'), comm=mycomm).T

    frf.append(scipy.sqrt(Q[(1-1)*3]**2+Q[(1-1)*3+1]**2+Q[(1-1)*3+2]**2))
        
    disp=scipy.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,3))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,3))]
    disp[range(nnodes),2]=Q[list(range(2,ndof,3))]
    disp_save.append(disp)
    
frfsave=[frequencies,frf]
silex_lib_gmsh.WriteResults2(ResultsFileName+'_disp_frf',nodes,elements,eltype,[[disp_save,'nodal',3,'displacement']])

print (" time at the end of the FRF:",time.ctime())

# Save the FRF problem
f=open(ResultsFileName+'_no_damping.frf','wb')
pickle.dump(frfsave, f)
f.close()

print ("----- END -----")

