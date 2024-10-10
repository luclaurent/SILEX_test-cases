#############################################################################
#      Import libraries
#############################################################################
#      FRF DYNAMIC / WITH DAMPING
#############################################################################
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import sys
sys.path.append('../../librairies')
import silex_lib_gmsh

import pickle

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
MeshFileName='truss'

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_truss_3_damping'

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

# Number of eigenmodes
nbmodes=30

# FRF: frequency range
frequencies=scipy.linspace(1,1000,1000)

# Define material
Young  = 2e11
rho    = 7700.0

# define damping
alpha = 0.0001
beta  = 0.00001

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



# Choose initial "unit" xe
#xe=scipy.ones(nelem)
# or load static optimum xe
f=open('Results_truss_optim_xe_static.pkl','rb')
xe=pickle.load(f)
f.close()



penal      = 3.0
YoungMin   = 2e9
Youngelem  = YoungMin+xe**penal*(Young-YoungMin)
rhoMin     = 77.0
rhoelem    = rhoMin+xe**penal*(rho-rhoMin)


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

Ik,Jk,Vk=silex_lib_elt.stiffnessmatrixoptim(nodes,elements,[Youngelem,Section])
K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof,ndof) )

Im,Jm,Vm=silex_lib_elt.massmatrixoptim(nodes,elements,[rhoelem,Section])
M=scipy.sparse.csc_matrix( (Vm,(Im,Jm)), shape=(ndof,ndof) )

C=alpha*M+beta*K

#############################################################################
#       Eigen value problem
#############################################################################

eigen_values,eigen_vectors= scipy.sparse.linalg.eigsh(K[SolvedDofs,:][:,SolvedDofs],nbmodes,M[SolvedDofs,:][:,SolvedDofs],sigma=0,which='LM')

freq_eigv=list(scipy.sqrt(eigen_values)/(2*scipy.pi))

eigen_vector_list=[]
for i in range(eigen_values.shape[0]):
    Q=scipy.zeros(ndof)
    Q[SolvedDofs]=eigen_vectors[:,i]
    disp=scipy.zeros((nnodes,3))
    disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,2))]
    disp[range(nnodes),2]=scipy.zeros(nnodes)
    eigen_vector_list.append(disp)

silex_lib_gmsh.WriteResults2(ResultsFileName+'_modes',nodes,elements,eltype,[[eigen_vector_list,'nodal',3,'modes']])

#############################################################################
#       Solve the FRF problem
#############################################################################

frf=[]

disp_save=[]
Q=scipy.zeros(ndof, dtype=complex)

for i in range(len(frequencies)):

    freq = frequencies[i]
    omega=2*scipy.pi*freq

    print ("frequency=",freq)

    #sol = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs], F[SolvedDofs].T)
    #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype='d') , scipy.array(F.todense() , dtype='d'), comm=comm_mumps_one_proc()).T
    Q[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]+omega*1j*C[SolvedDofs,:][:,SolvedDofs]-(omega**2)*M[SolvedDofs,:][:,SolvedDofs],dtype='c16') , scipy.array(F[SolvedDofs],dtype='c16'), comm=mycomm).T
    
    #frf.append(scipy.sqrt(Q[(24-1)*2]**2+Q[(24-1)*2+1]**2))
    frf.append(    scipy.sqrt(     (Q[(24-1)*2].real)**2+     (Q[(24-1)*2].imag)**2   )    )

    # displacement written on 2 columns:
    disp=scipy.zeros((nnodes,2))
    disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,2))]
    disp_save.append(disp)
    
frfsave=[frequencies,frf]
silex_lib_gmsh.WriteResults2(ResultsFileName+'_static_optim_with_damping_disp_frf',nodes,elements,eltype,[[disp_save,'nodal',2,'displacement']])

print (" time at the end of the FRF:",time.ctime())

# Save the FRF problem
f=open(ResultsFileName+'_static_optim_damping.frf','wb')
pickle.dump(frfsave, f)
f.close()

print ("structure eigen frequencies : ",freq_eigv)

print ("----- END -----")

