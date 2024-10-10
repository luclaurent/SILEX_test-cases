#############################################################################
#      Import libraries
#############################################################################
#      FRF DYNAMIC / NO DAMPING / WITH OR WITHOUT OPTIMIZATION
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

flag_xe_optim=0

# Input mesh: define the name of the mesh file (*.msh)
if flag_xe_optim==1:
    MeshFileName='Results_truss_static_optim_1_xe_static_optim'
    ResultsFileName='Results_truss_2_XE_optim_NO_damping'
    # Load optimized information
    f=open('Results_truss_static_optim_1_xe_static_optim.pkl','rb')
    [xe,dico]=pickle.load(f)
    f.close()
else:
    MeshFileName='truss'
    ResultsFileName='Results_truss_2_XE_1_NO_damping'
    f=open('Results_truss_static_optim_1_xe_static_optim.pkl','rb')
    [xe,dico]=pickle.load(f)
    f.close()


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
nbmodes=20

# FRF: frequency range
frequencies=scipy.linspace(1,1000,1000)

# Define material
Young  = 2e11
rho    = 7700.0

# define geometry of the cross section
#      
#         Cross section of the rods
Section = 0.05**2 # area of the section
# Inertia
Inertia = 0.05**4/12

# Boundary conditions
if flag_xe_optim==1:
    IdNodesFixed_x=scipy.array([dico[1],dico[19]],dtype=int)
    IdNodesFixed_y=scipy.array([dico[1],dico[19]],dtype=int)
else:
    IdNodesFixed_x=scipy.array([1,7,13,19],dtype=int)
    IdNodesFixed_y=scipy.array([1,7,13,19],dtype=int)

# If the user wants to have only the mesh for gmsh, uncomment next line
#silex_lib_gmsh.WriteResults('maillage_seul',nodes,elements,eltype)

# Load on x direction
LoadX=[]

# Load on y direction
if flag_xe_optim==1:
    LoadY=[[dico[24],-100]
           ]
else:
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

penal      = 3.0
YoungMin   = 2e9
Youngelem  = YoungMin+xe**penal*(Young-YoungMin)
penalrho   = 1.0
rhoMin     = 77.0
rhoelem    = rhoMin+xe**penalrho*(rho-rhoMin)


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

#############################################################################
#       Eigen value problem
#############################################################################

if 1==0:
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

    toc = time.clock()

    silex_lib_gmsh.WriteResults2(ResultsFileName+'_modes',nodes,elements,eltype,[[eigen_vector_list,'nodal',3,'modes']])

#############################################################################
#       Solve the FRF problem
#############################################################################

frf=[]

disp_save=[]

for i in range(len(frequencies)):

    freq = frequencies[i]
    omega=2*scipy.pi*freq

    print ("frequency=",freq)

    #Q[SolvedDofs] = scipy.sparse.linalg.spsolve(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs], F[SolvedDofs].T)
    #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype='d') , scipy.array(F.todense() , dtype='d'), comm=comm_mumps_one_proc()).T
    Q[SolvedDofs] = mumps.spsolve( scipy.sparse.csc_matrix(K[SolvedDofs,:][:,SolvedDofs]-(omega*omega)*M[SolvedDofs,:][:,SolvedDofs],dtype='d') , scipy.array(F[SolvedDofs],dtype='d'), comm=mycomm).T

    if flag_xe_optim==1:
        frf.append(scipy.sqrt( Q[(dico[24]-1)*2]**2 + Q[(dico[24]-1)*2+1]**2) )
    else:
        frf.append(scipy.sqrt( Q[(24-1)*2]**2 + Q[(24-1)*2+1]**2) )
    #frf.append(scipy.absolute(Q[(24-1)*2+1]))

    # displacement written on 2 columns:
    disp=scipy.zeros((nnodes,2))
    disp[range(nnodes),0]=Q[list(range(0,ndof,2))]
    disp[range(nnodes),1]=Q[list(range(1,ndof,2))]
    disp_save.append(disp)
    
frfsave=[frequencies,frf]
silex_lib_gmsh.WriteResults2(ResultsFileName+'_disp_frf',nodes,elements,eltype,[[disp_save,'nodal',2,'displacement']])

print (" time at the end of the FRF:",time.ctime())

#print ("structure eigen frequencies : ",freq_eigv)

# Save the FRF problem
f=open(ResultsFileName+'_static_no_damping.frf','wb')
pickle.dump(frfsave, f)
f.close()

print ("----- END -----")

