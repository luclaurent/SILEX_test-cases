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

import silex_lib_rod as silex_lib_elt
import pylab as pl
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

#############################################################################
print("SILEX CODE - analyse frequentielle classique d'une barre en traction")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.process_time()

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_modal_classic'

L    = 1.0 #m
young  = 200000.0e6
section = 0.01**2 
k    = young*section
rho  = 7500.0
w1   = 100.0*2.0*np.pi
w2   = 10000.0*2.0*np.pi
load = 100.0
nb_node_x = 100
nb_elem_x = nb_node_x-1

nb_node_w = 500
nb_elem_w = nb_node_w-1


damping=1e-6j
#damping=1.0

if damping.imag==0.0:
    mytype='float'
else:
    mytype='c16'

# SPACE mesh
nodes_x=scipy.linspace(0.0, L , num=nb_node_x)
Idnodes_x=list(range(1,nb_node_x+1,1))

elements_x=[]
for e in range(nb_elem_x):
    elements_x.append([Idnodes_x[e],Idnodes_x[e+1]])
    
# PARAMETER mesh: Omega
nodes_w=scipy.linspace(w1, w2 , num=nb_node_w)

# SPACE : Boundary conditions
IdNodesFixed_x=np.array([1],dtype=int)
# initialize force vector
P=np.zeros((nb_node_x))
P[nb_node_x-1]=load


#############################################################################
#      expert part: define load and prescribed displacements
#############################################################################
#      initialisations
#############################################################################

ndof_x=nb_node_x

# define fixed dof
Fixed_Dofs_x = scipy.hstack([(IdNodesFixed_x-1)*1])

# define free dof
SolvedDofs_x = np.setdiff1d(range(ndof_x),Fixed_Dofs_x)

# initialize displacement vector
U=np.zeros(ndof_x,dtype=mytype)

#############################################################################
#      compute matrices
#############################################################################

# SPACE : stiffness matrix
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(np.array(nodes_x),np.array(elements_x),[young,section])
K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_x,ndof_x) )

# SPACE : mass matrix
Ik,Jk,Vk=silex_lib_elt.massmatrix(np.array(nodes_x),np.array(elements_x),rho*section)
M=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_x,ndof_x) )

# SPACE : damping matrix
D=damping*K

#############################################################################
#       Solve the problem
#############################################################################

sol=[]
for omega in nodes_w:
    
        Big_matrix_x = K[SolvedDofs_x,:][:,SolvedDofs_x]-omega**2*M[SolvedDofs_x,:][:,SolvedDofs_x]+omega*D[SolvedDofs_x,:][:,SolvedDofs_x]
        #U[SolvedDofs_x] = scipy.sparse.linalg.spsolve(Big_matrix_x,P[SolvedDofs_x])
        U[SolvedDofs_x] = mumps.spsolve(Big_matrix_x,np.array(P[SolvedDofs_x],dtype=mytype) , comm=mycomm)
        sol.append(U[nb_node_x-1])

u_ex=[]
for i in range(10):
    print("Freq. analytique",((2*(i+1)-1)*np.pi*scipy.sqrt(young/rho))/(2.0*L*2.0*np.pi))
    u_ex.append(np.sin((2*(i+1)-1)*np.pi*nodes_x/(2.0*L)))


f=open(ResultsFileName+'_frf','wb')
pickle.dump([sol,nodes_w], f)
f.close()

pl.figure(3)
pl.plot(nodes_w/(2.0*np.pi),np.log(sol),label='dep. x=L', linewidth=2)
pl.xlabel('freq.')
pl.ylabel('log(dep.)')
pl.grid('on')
pl.legend(loc=1)

pl.show()

