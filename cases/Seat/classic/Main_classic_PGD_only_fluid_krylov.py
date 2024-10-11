import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen,cg
#from numpy.linalg import solve, norm

#import os
#from scipy.linalg import decomp
import pylab as pl
import pickle
import scipy.optimize
import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt
import silex_lib_gmsh
import silex_lib_pgd

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
##############################################################
##############################################################

def NLsystem(X):
    F=X[SolvedDofs_x].copy()
    G=X[SolvedDofs_w].copy()

    #Big_matrix=scipy.sparse.construct.bmat([
    #    [scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  ,  K_sum_fi_Gi_A*fluid_damping-M_sum_fi_Gi_B],
    #    [K_sum_fi_Gi_A.T*fluid_damping-M_sum_fi_Gi_B.T        ,  scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B]] )
    #Second_member = scipy.hstack([  P*scipy.dot(cc,G)  ,  cc*scipy.dot(P,F)  ])
    Big_matrix=scipy.sparse.construct.bmat([
        [scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  ,  K_sum_fi_Gi_A*fluid_damping-M_sum_fi_Gi_B+scipy.tensordot(P,cc,0)],
        [K_sum_fi_Gi_A.T*fluid_damping-M_sum_fi_Gi_B.T+scipy.tensordot(cc,P,0)        ,  scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B]] )
    #Second_member = scipy.hstack([  P*scipy.dot(cc,G)  ,  cc*scipy.dot(P,F)  ])

    #residue = Big_matrix*scipy.sparse.coo_matrix(X).T-scipy.sparse.coo_matrix(Second_member).T

    residue = Big_matrix*scipy.sparse.coo_matrix(X).T

##    Big_matrix_x = scipy.sparse.csc_matrix(  scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  , dtype=mytype)
##    second_member_x = scipy.array(  scipy.dot(cc,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G) , dtype=mytype)
##    F_new = mumps.spsolve( Big_matrix_x , second_member_x , comm=mycomm).T
##    F=F_new.copy()
##
##    Big_matrix_w = scipy.sparse.csc_matrix(  scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B  , dtype=mytype)
##    second_member_w = scipy.array(  cc*scipy.dot(P,F)-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F)+scipy.dot(M_sum_fi_Gi_B.T,F) , dtype=mytype)
##    G_new = mumps.spsolve( Big_matrix_w , second_member_w , comm=mycomm).T
##    G=G_new.copy()


    return scipy.array(residue.todense())

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='geom/sieges_classic_fluid_cavity'
results_file='results/classic_pgd_fluid_cavity'


freq_ini     = 60.0
freq_end     = 80.0

nb_freq_step = 200

# air
celerity=340.0
rho=1.2
#fluid_damping=1.0+0.01j
fluid_damping=1.0

if fluid_damping.imag==0.0:
    mytype='float'
else:
    mytype='c16'

print(mytype)

##############################################################
# PARAMETER mesh: Omega
##############################################################

nb_fcts_PGD = 1

nb_node_w  = nb_freq_step
nodes_w    = scipy.linspace(freq_ini*2.0*scipy.pi, freq_end*2.0*scipy.pi , num=nb_freq_step)
Idnodes_w  = list(range(1,nb_node_w+1,1))
omega_ndof = nb_node_w
nb_elem_w  = nb_node_w-1

SolvedDofW = list(range(omega_ndof))

elements_w = []
for e in range(nb_elem_w):
    elements_w.append([Idnodes_w[e],Idnodes_w[e+1]])

##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes
silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=list(range(fluid_ndof))

##################################################################
# PGD matrices for Omega
##################################################################

IIa,JJa,Va = silex_lib_pgd.frequency_a_matrix(nodes_w,elements_w)
A=scipy.sparse.csc_matrix( (Va,(IIa,JJa)), shape=(omega_ndof,omega_ndof) )
IIb,JJb,Vb = silex_lib_pgd.frequency_b_matrix(nodes_w,elements_w)
B=scipy.sparse.csc_matrix( (Vb,(IIb,JJb)), shape=(omega_ndof,omega_ndof) )

#bb=silex_lib_pgd.frequency_om0_second_member(nodes_w,elements_w)
cc=silex_lib_pgd.frequency_om2_second_member(nodes_w,elements_w)
#cc=silex_lib_pgd.frequency_om0_second_member(nodes_w,elements_w)

##################################################################
# Construct the whole system
##################################################################


K=KFF[SolvedDofF,:][:,SolvedDofF]
M=MFF[SolvedDofF,:][:,SolvedDofF]

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )

P=scipy.zeros((fluid_ndof))
P[0]=1.0

#############################################################################
#       Solve the problem
#############################################################################

ndof_x=fluid_ndof
ndof_w=omega_ndof
SolvedDofs_x=SolvedDofF
SolvedDofs_w=scipy.array(list(range(ndof_w)))
F_i=[]
G_i=[]

#F=scipy.zeros(ndof_x)+1.0
#G=scipy.zeros(ndof_w)+1.0

sum_fi_Gi=scipy.zeros((ndof_x,ndof_w))
niter=0
for i in range(nb_fcts_PGD):
    print("------------------------------")
    print("Recherche du couple : ",i)
    niter   = 0
    critere = 9.9

    if i!=0:
        sum_fi_Gi=sum_fi_Gi+scipy.tensordot(F_i[i-1],G_i[i-1],0)

    K_sum_fi_Gi_A=K*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*A
    M_sum_fi_Gi_B=M*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*B

    # initialize displacement vector
##    F=scipy.array(  scipy.random.random(ndof_x) , dtype=mytype  )
##    G=scipy.array(  scipy.random.random(ndof_w) , dtype=mytype  )
    F=scipy.array(  scipy.zeros(ndof_x)+1.0 , dtype=mytype  )
    G=scipy.array(  scipy.zeros(ndof_w)+1.0 , dtype=mytype  )
    X=scipy.hstack([F,G])

    X_old=X.copy()

##        Big_matrix_w = scipy.sparse.csc_matrix(  scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B  , dtype=mytype)
##        second_member_w = scipy.array(  cc*scipy.dot(P,F)-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F)+scipy.dot(M_sum_fi_Gi_B.T,F) , dtype=mytype)
##        G_new = mumps.spsolve( Big_matrix_w , second_member_w , comm=mycomm).T
##        #G_new = G_new.real
##        #G_new = G_new/scipy.linalg.norm(G_new)
##        #F = F*scipy.linalg.norm(G_new)
##        #G_new[SolvedDofs_w] = scipy.sparse.linalg.spsolve( Big_matrix_w , second_member_w )
##        G=G_new.copy()
##
##        Big_matrix_x = scipy.sparse.csc_matrix(  scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  , dtype=mytype)
##        second_member_x = scipy.array(  scipy.dot(cc,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G) , dtype=mytype)
##        F_new = mumps.spsolve( Big_matrix_x , second_member_x , comm=mycomm).T
##        #F_new[SolvedDofs_x] = scipy.sparse.linalg.spsolve( Big_matrix_x , second_member_x )
##        F=F_new.copy()

##        Big_matrix=scipy.sparse.construct.bmat( [[scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  ,  K_sum_fi_Gi_A*fluid_damping-M_sum_fi_Gi_B],
##                [K_sum_fi_Gi_A.T*fluid_damping-M_sum_fi_Gi_B.T,scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B]
##                                                 ] )
##        Second_member = scipy.hstack([  scipy.dot(cc,G)*P  ,  cc*scipy.dot(P,F)  ])
##        X = mumps.spsolve( Big_matrix , Second_member , comm=mycomm).T
    
##    X_new = scipy.optimize.newton_krylov( NLsystem , X_old, verbose=1)
    X_new = scipy.optimize.anderson( NLsystem , X_old, verbose=1)

    F=X[SolvedDofs_x].copy()
    G=X[SolvedDofs_w].copy()
    
    F_i.append(F.copy())
    G_i.append(G.copy())

##############################################################
# FRF construction
##############################################################
frf=[]
frequencies=[]
press_save=[]
for omi in range(nb_node_w):
    press = scipy.zeros(fluid_ndof, dtype=mytype)
    freq  = nodes_w[omi]/(2.0*scipy.pi)
    frequencies.append(freq)
    for i in range(nb_fcts_PGD):
        press[SolvedDofF]=press[SolvedDofF]+F_i[i][list(range(len(SolvedDofF)))]*G_i[i][omi]

    press_save.append(press)
    #frf.append(silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press))
    frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,press))

frfsave=[frequencies,frf]
f=open(results_file+'.frf','wb')
pickle.dump(frfsave, f)
f.close()

#silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

