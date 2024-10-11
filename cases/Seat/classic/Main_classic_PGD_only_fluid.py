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
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='geom/sieges_classic_fluid_cavity'
results_file='results/classic_pgd_fluid_cavity'


freq_ini     = 20.0
freq_end     = 70.0

nb_freq_step = 40

# air
celerity=340.0
rho=1.2
fluid_damping=1.0+0.002j
#fluid_damping=1.0

if fluid_damping.imag==0.0:
    mytype='float'
else:
    mytype='c16'

print("type de resolution: ",mytype)

##############################################################
# PARAMETER mesh: Omega
##############################################################

#nb_fcts_PGD = 19

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
fluid_elements_S2,IdNodesS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,2)

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

bb=silex_lib_pgd.frequency_om0_second_member(nodes_w,elements_w)
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

#P=scipy.zeros((fluid_ndof))
#P[0]=1.0

P = silex_lib_xfem_acou_tet4.forceonsurface(fluid_nodes,fluid_elements_S2,1.0)

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

#sum_fi_Gi=scipy.zeros((ndof_x,ndof_w))

alpha_i = scipy.array(  scipy.zeros(5) , dtype=mytype  )

##F=scipy.array(  scipy.random.random(ndof_x) , dtype=mytype  )
##G=scipy.array(  scipy.random.random(ndof_w) , dtype=mytype  )

niter=0
#for i in range(nb_fcts_PGD):

nb_fcts_PGD=0
residu=99.9
while (residu>1e-2):
    i=nb_fcts_PGD
    print("------------------------------")
    print("Recherche du couple : ",i)
    niter   = 0
    critere = 9.9

    #if i!=0:
        #sum_fi_Gi=sum_fi_Gi+scipy.tensordot(F_i[i-1],G_i[i-1],0)
    sum_fi_Gi=scipy.zeros((ndof_x,ndof_w))
    for k in range(i):
        sum_fi_Gi=sum_fi_Gi+scipy.tensordot(F_i[k],G_i[k],0)*alpha_i[k]


    K_sum_fi_Gi_A=K*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*A
    M_sum_fi_Gi_B=M*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*B
    #print("hello 1")
    # initialize displacement vector
##    F=scipy.array(  scipy.random.random(ndof_x) , dtype=mytype  )
##    G=scipy.array(  scipy.random.random(ndof_w) , dtype=mytype  )
##    F=scipy.array(  scipy.random.random(ndof_x)+scipy.random.random(ndof_x)*1j , dtype='c16'  )
##    G=scipy.array(  scipy.random.random(ndof_w)+scipy.random.random(ndof_w)*1j , dtype='c16'  )
    #X=scipy.hstack([F,G])
##    if i==0:
##        F=scipy.zeros(ndof_x)+1.0
##        G=scipy.zeros(ndof_w)+1.0
##    F=F/scipy.linalg.norm(F)
##    G=G/scipy.linalg.norm(G)
    
##    if i!=0:
##        #F=F_i[i-1].copy()+scipy.random.random(ndof_x)*max(F_i[i-1])*1e-1
##        #G=G_i[i-1].copy()+scipy.random.random(ndof_w)*max(G_i[i-1])*1e-1
##        F=F_i[i-1].copy()+scipy.random.random(ndof_x)*max(F_i[i-1])*1e-1
##        G=G_i[i-1].copy()+scipy.random.random(ndof_w)*max(G_i[i-1])*1e-1

    F=scipy.array(  scipy.zeros(ndof_x)+1.0 , dtype=mytype  )
    G=scipy.array(  scipy.zeros(ndof_w)+1.0 , dtype=mytype  )


    while (critere>3e-2):
        F_new=scipy.zeros(ndof_x, dtype=mytype)
        G_new=scipy.zeros(ndof_w, dtype=mytype)

        G_old=G.copy()
        F_old=F.copy()

        #X_old=X.copy()
        

        Big_matrix_w = scipy.sparse.csc_matrix(  scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B  , dtype=mytype)
        second_member_w = scipy.array(  cc*scipy.dot(P,F)-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F)+scipy.dot(M_sum_fi_Gi_B.T,F) , dtype=mytype)
        G_new = mumps.spsolve( Big_matrix_w , second_member_w , comm=mycomm).T
        #G_new = G_new.real
        #G_new = G_new/scipy.linalg.norm(G_new)
        #F = F*scipy.linalg.norm(G_new)
        #G_new[SolvedDofs_w] = scipy.sparse.linalg.spsolve( Big_matrix_w , second_member_w )
        G=G_new.copy()

        Big_matrix_x = scipy.sparse.csc_matrix(  scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M  , dtype=mytype)
        second_member_x = scipy.array(  scipy.dot(cc,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G) , dtype=mytype)
        F_new = mumps.spsolve( Big_matrix_x , second_member_x , comm=mycomm).T
        #F_new[SolvedDofs_x] = scipy.sparse.linalg.spsolve( Big_matrix_x , second_member_x )
        F=F_new.copy()
        
        #critere = max( max( abs( (G_old-G)/G) ) , max( abs( (F_old-F)/F ) ) )
        critere = max( scipy.linalg.norm(G_old-G)/scipy.linalg.norm(G_old+G) , scipy.linalg.norm(F_old-F)/scipy.linalg.norm(F_old+F) )
        #critere = scipy.linalg.norm(F_old-F)/scipy.linalg.norm(F)

        #F=0.5*F.copy()+0.5*F_old
        #G=0.5*G.copy()+0.5*G_old
        
        niter=niter+1
        if niter%15==0:
            print("forcage!")
            critere=0.0
##            print("re-initialisation")
##            F=scipy.array(  scipy.random.random(ndof_x) , dtype=mytype  )
##            G=scipy.array(  scipy.random.random(ndof_w) , dtype=mytype  )
  
        print("critere = ",critere,      "niter = ",niter)

    second_member_w = scipy.array(  cc*scipy.dot(P,F)-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F)+scipy.dot(M_sum_fi_Gi_B.T,F) , dtype=mytype)
    second_member_x = scipy.array(  scipy.dot(cc,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G) , dtype=mytype)
##    residu=max( scipy.linalg.norm(second_member_w)/scipy.linalg.norm(cc*scipy.dot(P,F))
##                , scipy.linalg.norm(second_member_x)/scipy.linalg.norm(scipy.dot(cc,G)*P) )
##    residu=max( scipy.linalg.norm(second_member_w)/(scipy.linalg.norm(cc*scipy.dot(P,F))
##                                                    +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F))
##                                                    +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B.T,F)) )
##                , scipy.linalg.norm(second_member_x)/(scipy.linalg.norm(scipy.dot(cc,G)*P)
##                                                       +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping)
##                                                       +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B,G)) )
##                )
    residu=max( scipy.linalg.norm(second_member_w)/(scipy.linalg.norm(cc*scipy.dot(P,F))
                                                    +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F))
                                                    +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B.T,F)) )
                , scipy.linalg.norm(second_member_x)/(scipy.linalg.norm(scipy.dot(cc,G)*P)
                                                       +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping)
                                                       +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B,G)) )
                )
    print("residu global = ",residu)

    F=F.copy()/scipy.linalg.norm(F.copy()) # normalisation
    G=G.copy()/scipy.linalg.norm(G.copy())

    F_i.append(F.copy())
    G_i.append(G.copy())

    # Projection on the new basis
    Ktilde=scipy.zeros((i+1,i+1), dtype=mytype)
    Mtilde=scipy.zeros((i+1,i+1), dtype=mytype)
    second_member_alpha = scipy.zeros( i+1 , dtype=mytype)

    for l in range(i+1):
        second_member_alpha[l] = scipy.dot(cc,G_i[l])*scipy.dot(P,F_i[l])
        for m in range(i+1):
            Ktilde[l,m]=scipy.dot(F_i[l],K*F_i[m])*scipy.dot(G_i[l],A*G_i[m])*fluid_damping
            Mtilde[l,m]=scipy.dot(F_i[l],M*F_i[m])*scipy.dot(G_i[l],B*G_i[m])

    alpha_i = mumps.spsolve( scipy.sparse.csc_matrix( Ktilde-Mtilde ) , second_member_alpha , comm=mycomm).T
    #print("     Solution alpha_i:",alpha_i)

    nb_fcts_PGD=nb_fcts_PGD+1

##    # calcul du residu
##    residu_x = scipy.array(  scipy.zeros(ndof_x) , dtype=mytype )
##    residu_w = scipy.array(  scipy.zeros(ndof_w) , dtype=mytype )
##
##    for k in range(i+1):
##        residu_w = residu_w+scipy.dot( scipy.dot(F,K*F)*A*fluid_damping-scipy.dot(F,M*F)*B , )
    
# Save the space modes

F_list=[]
for i in range(nb_fcts_PGD):
    tmp=F_i[i].real
    F_list.append(tmp)

silex_lib_gmsh.WriteResults2(results_file+'_fluid_space_modes',fluid_nodes,fluid_elements,4,[[F_list,'nodal',1,'pressure']])

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
        press[SolvedDofF]=press[SolvedDofF]+F_i[i][list(range(len(SolvedDofF)))]*G_i[i][omi]*alpha_i[i]

    press_save.append(press.real)
    #frf.append(silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press))
    frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,press))

frfsave=[frequencies,frf]
f=open(results_file+'.frf','wb')
pickle.dump(frfsave, f)
f.close()

silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

##############################################################
import pylab as pl
prefsquare=20e-6*20e-6

pl.figure(1)
for i in range(nb_fcts_PGD):
    pl.plot(frequencies,scipy.array(G_i[i]),label='G'+str(i), linewidth=2)

pl.show()
