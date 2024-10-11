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

import silex_lib_pgd_fortran as silex_lib_pgd

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
print("SILEX CODE - analyse frequentielle PGD d'une barre en traction")
#############################################################################

#############################################################################
#      USER PART: Import mesh, boundary conditions and material
#############################################################################
tic = time.clock()

# Output result file: define the name of the result file (*.msh)
ResultsFileName='Results_modal_pgd'

L    = 1.0 #m
young  = 200000.0e6
section = 0.01**2 
k    = young*section
rho  = 7500.0
w1   = 100.0*2.0*scipy.pi
w2   = 10000.0*2.0*scipy.pi
load = 100.0
nb_node_x = 100
nb_elem_x = nb_node_x-1

#nb_fcts_PGD = 15

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
Idnodes_w=list(range(1,nb_node_w+1,1))

elements_w=[]
for e in range(nb_elem_w):
    elements_w.append([Idnodes_w[e],Idnodes_w[e+1]])

#silex_lib_gmsh.WriteResults('toto',nodes,elements,1)

# SPACE : Boundary conditions
IdNodesFixed_x=scipy.array([1],dtype=int)
# initialize force vector
P=scipy.zeros((nb_node_x))
P[nb_node_x-1]=load

mytype='float'

#############################################################################
#      expert part: define load and prescribed displacements
#############################################################################
#      initialisations
#############################################################################

ndof_x=nb_node_x
ndof_w=nb_node_w

# define fixed dof
Fixed_Dofs_x = scipy.hstack([(IdNodesFixed_x-1)*1])

# define free dof
SolvedDofs_x = scipy.setdiff1d(range(ndof_x),Fixed_Dofs_x)

SolvedDofs_w = list(range(nb_node_w))

# initialize displacement vector
F=scipy.zeros(ndof_x)+1.0
G=scipy.zeros(ndof_w)+1.0

#F=nodes_x*load
#G=1.0/nodes_k

#############################################################################
#      compute matrices
#############################################################################

# SPACE : stiffness matrix
Ik,Jk,Vk=silex_lib_elt.stiffnessmatrix(scipy.array(nodes_x),scipy.array(elements_x),[young,section])
K=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_x,ndof_x) )

# SPACE : mass matrix
Ik,Jk,Vk=silex_lib_elt.massmatrix(scipy.array(nodes_x),scipy.array(elements_x),rho*section)
M=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_x,ndof_x) )

# PARAMETER matrix A
Ik,Jk,Vk=silex_lib_pgd.frequency_a_matrix(scipy.array(nodes_w),scipy.array(elements_w))
A=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_w,ndof_w) )

# PARAMETER matrix B
Ik,Jk,Vk=silex_lib_pgd.frequency_b_matrix(scipy.array(nodes_w),scipy.array(elements_w))
B=scipy.sparse.csc_matrix( (Vk,(Ik,Jk)), shape=(ndof_w,ndof_w) )

# PARAMETER second member b
b=silex_lib_pgd.frequency_om0_second_member(scipy.array(nodes_w),scipy.array(elements_w))


#############################################################################
#       Solve the problem
#############################################################################


F_i=[]
G_i=[]

alpha_i = scipy.array(  scipy.zeros(5) , dtype=mytype  )
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

    #for j in range(i):
        #print j
        #for k in range(ndof_x):
            #for l in range(ndof_w):
                #if i==2:
                #    stop
                #sum_fi_GiT[k,l]=sum_fi_GiT[k,l]+F_i[j][k]*G_i[j][l]
        #sum_fi_GiT=sum_fi_GiT+scipy.tensordot(F_i[j],G_i[j],0)
    #if i!=0:
    #    sum_fi_GiT=sum_fi_GiT+scipy.tensordot(F_i[i-1],G_i[i-1],0)

    sum_fi_Gi=scipy.zeros((ndof_x,ndof_w))
    for k in range(i):
        sum_fi_Gi=sum_fi_Gi+scipy.tensordot(F_i[k],G_i[k],0)*alpha_i[k]

    K_sum_fi_Gi_A=K[SolvedDofs_x,:][:,SolvedDofs_x]*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*A[SolvedDofs_w,:][:,SolvedDofs_w]
    M_sum_fi_Gi_B=M[SolvedDofs_x,:][:,SolvedDofs_x]*sum_fi_Gi[SolvedDofs_x,:][:,SolvedDofs_w]*B[SolvedDofs_w,:][:,SolvedDofs_w]
    # initialize displacement vector
    F=scipy.random.random(ndof_x)
    G=scipy.random.random(ndof_w)
    #F=scipy.array(  scipy.zeros(ndof_x)+1.0 , dtype=mytype  )
    #G=scipy.array(  scipy.zeros(ndof_w)+1.0 , dtype=mytype  )

    while (critere>1e-2):
        F_new=scipy.zeros(ndof_x, dtype=mytype)
        G_new=scipy.zeros(ndof_w, dtype=mytype)
        G_old=G.copy()
        F_old=F.copy()

        Big_matrix_w = scipy.dot(F[SolvedDofs_x].T,K[SolvedDofs_x,:][:,SolvedDofs_x]*F[SolvedDofs_x])*A[SolvedDofs_w,:][:,SolvedDofs_w]
        Big_matrix_w = Big_matrix_w-scipy.dot(F[SolvedDofs_x].T,M[SolvedDofs_x,:][:,SolvedDofs_x]*F[SolvedDofs_x])*B[SolvedDofs_w,:][:,SolvedDofs_w]
        second_member_w = b*scipy.dot(P[SolvedDofs_x],F[SolvedDofs_x])-scipy.dot(K_sum_fi_Gi_A.T,F[SolvedDofs_x])+scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])
        G_new[SolvedDofs_w] = mumps.spsolve(Big_matrix_w,second_member_w , comm=mycomm).T
        G=G_new.copy()

        Big_matrix_x = scipy.dot(G[SolvedDofs_w].T,A[SolvedDofs_w,:][:,SolvedDofs_w]*G[SolvedDofs_w])*K[SolvedDofs_x,:][:,SolvedDofs_x]
        Big_matrix_x = Big_matrix_x-scipy.dot(G[SolvedDofs_w].T,B[SolvedDofs_w,:][:,SolvedDofs_w]*G[SolvedDofs_w])*M[SolvedDofs_x,:][:,SolvedDofs_x]
        second_member_x = scipy.dot(b,G)*P[SolvedDofs_x]-scipy.dot(K_sum_fi_Gi_A,G[SolvedDofs_w])+scipy.dot(M_sum_fi_Gi_B,G[SolvedDofs_w])
        F_new[SolvedDofs_x] = mumps.spsolve(Big_matrix_x,second_member_x , comm=mycomm).T
        F=F_new.copy()
        #G=G_new.copy()

        #critere = scipy.dot(G_old-G_new,G_old-G_new)**2+scipy.dot(F-F_new,F-F_new)**2
        critere = max( scipy.linalg.norm(G_old-G)/scipy.linalg.norm(G_old+G) , scipy.linalg.norm(F_old-F)/scipy.linalg.norm(F_old+F) )
        niter=niter+1

        if niter%15==0:
            print("forcage!")
            critere=0.0

   
        print("critere = ",critere,      "niter = ",niter)
        #c=input('?')

    second_member_w = scipy.array(  b*scipy.dot(P[SolvedDofs_x],F[SolvedDofs_x])-scipy.dot(K_sum_fi_Gi_A.T,F[SolvedDofs_x])+scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x]) , dtype=mytype)
    second_member_x = scipy.array(  scipy.dot(b,G[SolvedDofs_w])*P[SolvedDofs_x]-scipy.dot(K_sum_fi_Gi_A,G[SolvedDofs_w])+scipy.dot(M_sum_fi_Gi_B,G[SolvedDofs_w]) , dtype=mytype)
    residu=max( scipy.linalg.norm(second_member_w)/(scipy.linalg.norm(b*scipy.dot(P,F))
                                                    +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A.T,F[SolvedDofs_x]))
                                                    +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])) )
                , scipy.linalg.norm(second_member_x)/(scipy.linalg.norm(scipy.dot(b,G[SolvedDofs_w])*P[SolvedDofs_x])
                                                       +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A,G[SolvedDofs_w]))
                                                       +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B,G[SolvedDofs_w])) )
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
        second_member_alpha[l] = scipy.dot(b,G_i[l])*scipy.dot(P[SolvedDofs_x],F_i[l][SolvedDofs_x])
        for m in range(i+1):
            Ktilde[l,m]=scipy.dot(F_i[l][SolvedDofs_x],K[SolvedDofs_x,:][:,SolvedDofs_x]*F_i[m][SolvedDofs_x])*scipy.dot(G_i[l],A*G_i[m])
            Mtilde[l,m]=scipy.dot(F_i[l][SolvedDofs_x],M[SolvedDofs_x,:][:,SolvedDofs_x]*F_i[m][SolvedDofs_x])*scipy.dot(G_i[l],B*G_i[m])

    alpha_i = mumps.spsolve( scipy.sparse.csc_matrix( Ktilde-Mtilde ) , second_member_alpha , comm=mycomm).T

    nb_fcts_PGD=nb_fcts_PGD+1

    # plot the solution
    sol=scipy.zeros(ndof_w)
    for i in range(nb_fcts_PGD):
        sol=sol+F_i[i][nb_node_x-1]*G_i[i]*alpha_i[i]
    
##    pl.figure(1)
##    pl.plot(nodes_x,F/max(F),label='F', linewidth=2)
##    pl.xlabel('x')
##    pl.ylabel('Displacement (m)')
##    pl.grid('on')
##    pl.legend(loc=4)
##
##    pl.figure(2)
##    pl.plot(nodes_w/(2.0*scipy.pi),scipy.log(G),label='G', linewidth=2)
##    pl.xlabel('freq.')
##    pl.ylabel('log(G(w))')
##    pl.grid('on')
##    pl.legend(loc=1)
##
##    pl.figure(3)
##    pl.plot(nodes_w/(2.0*scipy.pi),scipy.log(sol),label='dep. x=L', linewidth=2)
##    pl.xlabel('freq.')
##    pl.ylabel('log(dep.)')
##    pl.grid('on')
##    pl.legend(loc=1)
##    #pl.ion()
##    #pl.show()
##    pl.draw()
##    #pl.pause(0.0001)


u_ex=[]
for i in range(10):
    print("Freq. analytique",((2*(i+1)-1)*scipy.pi*scipy.sqrt(young/rho))/(2.0*L*2.0*scipy.pi))
    u_ex.append(scipy.sin((2*(i+1)-1)*scipy.pi*nodes_x/(2.0*L)))

sol=scipy.zeros(ndof_w)
for i in range(nb_fcts_PGD):
    sol=sol+F_i[i][nb_node_x-1]*G_i[i]*alpha_i[i]

f=open(ResultsFileName+'_frf','wb')
pickle.dump([sol,nodes_w], f)
f.close()

### Save the space modes
##
##F_list=[]
##for i in range(nb_fcts_PGD):
##    tmp=F_i[i].real
##    F_list.append(tmp)
##
##silex_lib_gmsh.WriteResults2(results_file+'_space_modes',fluid_nodes,fluid_elements,4,[[F_list,'nodal',1,'pressure']])



pl.figure(1)
for i in range(nb_fcts_PGD):
    pl.plot(nodes_x,F_i[i]/max(F_i[i]),label='F_'+str(i), linewidth=2)
    #pl.plot(nodes_x,u_ex[8-1]*max(F))
pl.xlabel('x')
pl.ylabel('Displacement (m)')
pl.grid('on')
pl.legend(loc=4)

pl.figure(2)
for i in range(nb_fcts_PGD):
    pl.plot(nodes_w/(2.0*scipy.pi),scipy.log(G_i[i]),label='G_'+str(i), linewidth=2)
pl.xlabel('freq.')
pl.ylabel('log(G(w))')
pl.grid('on')
pl.legend(loc=1)

f=open('Results_modal_classic_frf','rb')
sol_ref,nodes_w_ref=pickle.load(f)
f.close()

pl.figure(3)
pl.plot(nodes_w_ref/(2.0*scipy.pi),scipy.log(sol_ref),label='dep. x=L; ref.', linewidth=2)
pl.plot(nodes_w/(2.0*scipy.pi),scipy.log(sol),label='dep. x=L; PGD', linewidth=2)
pl.xlabel('freq.')
pl.ylabel('log(dep.)')
pl.grid('on')
pl.legend(loc=1)


pl.show()

pl.show()

