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
mesh_file='geom/sieges_classic'
results_file='results/classic_pgd'


freq_ini     = 40.0
freq_end     = 60.0

nb_freq_step = 100
#deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)


# air
celerity=340.0
rho=1.2
#fluid_damping=1.0+0.005j
fluid_damping=1.0

# shell structure
material=[]
material.append(75000.0e6)
material.append(0.33)
material.append(15.0e-3)
material.append(2700.0)

if fluid_damping.imag==0.0:
    mytype='float'
else:
    mytype='c16'

print("type de resolution: ",mytype)

val_residu=0.1
val_critere=0.1

##############################################################
# PARAMETER mesh: Omega
##############################################################

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
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,3)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes


##############################################################
# Load structure mesh
##############################################################


struc_elements_old,struc_node_id = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,2)
struc_nnodes   = len(scipy.unique(struc_elements_old))
struc_nelem    = struc_elements_old.shape[0]
struc_ndof     = struc_nnodes*6

struc_boun,struc_boun_id = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',1,1)

# renumbering of structure nodes
#struc_node_id=scipy.unique(struc_elements_old)
struc_nodes=scipy.zeros((struc_nnodes,3))
for i in range(struc_nnodes):
    struc_nodes[i,0]=fluid_nodes[struc_node_id[i]-1,0]
    struc_nodes[i,1]=fluid_nodes[struc_node_id[i]-1,1]
    struc_nodes[i,2]=fluid_nodes[struc_node_id[i]-1,2]

struc_elements=silex_lib_xfem_acou_tet4.changestructureconnectivity(struc_node_id,struc_elements_old)

toc = time.clock()
if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("time for the reading data part:",toc-tic)
    silex_lib_gmsh.WriteResults2(results_file+'_struc_mesh',struc_nodes,struc_elements,2)
    silex_lib_gmsh.WriteResults2(results_file+'_struc_boun_mesh',fluid_nodes,struc_boun,1)

##################################################################
# compute level set: just to make the compatible mesh
##################################################################

tic = time.clock()

LevelSet,LevelSetDist = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)

toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)
    silex_lib_gmsh.WriteResults2(results_file+'_signed_distance',fluid_nodes,fluid_elements,4,[[[LevelSet],'nodal',1,'Level set']])


##################################################################
# Make the compatible fsi mesh at the interface
##################################################################

#print xvibacoufo.makecompatiblefsimesh.__doc__

#struc_boun_id=scipy.unique(struc_boun)
interfaceIdnodes = scipy.setdiff1d(struc_node_id,struc_boun_id)

fluid_elements_new,fluid_nodes_new,interface_elements=silex_lib_xfem_acou_tet4.makecompatiblefsimesh(fluid_nodes,
                                                                            fluid_elements,
                                                                            struc_elements_old,
                                                                            struc_elements,
                                                                            interfaceIdnodes,
                                                                            LevelSet)

fluid_elements = fluid_elements_new
fluid_nodes    = fluid_nodes_new
fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

if rank==0:
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,4)
    silex_lib_gmsh.WriteResults2(results_file+'_interface_mesh_1',fluid_nodes,interface_elements[:,0:3],2)
    silex_lib_gmsh.WriteResults2(results_file+'_interface_mesh_2',struc_nodes,interface_elements[:,3:6],2)
    print ("nnodes for fluid=",fluid_nnodes)
    print ("nelem for fluid=",fluid_nelem)
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)

##############################################################
# Material, Boundary conditions
##############################################################

# Find the fixed dofs and the free dofs

tmp4=scipy.sparse.find(struc_nodes[:,2]==0.0) # z=0
FixedStrucNodes=scipy.unique(scipy.hstack([tmp4[1]+1]))


FixedStrucDofUx=(FixedStrucNodes-1)*6
FixedStrucDofUy=(FixedStrucNodes-1)*6+1
FixedStrucDofUz=(FixedStrucNodes-1)*6+2
FixedStrucDofRx=(FixedStrucNodes-1)*6+3
FixedStrucDofRy=(FixedStrucNodes-1)*6+4
FixedStrucDofRz=(FixedStrucNodes-1)*6+5

##FixedStrucDofUx=list(range(1,struc_nnodes+1,1))*6
##FixedStrucDofUy=list(range(1,struc_nnodes+1,1))*6+1
##FixedStrucDofUz=list(range(1,struc_nnodes+1,1))*6+2
##FixedStrucDofRx=list(range(1,struc_nnodes+1,1))*6+3
##FixedStrucDofRy=list(range(1,struc_nnodes+1,1))*6+4
##FixedStrucDofRz=list(range(1,struc_nnodes+1,1))*6+5

FixedStrucDof=scipy.hstack([FixedStrucDofUx,FixedStrucDofUy,FixedStrucDofUz,FixedStrucDofRx,FixedStrucDofRy,FixedStrucDofRz])
SolvedDofS=scipy.setdiff1d(list(range(struc_ndof)),FixedStrucDof)
#SolvedDofS=struc_ndof

# To impose the load on the structure
IdNodeLoadStructure=13


FS=scipy.zeros(struc_ndof)
#IddofLoadStructure=[(IdNodeLoadStructure-1)*6 , (IdNodeLoadStructure-1)*6+1 , (IdNodeLoadStructure-1)*6+2]
IddofLoadStructure=[(IdNodeLoadStructure-1)*6+2]
FS[IddofLoadStructure]=1.0


##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

#FixedFluidDof=[0]
SolvedDofF=list(range(fluid_ndof))
#SolvedDofF=setdiff1d(range(fluid_ndof),FixedFluidDof)
#pl.spy(KFF,markersize=0.08)
#pl.show()


##############################################################
# Compute structure matrices
##############################################################
tic = time.clock()


IIks,JJks,Vks,Vms=silex_lib_dkt.stiffnessmatrix(struc_nodes,struc_elements,material)

KSS = scipy.sparse.csc_matrix( (Vks,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )
MSS = scipy.sparse.csc_matrix( (Vms,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )

toc = time.clock()
if rank==0:
    print ("time for computing the structure:",toc-tic)

##################################################################
# Compute coupling terms on interface
##################################################################
tic = time.clock()

IIc,JJc,Vc=silex_lib_xfem_acou_tet4.computecoupling(fluid_nodes,interface_elements)

CSF=scipy.sparse.csc_matrix( (Vc,(IIc,JJc)), shape=(struc_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time for computing the coupling:",toc-tic)


##################################################################
# PGD matrices for Omega
##################################################################

IIa,JJa,Va = silex_lib_pgd.frequency_a_matrix(nodes_w,elements_w)
A=scipy.sparse.csc_matrix( (Va,(IIa,JJa)), shape=(omega_ndof,omega_ndof) )
IIb,JJb,Vb = silex_lib_pgd.frequency_b_matrix(nodes_w,elements_w)
B=scipy.sparse.csc_matrix( (Vb,(IIb,JJb)), shape=(omega_ndof,omega_ndof) )

bb=silex_lib_pgd.frequency_om0_second_member(nodes_w,elements_w)

##################################################################
# Construct the whole system
##################################################################


K=scipy.sparse.construct.bmat( [[KFF[SolvedDofF,:][:,SolvedDofF],None],
                                [-CSF[SolvedDofS,:][:,SolvedDofF],KSS[SolvedDofS,:][:,SolvedDofS]]
                                ] )

M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],CSF[SolvedDofS,:][:,SolvedDofF].T],
                                [None,MSS[SolvedDofS,:][:,SolvedDofS]]
                                ] )

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )


PP=scipy.zeros((fluid_ndof))
P=PP[SolvedDofF]
P=scipy.append(P,FS[SolvedDofS])
#P=scipy.sparse.csc_matrix(P)

#############################################################################
#       Solve the problem
#############################################################################

ndof_x=fluid_ndof+struc_ndof
ndof_w=omega_ndof
SolvedDofs_x=scipy.hstack([SolvedDofF,fluid_ndof+SolvedDofS])
SolvedDofs_w=scipy.array(list(range(ndof_w)))
##SolvedDofs_X=scipy.hstack([SolvedDofF,fluid_ndof+SolvedDofS,SolvedDofs_w+fluid_ndof+struc_ndof])
F_i=[]
G_i=[]
alpha_i = scipy.array(  scipy.zeros(5) , dtype=mytype  )

#F=scipy.zeros(ndof_x)+1.0
#G=scipy.zeros(ndof_w)+1.0

sum_fi_Gi=scipy.zeros((ndof_x,ndof_w))
niter=0
#for i in range(nb_fcts_PGD):
residu=99.9
nb_fcts_PGD=0
while (residu>val_residu):
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
    #F=scipy.random.random(ndof_x)
    #G=scipy.random.random(ndof_w)
    F=scipy.array(  scipy.zeros(ndof_x)+1.0 , dtype=mytype  )
    G=scipy.array(  scipy.zeros(ndof_w)+1.0 , dtype=mytype  )
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

    while (critere>val_critere):
        F_new=scipy.zeros(ndof_x, dtype=mytype)
        G_new=scipy.zeros(ndof_w, dtype=mytype)

        G_old=G.copy()
        F_old=F.copy()

        #X_old=X.copy()

        Big_matrix_w = scipy.dot(F[SolvedDofs_x],K*F[SolvedDofs_x])*A*fluid_damping-scipy.dot(F[SolvedDofs_x],M*F[SolvedDofs_x])*B
        second_member_w = bb*scipy.dot(P,F[SolvedDofs_x])-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F[SolvedDofs_x])+scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])
        G_new[SolvedDofs_w] = mumps.spsolve( Big_matrix_w , second_member_w , comm=mycomm).T
        #G_new[SolvedDofs_w] = scipy.sparse.linalg.spsolve( Big_matrix_w , second_member_w )
        G=G_new.copy()

        Big_matrix_x = scipy.dot(G,A*G)*K*fluid_damping-scipy.dot(G,B*G)*M
        second_member_x = scipy.dot(bb,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G)
        F_new[SolvedDofs_x] = mumps.spsolve( Big_matrix_x , second_member_x , comm=mycomm).T
        #F_new[SolvedDofs_x] = scipy.sparse.linalg.spsolve( Big_matrix_x , second_member_x )
        F=F_new.copy()
        #stop


        #Big_matrix=scipy.sparse.construct.bmat( [[scipy.dot(G,A*G)*K-scipy.dot(G,B*G)*M  ,  K_sum_fi_Gi_A-M_sum_fi_Gi_B],
        #        [K_sum_fi_Gi_A.T+M_sum_fi_Gi_B.T,scipy.dot(F[SolvedDofs_x]  ,  K*F[SolvedDofs_x])*A-scipy.dot(F[SolvedDofs_x],M*F[SolvedDofs_x])*B]
        #                                         ] )
        #Second_member = scipy.hstack([  scipy.dot(bb,G)*P  ,  bb*scipy.dot(P,F[SolvedDofs_x])  ])

        #print("hello 2")
        #Big_matrix=scipy.sparse.construct.bmat( [[scipy.dot(G,A*G)*K-scipy.dot(G,B*G)*M  ,  K_sum_fi_Gi_A-M_sum_fi_Gi_B-scipy.tensordot(P,bb,0) ],
        #        [K_sum_fi_Gi_A.T+M_sum_fi_Gi_B.T-scipy.tensordot(bb,P,0) , scipy.dot(F[SolvedDofs_x],K*F[SolvedDofs_x])*A-scipy.dot(F[SolvedDofs_x],M*F[SolvedDofs_x])*B]
        #                                         ] )
        #print("hello 3")
        #Second_member = scipy.hstack([  1e-10*scipy.dot(bb,G)*P  ,  1e-10*bb*scipy.dot(P,F[SolvedDofs_x])  ])
        #print("hello 4")
        
        #X[SolvedDofs_X] = mumps.spsolve( Big_matrix , Second_member , comm=mycomm).T
        #print("hello 5")

        #critere = scipy.dot(X_old-X,X_old-X)
        #critere = scipy.dot(G_old-G,G_old-G)**2 + scipy.dot(F_old-F,F_old-F)**2

        critere = max( scipy.linalg.norm(G_old-G)/scipy.linalg.norm(G_old+G) , scipy.linalg.norm(F_old-F)/scipy.linalg.norm(F_old+F) )
        #critere = max( max( abs( (G_old-G)/(G+G_old)) ) , max( abs( (F_old-F)/(F+F_old) ) ) )

##        Big_matrix_x = scipy.dot(G,A*G)*K-scipy.dot(G,B*G)*M
##        second_member_x = scipy.dot(bb,G)*P-scipy.dot(K_sum_fi_Gi_A,G)+scipy.dot(M_sum_fi_Gi_B,G)
##        Big_matrix_w = scipy.dot(F[SolvedDofs_x],K*F[SolvedDofs_x])*A-scipy.dot(F[SolvedDofs_x],M*F[SolvedDofs_x])*B
##        second_member_w = bb*scipy.dot(P,F[SolvedDofs_x])-scipy.dot(K_sum_fi_Gi_A.T,F[SolvedDofs_x])+scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])
##
##        error_x = scipy.linalg.norm(Big_matrix_x*F[SolvedDofs_x]-second_member_x)/scipy.linalg.norm(second_member_x)
##        error_w = scipy.linalg.norm(Big_matrix_w*G-second_member_w)/scipy.linalg.norm(second_member_w)
##
##        critere = error_x+error_w

        #critere = max( max(abs(G_old-G))/max(abs(G)) , max(abs(F_old-F))/max(abs(F)) )
        niter=niter+1

        #F=F_new
        if niter%15==0:
            print("forcage!")
            critere=0.0
   
        print("critere = ",critere,      "niter = ",niter)
        #c=input('?')
        #if niter==5:
        #    critere=1e-8

        #F[SolvedDofs_x]=X[SolvedDofs_x].copy()
        #G=X[SolvedDofs_w+fluid_ndof+struc_ndof].copy()
    second_member_x = scipy.dot(bb,G)*P-scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping+scipy.dot(M_sum_fi_Gi_B,G)
    second_member_w = bb*scipy.dot(P,F[SolvedDofs_x])-scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F[SolvedDofs_x])+scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])
    residu=max( scipy.linalg.norm(second_member_w)/scipy.linalg.norm(bb*scipy.dot(P,F[SolvedDofs_x]))
                , scipy.linalg.norm(second_member_x)/scipy.linalg.norm(scipy.dot(bb,G)*P) )
##    residu=max( scipy.linalg.norm(second_member_w)/(scipy.linalg.norm(bb*scipy.dot(P,F[SolvedDofs_x]))
##                                                    +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A.T*fluid_damping,F[SolvedDofs_x]))
##                                                    +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B.T,F[SolvedDofs_x])) )
##                , scipy.linalg.norm(second_member_x)/(scipy.linalg.norm(scipy.dot(bb,G)*P)
##                                                       +scipy.linalg.norm(scipy.dot(K_sum_fi_Gi_A,G)*fluid_damping)
##                                                       +scipy.linalg.norm(scipy.dot(M_sum_fi_Gi_B,G)) )
##                )
    if nb_fcts_PGD==5:
        residu=0.0
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
        second_member_alpha[l] = scipy.dot(bb,G_i[l])*scipy.dot(P,F_i[l][SolvedDofs_x])
        for m in range(i+1):
            Ktilde[l,m]=scipy.dot(F_i[l][SolvedDofs_x],K*F_i[m][SolvedDofs_x])*scipy.dot(G_i[l],A*G_i[m])*fluid_damping
            Mtilde[l,m]=scipy.dot(F_i[l][SolvedDofs_x],M*F_i[m][SolvedDofs_x])*scipy.dot(G_i[l],B*G_i[m])

    alpha_i = mumps.spsolve( scipy.sparse.csc_matrix( Ktilde-Mtilde ) , second_member_alpha , comm=mycomm).T
    #print("     Solution alpha_i:",alpha_i)

    nb_fcts_PGD=nb_fcts_PGD+1

##############################################################
# FRF construction
##############################################################
frf=[]
frequencies=[]
press_save=[]
for omi in range(nb_node_w):
    press = scipy.zeros(fluid_ndof)
    freq  = nodes_w[omi]/(2.0*scipy.pi)
    frequencies.append(freq)
    for i in range(nb_fcts_PGD):
        press[SolvedDofF]=press[SolvedDofF]+F_i[i][list(range(len(SolvedDofF)))]*G_i[i][omi]*alpha_i[i]

    press_save.append(press)
    frf.append(silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press))

frfsave=[frequencies,frf]
f=open(results_file+'.frf','wb')
pickle.dump(frfsave, f)
f.close()

#silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])


##############################################################
import pylab as pl
prefsquare=20e-6*20e-6

pl.figure(1)
for i in range(nb_fcts_PGD):
    pl.plot(frequencies,scipy.array(G_i[i]),label='G'+str(i), linewidth=2)

pl.show()

