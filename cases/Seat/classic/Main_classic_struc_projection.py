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
results_file='results/classic'

freq_ini     = 10.0
freq_end     = 200.0

nb_freq_step_per_proc=40

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=340.0
rho=1.2

# shell structure
material=[]
material.append(75000.0e6)
material.append(0.33)
material.append(15.0e-3)
material.append(2700.0)

# number of structure modes
nb_mode_S=44
modal_damping_S=0.02

# fluid damping
#fluid_damping=1.0
fluid_damping=(1.0+0.01j)

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
#FixedStrucDofRx=[]
#FixedStrucDofRy=[]
#FixedStrucDofRz=[]

FixedStrucDof=scipy.hstack([FixedStrucDofUx,FixedStrucDofUy,FixedStrucDofUz,FixedStrucDofRx,FixedStrucDofRy,FixedStrucDofRz])
SolvedDofS=scipy.setdiff1d(range(struc_ndof),FixedStrucDof)
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
# Compute eigen modes of the structure
##################################################################
tic = time.clock()

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],nb_mode_S,MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

##################################################################
# Add static solution to the structure basis
##################################################################
tic = time.clock()

Static_mode_S = mumps.spsolve( KSS[SolvedDofS,:][:,SolvedDofS] , FS[SolvedDofS] ,comm=mycomm ).T
#Static_mode_S = scipy.sparse.linalg.spsolve( KSS[SolvedDofS,:][:,SolvedDofS] , FS[SolvedDofS] ).T

toc = time.clock()
if rank==0:
    print ("time for computing the structure static mode:",toc-tic)

##################################################################
# Orthogonalisation of the static mode
##################################################################
tic = time.clock()
eigen_vectors_S=scipy.sparse.csc_matrix(eigen_vectors_S)

lines=list(range(len(SolvedDofS)))
S=scipy.sparse.csc_matrix(Static_mode_S)

MSS_static_S1 = scipy.dot(MSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT__MSS_static_11 = scipy.array(S*MSS_static_S1.todense())[0][0]
S=S/(scipy.sqrt(staticT__MSS_static_11))

for i in range(nb_mode_S):
    a=eigen_vectors_S[lines,:][:,i]
    b=scipy.dot(KSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
    tmp = scipy.array(a.T*b.todense())[0][0]/eigen_values_S[i]
    S=S-tmp*a.T
    MSS_static_S1 = scipy.dot(MSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
    staticT__MSS_static_11 = scipy.array(S*MSS_static_S1.todense())[0][0]
    S=S/(scipy.sqrt(staticT__MSS_static_11))

toc = time.clock()
if rank==0:
    print ("time for orthogonalization of the static mode:",toc-tic)

##################################################################
# Build and save the Structure Basis
##################################################################
tic = time.clock()

PSn = scipy.sparse.construct.bmat( [ [scipy.sparse.coo_matrix(eigen_vectors_S),scipy.sparse.coo_matrix(S).T] ] )

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))
freq_eigv_S.append(0.0)

eigen_vector_S_list=[]
for i in range(PSn.shape[1]):
    Q=scipy.zeros(struc_ndof)
    Q[SolvedDofS]=PSn.todense()[:,i]
    disp=scipy.zeros((struc_nnodes,3))
    disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
    disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
    disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
    eigen_vector_S_list.append(disp)

toc = time.clock()
if rank==0:
    print ("structure eigen frequencies : ",freq_eigv_S)
    print ("time for computing the structure modes:",toc-tic)
    silex_lib_gmsh.WriteResults2(results_file+'_structure_modes',struc_nodes,struc_elements,2,[[eigen_vector_S_list,'nodal',3,'modes']])

##################################################################
# Compute structure damping matrix
##################################################################
VDnn = 2.0*modal_damping_S*scipy.sqrt(eigen_values_S)
IIDnn = list(range(len(SolvedDofF),len(SolvedDofF)+nb_mode_S))
JJDnn = list(range(len(SolvedDofF),len(SolvedDofF)+nb_mode_S))

D = scipy.sparse.coo_matrix( (VDnn,(IIDnn,JJDnn)), shape=(len(SolvedDofF)+nb_mode_S+1,len(SolvedDofF)+nb_mode_S+1) )

##################################################################
# Compute eigen modes of the fluid
##################################################################
tic = time.clock()
nb_mode_F=15
eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(KFF[SolvedDofF,:][:,SolvedDofF],nb_mode_F,MFF[SolvedDofF,:][:,SolvedDofF],sigma=0,which='LM')

freq_eigv_F=list(scipy.sqrt(eigen_values_F)/(2*scipy.pi))
eigen_vector_F_list=[]
for i in range(nb_mode_F):
    tmp=eigen_vectors_F[:,i].real
    eigen_vector_F_list.append(tmp[SolvedDofF])

toc = time.clock()
if rank==0:
    print ("fluid eigen frequencies : ",freq_eigv_F)
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes,fluid_elements,4,[[eigen_vector_F_list,'nodal',1,'pressure']])
    print ("time for computing the fluid modes:",toc-tic)


##################################################################
# Construct the whole system
##################################################################
VK_diag_nn = eigen_values_S
VM_diag_nn = eigen_values_S/eigen_values_S
IIDnn = list(range(nb_mode_S))
JJDnn = list(range(nb_mode_S))

K_diag_nn= scipy.sparse.csc_matrix( (VK_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )
M_diag_nn= scipy.sparse.csc_matrix( (VM_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )

KSS_static_S1 = scipy.dot(KSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT__KSS_static_11 = scipy.array(S*KSS_static_S1.todense())[0][0]

Knn = scipy.sparse.construct.bmat( [[K_diag_nn,None],
                                    [None,staticT__KSS_static_11]
                                    ] )

Mnn = scipy.sparse.construct.bmat( [[M_diag_nn,None],
                                    [None,1.0]
                                    ] )

CnF = scipy.dot(PSn.T,scipy.sparse.coo_matrix(CSF[SolvedDofS,:][:,SolvedDofF]))

K=scipy.sparse.construct.bmat( [[fluid_damping*KFF[SolvedDofF,:][:,SolvedDofF],None],
                                [-CnF,Knn]
                                ] )

M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],CnF.T],
                                [None,Mnn]
                                ] )

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )


FF=scipy.zeros((fluid_ndof))
F=FF[SolvedDofF]
Fn = scipy.dot(PSn.T,scipy.sparse.coo_matrix(FS[SolvedDofS]).T)
F=scipy.append(F,scipy.array(Fn.todense()))
F=scipy.sparse.coo_matrix(F).T

print("hello 9")

##############################################################
# FRF computation of the FSI problem
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
    #for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("proc number",rank,"frequency=",freq)

        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , scipy.array(F.todense() , dtype=complex) )
        sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype='c16') , scipy.array(F.todense() , dtype='c16'), comm=mycomm )
        

        press = scipy.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,press))
        #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
        i=i+1

        if rank==0:
            Q=scipy.zeros((struc_ndof),dtype=float)
            tmp=sol[list(range(len(SolvedDofF),len(SolvedDofF)+PSn.shape[1],1))].real
            Q[SolvedDofS]=PSn*tmp
            disp=scipy.zeros((struc_nnodes,3))
            disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
            disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
            disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
            disp_save.append(disp)
            press_save.append(press.real)

    #frfsave=[scipy.linspace(freq_ini,freq_end,nb_freq_step),frf]
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if rank==0:
        silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])
        silex_lib_gmsh.WriteResults2(results_file+'_results_struct_frf',struc_nodes,struc_elements,2,[[disp_save,'nodal',3,'displacement']])

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())



    # save the FRF problem
    Allfrequencies=scipy.zeros(nb_freq_step)
    Allfrf=scipy.zeros(nb_freq_step)
    k=0
    if rank==0:
        for i in range(nproc):
            data = comm.recv(source=i, tag=11)
            for j in range(len(data[0])):
                Allfrequencies[k]=data[0][j]
                Allfrf[k]=data[1][j]
                k=k+1

        Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
        f=open(results_file+'_results_with_damping.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()

