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

import mumps

import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt
import silex_lib_gmsh

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

if rank==0:
    print ("time at the beginning of the computation:",time.ctime())

# parallepipedic cavity with plane structure
mesh_file='geom/cyl'
results_file='results/xfem_struc_projection_fluid_damping'

# Fluid: air
celerity=340.0
rho=1.2

# Fluid: water
##celerity=1500.0
##rho=1000.0

# shell structure
material=[]
material.append(70000.0e6)
material.append(0.27)
material.append(6.0e-3)
material.append(2700.0)

freq_ini     = 5.0
freq_end     = 1000.0
nb_freq_step_per_proc=80

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=1

# number of structure modes
nb_mode_S=17*5
modal_damping_S=0.02

# fluid damping
#fluid_damping=1.0
fluid_damping=(1.0+0.01j)

##############################################################
# Load fluid mesh
##############################################################
tic = time.clock()
tic0=time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_fluid.msh',3)
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',4,1)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,4)
if rank==0:
    print ("nnodes for fluid=",fluid_nnodes)
    print ("nelem for fluid=",fluid_nelem)

##############################################################
# Load structure mesh
##############################################################
struc_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh',3)
struc_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',2,2)

struc_nnodes   = struc_nodes.shape[0]
struc_nelem    = struc_elements.shape[0]
struc_ndof     = struc_nnodes*6

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_struc_mesh',struc_nodes,struc_elements,2)

if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)

##############################################################
# Material, Boundary conditions
##############################################################

SolvedDofF=list(range(fluid_ndof))

# Find the fixed dofs and the free dofs of the structure

tmp=scipy.sparse.find(struc_nodes[:,0]==0.0)# x=0
FixedStrucNodes=tmp[1]+1
FixedStrucDofUx=(FixedStrucNodes-1)*6
FixedStrucDofUy=(FixedStrucNodes-1)*6+1
FixedStrucDofUz=(FixedStrucNodes-1)*6+2
FixedStrucDofRx=(FixedStrucNodes-1)*6+3
FixedStrucDofRy=(FixedStrucNodes-1)*6+4
FixedStrucDofRz=(FixedStrucNodes-1)*6+5

FixedStrucDof=scipy.hstack([FixedStrucDofUx,FixedStrucDofUy,FixedStrucDofUz,FixedStrucDofRx,FixedStrucDofRy,FixedStrucDofRz])

SolvedDofS=scipy.setdiff1d(range(struc_ndof),FixedStrucDof)

# To impose the load on the structure
tmp=scipy.sparse.find(struc_nodes[:,1]==0.0)# y=0.0
IdNodeLoadStructure=tmp[1]+1

FS=scipy.zeros(struc_ndof)

for i in range(len(IdNodeLoadStructure)-1):
    IddofLoadStructure=[(IdNodeLoadStructure[i]-1)*6]
    FS[IddofLoadStructure]=FS[IddofLoadStructure]+1.0/len(IdNodeLoadStructure)
    IddofLoadStructure=[(IdNodeLoadStructure[i+1]-1)*6]
    FS[IddofLoadStructure]=FS[IddofLoadStructure]+1.0/len(IdNodeLoadStructure)

##################################################################
# compute level set
##################################################################

tic = time.clock()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)

toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_signed_distance',fluid_nodes,fluid_elements,4,[[[LevelSet],'nodal',1,'Level set']])

##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()

LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements,LevelSet)
LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]

EnrichedElements,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements[LSEnrichedElements])
EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.clock()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes,fluid_elements[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],4)

##############################################################
# Compute Standard Fluid Matrices
##############################################################
tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time to compute fluid matrices:",toc-tic)

##############################################################
# Compute structure matrices
##############################################################
tic = time.clock()

IIks,JJks,Vks,Vms=silex_lib_dkt.stiffnessmatrix(struc_nodes,struc_elements,material)

KSS = scipy.sparse.csc_matrix( (Vks,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )
MSS = scipy.sparse.csc_matrix( (Vms,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )

toc = time.clock()
if rank==0:
    print ("time for computing structure:",toc-tic)

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.clock()

Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,celerity,rho)

KAAheaviside = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
MAAheaviside = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
KAFheaviside = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
MAFheaviside = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofA=Enrichednodes-1

KAA=KAAheaviside
MAA=MAAheaviside
KAF=KAFheaviside
MAF=MAFheaviside

toc = time.clock()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

##################################################################
# Compute coupling terms on interface
##################################################################
tic = time.clock()

IIc1,JJc1,Vc1=silex_lib_xfem_acou_tet4.computexfemcoupling1(fluid_nodes,struc_nodes,fluid_elements,struc_elements,EnrichedElements)
IIc2,JJc2,Vc2=silex_lib_xfem_acou_tet4.computexfemcoupling2(fluid_nodes,struc_nodes,fluid_elements,struc_elements,EnrichedElements,LevelSet)

CSA=0.5*scipy.sparse.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(struc_ndof,fluid_ndof) )+0.5*scipy.sparse.csc_matrix( (Vc2,(IIc2,JJc2)), shape=(struc_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time to compute coupling matrices:",toc-tic)

##################################################################
# Compute eigen modes of the structure
##################################################################
tic = time.clock()

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],nb_mode_S,MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

toc = time.clock()
if rank==0:
    print ("time for computing the structure modal basis:",toc-tic)

##################################################################
# Add static solution to the structure basis
##################################################################
tic = time.clock()
Static_mode_S = mumps.spsolve( KSS[SolvedDofS,:][:,SolvedDofS] , FS[SolvedDofS] , comm=mycomm ).T
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
S=scipy.sparse.coo_matrix(Static_mode_S)

MSS_static_S1 = scipy.dot(MSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT__MSS_static_11 = scipy.array(S*MSS_static_S1.todense())[0][0]
S=S/(scipy.sqrt(staticT__MSS_static_11))

for i in range(nb_mode_S):
    a=eigen_vectors_S[lines,:][:,i]
    b=KSS[SolvedDofS,:][:,SolvedDofS]*S.T
    tmp = scipy.array(a.T*b.todense())[0][0]/eigen_values_S[i]
    S=S-tmp*a.T
    MSS_static_S1 = MSS[SolvedDofS,:][:,SolvedDofS]*S.T
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

if rank==0:
    print ("structure eigen frequencies : ",freq_eigv_S)
if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_structure_modes',struc_nodes,struc_elements,2,[[eigen_vector_S_list,'nodal',3,'modes']])


##################################################################
# Compute structure damping matrix
##################################################################
VDnn = 2.0*modal_damping_S*scipy.sqrt(eigen_values_S)
IIDnn = list(range(len(SolvedDofF)+len(SolvedDofA),len(SolvedDofF)+len(SolvedDofA)+nb_mode_S))
JJDnn = list(range(len(SolvedDofF)+len(SolvedDofA),len(SolvedDofF)+len(SolvedDofA)+nb_mode_S))

D = scipy.sparse.coo_matrix( (VDnn,(IIDnn,JJDnn)), shape=(len(SolvedDofF)+len(SolvedDofA)+nb_mode_S+1,len(SolvedDofF)+len(SolvedDofA)+nb_mode_S+1) )

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


CnA = PSn.T*CSA[SolvedDofS,:][:,SolvedDofA]


K=scipy.sparse.construct.bmat( [[fluid_damping*KFF[SolvedDofF,:][:,SolvedDofF],fluid_damping*KAF[SolvedDofF,:][:,SolvedDofA],None],
                                [fluid_damping*KAF[SolvedDofA,:][:,SolvedDofF],fluid_damping*KAA[SolvedDofA,:][:,SolvedDofA],None],
                                [None,-CnA,Knn]
                                ] )
M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA],None],
                                [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA],CnA.T],
                                [None,None,Mnn]
                                ] )


FF = scipy.zeros(fluid_ndof)
FA = scipy.zeros(fluid_ndof)
F  = FF[SolvedDofF]
F  = scipy.append(F,FA[SolvedDofA])
Fn = scipy.dot(PSn.T,scipy.sparse.coo_matrix(FS[SolvedDofS]).T)
F  = scipy.append(F,scipy.array(Fn.todense()))
F  = scipy.sparse.coo_matrix(F).T
##############################################################
# FRF computation of the FSI problem
##############################################################
toc0=time.clock()

if rank==0:
    print ("fixed time before FRF loop:",toc0-tic0)


Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    if rank==0:
        print ("time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq
        if rank==0:
            print ("proc number",rank,"frequency=",freq)

        #sol = scipy.sparse.linalg.spsolve(K-(omega*omega)*M+omega*D*1j, F)
        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype='c16')  , scipy.array(F.todense() , dtype='c16'), comm=mycomm )

        press = scipy.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        enrichment=scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=press
        CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
        #frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press,enrichment,LevelSet,LevelSetTangent))

        if rank==0:
            Q=scipy.zeros((struc_ndof),dtype=float)
            tmp=sol[list(range(len(SolvedDofF)+len(SolvedDofA),len(SolvedDofF)+len(SolvedDofA)+PSn.shape[1],1))].real
            Q[SolvedDofS]=PSn*tmp
            disp=scipy.zeros((struc_nnodes,3))
            disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
            disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
            disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
            disp_save.append(disp)
            press_save.append(CorrectedPressure.real)


    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_struct_frf',struc_nodes,struc_elements,2,[[disp_save,'nodal',3,'displacement']])

    if rank==0:
        print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    # Save the FRF problem
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
        f=open(results_file+'_results.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()
