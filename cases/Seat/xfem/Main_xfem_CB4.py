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
time_init=time.ctime()
tic00=time.clock()
if rank==0:
    print ("time at the beginning of the computation:",time.ctime())

# parallepipedic cavity with plane structure
mesh_file='geom/sieges_xfem'

# number of modes

##nb_mode_F=38
##nb_mode_S=24
##results_file='results/xfem_struc_CB_projection_fluid_damping_38F-24S'

##nb_mode_F=102
##nb_mode_S=32
##results_file='results/xfem_struc_CB_projection_fluid_damping_102F-32S'

nb_mode_F=250
nb_mode_S=44
nb_mode_A=88 # must be > nb_mode_S
results_file='results/xfem_struc_CB4_'+str(nb_mode_F)+'F-'+str(nb_mode_S)+'S-'+str(nb_mode_A)+'A'

# Fluid: air
celerity=340.0
rho=1.2

# shell structure
material=[]
material.append(75000.0e6)
material.append(0.33)
material.append(15.0e-3)
material.append(2700.0)

freq_ini     = 10.0
freq_end     = 200.0
nb_freq_step_per_proc=400

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=0

# structure damping
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
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',4,3)

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
struc_boun,tmp     = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',1,1)

struc_nnodes   = struc_nodes.shape[0]
struc_nelem    = struc_elements.shape[0]
struc_ndof     = struc_nnodes*6

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_struc_mesh',struc_nodes,struc_elements,2)
    silex_lib_gmsh.WriteResults2(results_file+'_struc_boun_mesh',struc_nodes,struc_boun,1)

if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)


##############################################################
# Material, Boundary conditions
##############################################################
SolvedDofF=list(range(fluid_ndof))

# Find the fixed dofs and the free dofs of the structure

tmp=scipy.sparse.find(struc_nodes[:,2]==0.0)# z=0
FixedStrucNodes=tmp[1]+1
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

# To impose the load on the structure
IdNodeLoadStructure=13

FS=scipy.zeros(struc_ndof)
#IddofLoadStructure=[(IdNodeLoadStructure-1)*6 , (IdNodeLoadStructure-1)*6+1 , (IdNodeLoadStructure-1)*6+2]
IddofLoadStructure=[(IdNodeLoadStructure-1)*6+2]
FS[IddofLoadStructure]=1.0

##################################################################
# compute level set
##################################################################
tic = time.clock()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)


toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)


tic = time.clock()

tangent_nodes,tangent_mesh=silex_lib_xfem_acou_tet4.buildtangentedgemesh(struc_nodes,struc_elements,struc_boun)

LevelSetTangent,tmp = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,tangent_nodes,tangent_mesh)

toc = time.clock()
if rank==0:
    print ("time to compute tangent level set:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_signed_distance',fluid_nodes,fluid_elements,4,[[[LevelSet],'nodal',1,'Level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,4,[[[LevelSetTangent],'nodal',1,'Tangent level set']])
    #silex_lib_gmsh.WriteResults2(results_file+'_distance',fluid_nodes,fluid_elements,4,[[[LevelSetDist],'nodal',1,'Distance']])
    silex_lib_gmsh.WriteResults2(results_file+'_tangent_mesh',tangent_nodes,tangent_mesh,2)



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


tic = time.clock()

EdgeEnrichedElements,nbenrelts = silex_lib_xfem_acou_tet4.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes,fluid_elements)
EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1


EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements,LevelSetTangent)
EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])


toc = time.clock()
if rank==0:
    print ("time to find edge enriched elements:",toc-tic)

HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

AllElementsExceptEdgeEnrichedElements=scipy.setdiff1d(range(fluid_nelem),EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes,fluid_elements[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements',fluid_nodes,fluid_elements[EdgeEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements_in_all_mesh',fluid_nodes,fluid_elements[EdgeEnrichedElementsInAllMesh],4)
    silex_lib_gmsh.WriteResults2(results_file+'_heaviside_enriched_elements',fluid_nodes,fluid_elements[HeavisideEnrichedElements],4)

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

#Enrichednodes = scipy.unique(fluid_elements[HeavisideEnrichedElements])
Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

NegativeLSelements,PositiveLSelements,NegativeLStgtElements,PositiveLStgtElements,nbNegLS,nbPosLS,nbNegLSt,nbPosLSt=silex_lib_xfem_acou_tet4.getpositivenegativeelts(fluid_elements,LevelSet,LevelSetTangent)

NegativeLSelements=NegativeLSelements[list(range(nbNegLS))]
PositiveLSelements=PositiveLSelements[list(range(nbPosLS))]
NegativeLStgtElements=NegativeLStgtElements[list(range(nbNegLSt))]
PositiveLStgtElements=PositiveLStgtElements[list(range(nbPosLSt))]


#IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements[AllElementsExceptEdge],fluid_nodes,LevelSet,celerity,rho)
IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements[NegativeLStgtElements],fluid_nodes,LevelSet,celerity,rho)
#IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,celerity,rho)

KAAheaviside = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
MAAheaviside = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
KAFheaviside = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
MAFheaviside = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofA=Enrichednodes-1

toc = time.clock()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

##################################################################
# Compute Edge enrichment
##################################################################
tic = time.clock()

PartiallyPositiveLStgtElements=scipy.hstack([PositiveLStgtElements,EdgeEnrichedElementsInAllMesh])

#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes,fluid_elements[EdgeEnrichedElements],LevelSet,LevelSetTangent,celerity,rho)
II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes,fluid_elements[PartiallyPositiveLStgtElements],LevelSet,LevelSetTangent,celerity,rho)
#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes,fluid_elements,LevelSet,LevelSetTangent,celerity,rho)

KAAedge = scipy.sparse.csc_matrix( (vkaa,(II,JJ)), shape=(fluid_ndof,fluid_ndof) )
MAAedge = scipy.sparse.csc_matrix( (vmaa,(II,JJ)), shape=(fluid_ndof,fluid_ndof) )
KAFedge = scipy.sparse.csc_matrix( (vkfa,(II,JJ)), shape=(fluid_ndof,fluid_ndof) )
MAFedge = scipy.sparse.csc_matrix( (vmfa,(II,JJ)), shape=(fluid_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time to compute edge enrichment:",toc-tic)

KAA=KAAheaviside+KAAedge
MAA=MAAheaviside+MAAedge
KAF=KAFheaviside+KAFedge
MAF=MAFheaviside+MAFedge

#KAA=KAAheaviside
#MAA=MAAheaviside
#KAF=KAFheaviside
#MAF=MAFheaviside

#KAA=KAAedge
#MAA=MAAedge
#KAF=KAFedge
#MAF=MAFedge

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_NegativeLSelements',fluid_nodes,fluid_elements[NegativeLSelements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PositiveLSelements',fluid_nodes,fluid_elements[PositiveLSelements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_NegativeLStgtElements',fluid_nodes,fluid_elements[NegativeLStgtElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PositiveLStgtElements',fluid_nodes,fluid_elements[PositiveLStgtElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PartiallyPositiveLStgtElements',fluid_nodes,fluid_elements[PartiallyPositiveLStgtElements],4)

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
    print ("time for computing the structure modal basis (fixed structure):",toc-tic)

tic = time.clock()

eigen_values_free_S,eigen_vectors_free_S= scipy.sparse.linalg.eigsh(KSS,nb_mode_A-nb_mode_S-1,MSS,sigma=0,which='LM')
toc = time.clock()
if rank==0:
    print ("time for computing the structure modal basis (free structure):",toc-tic)

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

PSn = scipy.sparse.construct.bmat( [ [scipy.sparse.csc_matrix(eigen_vectors_S),scipy.sparse.csc_matrix(S).T] ] )

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
IIDnn = list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+nb_mode_S))
JJDnn = list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+nb_mode_S))

D = scipy.sparse.coo_matrix( (VDnn,(IIDnn,JJDnn)), shape=(nb_mode_F+nb_mode_A+nb_mode_S+1,nb_mode_F+nb_mode_A+nb_mode_S+1) )

##################################################################
# Compute eigen modes of the fluid
##################################################################
tic = time.clock()

eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(KFF[SolvedDofF,:][:,SolvedDofF],nb_mode_F,MFF[SolvedDofF,:][:,SolvedDofF],sigma=0,which='LM')

freq_eigv_F=list(scipy.sqrt(eigen_values_F)/(2*scipy.pi))

if rank==0:
    print ("fluid eigen frequencies : ",freq_eigv_F)
 
if (flag_write_gmsh_results==1) and (rank==0):
    eigen_vector_F_list=[]
    for i in range(nb_mode_F):
        tmp=eigen_vectors_F[:,i].real
        eigen_vector_F_list.append(tmp[SolvedDofF])
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes,fluid_elements,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

toc = time.clock()
if rank==0:
    print ("time for computing the fluid modes:",toc-tic)


##################################################################
# Compute Psi_Fk and Psi_Ak for the fluid
##################################################################
tic = time.clock()

Psi_Fk=scipy.zeros((len(SolvedDofF),nb_mode_A))
Psi_Ak=scipy.zeros((len(SolvedDofA),nb_mode_A))
ModeS=scipy.zeros((struc_ndof,nb_mode_A))

K=scipy.sparse.construct.bmat( [[KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
                                [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]
                                ] )

M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                ] )

for i in range(nb_mode_S+1):
    ModeS[SolvedDofS,i]=scipy.array(PSn.todense()[:,i].T)[0]

for i in range(nb_mode_A-nb_mode_S-1):
    ModeS[:,i+nb_mode_S+1]=eigen_vectors_free_S[:,i]

OmegS=scipy.hstack([eigen_values_S,0.0,eigen_values_free_S])

toc = time.clock()
if rank==0:
    print ("time to form the system to solve PSI_Fk and Psi_Ak:",toc-tic)

omega_cst=200.0*2.0*scipy.pi

tic = time.clock()
MySolve = scipy.sparse.linalg.factorized( scipy.sparse.csc_matrix(K-omega_cst*omega_cst*M) ) # Makes LU decomposition.
toc = time.clock()
if rank==0:
    print ("time to factorized the system to solve PSI_Fk and Psi_Ak:",toc-tic)

tic = time.clock()
for k in range(nb_mode_A):
    tmp=scipy.hstack( [scipy.zeros(len(SolvedDofF)),CSA[:,SolvedDofA].T*ModeS[:,k]] )
    Xi=MySolve( tmp)
    #Xi = mumps.spsolve( scipy.sparse.csc_matrix(K-OmegS[k]*OmegS[k]*M) , tmp , comm=mycomm )
    #Xi = mumps.spsolve( scipy.sparse.csc_matrix(K-omega_cst*omega_cst*M) , tmp , comm=mycomm )
    #if rank==0:
    #    print("k=",k)

    Psi_Fk[:,k]=Xi[list(range(len(SolvedDofF)))]
    Psi_Ak[:,k]=Xi[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]


#Psi_Fk=scipy.sparse.csc_matrix(Psi_Fk)
#Psi_Ak=scipy.sparse.csc_matrix(Psi_Ak)

toc = time.clock()
if rank==0:
    print ("time to compute PSI_Fk and Psi_Ak:",toc-tic)

#if (flag_write_gmsh_results==1) and (rank==0):
#    silex_lib_gmsh.WriteResults2(results_file+'_Psi_Fk',fluid_nodes,fluid_elements,4,[[[Psi_FA[:,40].todense(),Psi_FA[:,50].todense(),Psi_FA[:,80].todense(),Psi_FA[:,300].todense()],'nodal',1,'pressure']])

##################################################################
# Construct the whole system
##################################################################

# Fluid part
tic = time.clock()

VK_diag_pp = eigen_values_F
VM_diag_pp = eigen_values_F/eigen_values_F
IIDpp = list(range(nb_mode_F))
JJDpp = list(range(nb_mode_F))

K_diag_pp= scipy.sparse.csc_matrix( (VK_diag_pp,(IIDpp,JJDpp)), shape=(nb_mode_F,nb_mode_F) )
M_diag_pp= scipy.sparse.csc_matrix( (VM_diag_pp,(IIDpp,JJDpp)), shape=(nb_mode_F,nb_mode_F) )


K_pk_1 = scipy.dot(scipy.array(eigen_vectors_F.T*KFF[SolvedDofF,:][:,SolvedDofF]),Psi_Fk)
K_pk_2 = scipy.dot(scipy.array(eigen_vectors_F.T*KAF[SolvedDofF,:][:,SolvedDofA]),Psi_Ak)
K_pk = K_pk_1+K_pk_2

K_kk_1 = scipy.dot(Psi_Fk.T*KFF[SolvedDofF,:][:,SolvedDofF],Psi_Fk)
K_kk_2 = scipy.dot(Psi_Fk.T*KAF[SolvedDofF,:][:,SolvedDofA],Psi_Ak)
K_kk_3 = K_kk_2.T
K_kk_4 = scipy.dot(Psi_Ak.T*KAA[SolvedDofA,:][:,SolvedDofA],Psi_Ak)

K_kk = K_kk_1+K_kk_2+K_kk_3+K_kk_4

CnA = PSn.T*CSA[SolvedDofS,:][:,SolvedDofA]
Cnk = CnA*Psi_Ak

toc = time.clock()
if rank==0:
    print ("time to compute K projections:",toc-tic)

tic = time.clock()

M_pk_1 = scipy.dot(scipy.array(eigen_vectors_F.T*MFF[SolvedDofF,:][:,SolvedDofF]),Psi_Fk)
M_pk_2 = scipy.dot(scipy.array(eigen_vectors_F.T*MAF[SolvedDofF,:][:,SolvedDofA]),Psi_Ak)
M_pk = M_pk_1+M_pk_2

M_kk_1 = scipy.dot(Psi_Fk.T*MFF[SolvedDofF,:][:,SolvedDofF],Psi_Fk)
M_kk_2 = scipy.dot(Psi_Fk.T*MAF[SolvedDofF,:][:,SolvedDofA],Psi_Ak)
M_kk_3 = M_kk_2.T
M_kk_4 = scipy.dot(Psi_Ak.T*MAA[SolvedDofA,:][:,SolvedDofA],Psi_Ak)

M_kk = M_kk_1+M_kk_2+M_kk_3+M_kk_4

toc = time.clock()
if rank==0:
    print ("time to compute M projections:",toc-tic)

tic = time.clock()

eigen_vectors_F=scipy.sparse.csc_matrix(eigen_vectors_F)

VK_diag_nn = eigen_values_S
VM_diag_nn = eigen_values_S/eigen_values_S
IIDnn = list(range(nb_mode_S))
JJDnn = list(range(nb_mode_S))

K_diag_nn= scipy.sparse.csc_matrix( (VK_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )
M_diag_nn= scipy.sparse.csc_matrix( (VM_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )

KSS_static_S1 = scipy.dot(KSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT_KSS_static_11 = scipy.array(S*KSS_static_S1.todense())[0][0]

Knn = scipy.sparse.construct.bmat( [[K_diag_nn,None],[None,staticT_KSS_static_11]] )

Mnn = scipy.sparse.construct.bmat( [[M_diag_nn,None],[None,1.0]] )



K=scipy.sparse.construct.bmat( [ [fluid_damping*K_diag_pp,K_pk,None],[K_pk.T,fluid_damping*K_kk,None],[None,-Cnk,Knn] ] )

M=scipy.sparse.construct.bmat( [ [M_diag_pp,M_pk,None],[M_pk.T,M_kk,Cnk.T],[None,None,Mnn] ] )

F  = scipy.zeros((nb_mode_F+nb_mode_A))
Fn = PSn.T*FS[SolvedDofS]
F  = scipy.append(F,Fn)
F  = scipy.sparse.csc_matrix(F).T

##############################################################

if rank==0:
    print ("nb. fluid modes = ",nb_mode_F)
    print ("nb. enriched nodes = ",len(SolvedDofA))
    print ("nb. structure modes = ",nb_mode_S)

##############################################################
# FRF computation of the FSI problem
##############################################################
toc0=time.clock()
time_before_frf=time.ctime()
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
        print ("proc number",rank,"frequency=",freq)

        #sol = scipy.sparse.linalg.spsolve(K-(omega*omega)*M+omega*D*1j, F)
        #sol = scipy.linalg.solve(scipy.array((K-(omega*omega)*M+omega*D*1j).todense()),scipy.array(F.todense()))
        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype='c16')  , scipy.array(F.todense() , dtype='c16'), comm=mycomm )

        alpha_p    = scipy.sparse.csc_matrix(sol[list(range(nb_mode_F))])
        alpha_k    = scipy.sparse.csc_matrix(sol[list(range(nb_mode_F,nb_mode_F+nb_mode_A,1))])
        P_A        = Psi_Ak*alpha_k
        press      = eigen_vectors_F*alpha_p+Psi_Fk*alpha_k
        enrichment = scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=P_A
        CorrectedPressure=press
        CorrectedPressure[SolvedDofA]=(CorrectedPressure[SolvedDofA].T+scipy.array(enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)).T
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
        #frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press.todense(),enrichment,LevelSet,LevelSetTangent))

        if rank==0:
            Q=scipy.zeros((struc_ndof),dtype=float)
            tmp=sol[list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+PSn.shape[1],1))].real
            Q[SolvedDofS]=PSn*tmp
            disp=scipy.zeros((struc_nnodes,3))
            disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
            disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
            disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
            disp_save.append(disp)
            press_save.append(CorrectedPressure.real)

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_struct_frf',struc_nodes,struc_elements,2,[[disp_save,'nodal',3,'displacement']])
    
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
        print("Real time at the beginning = ",time_init)
        print("Real time before FRF       = ",time_before_frf)
        print("Real time at the end       = ",time.ctime())
        print("Total time = ",time.clock()-tic00)

