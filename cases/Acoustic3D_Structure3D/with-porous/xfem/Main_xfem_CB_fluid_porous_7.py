###########################################################
## AIR CAVITY
## POROUS MATERIAL
## THIN RIGID STRUCTURE
## XFEM
## CB REDUCTION FOR AIR
###########################################################

## Wrong results because nothing is done for the air-porous interface

###########################################################
# Libraries
###########################################################

import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pickle

import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_gmsh

import silex_lib_porous_tet4_fortran

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

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

if rank==0:
    print ("time at the beginning of the computation:",time.ctime())

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

##############################################################
# Datas
##############################################################

mesh_file='geom/cavity7_with_porous_air'
results_file='results/cavity7_CB_with_porous_air'

flag_write_gmsh_results=1

nb_mode_F = 50

freq_ini     = 10.0
freq_end     = 120.0
nb_freq_step_per_proc=10

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=343.0 # ok
rho=1.21 # ok

# porous properties
E_sol = 286000.0*(3*1144.0+2*286.0)/(1144.0+286.0) #ok =800800.0
nu_sol = 1144.0/(2*(1144.0+286.0)) #ok =0.4
ro_sol = 30.0

to_por = 7.8 # ok
po_por = 0.9 # ok
sg_por = 25000.0 # ok
lambda_por = 226e-6 # ok
lambda_prime_por = 226e-6 # ok

ro_fl = 1.21 # ok
visco_fl = 0.0000184 # ok
pdtl_fl = 0.71 # ok
gamma_fl = 1.4 # ok
p0_fl = 101325.0
ce_fl = celerity 

##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
fluid_elements_S3,IdNodesS3 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,3)# air-porous interface surface
fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem1   = fluid_elements1.shape[0]
fluid_nelem2   = fluid_elements2.shape[0]
fluid_nelem3   = fluid_elements_S3.shape[0]
fluid_nelem4   = fluid_elements_S4.shape[0]
fluid_nelem5   = fluid_elements5.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
fluid_nnodes2 = IdNodes2.shape[0]
fluid_nnodes3 = IdNodesS3.shape[0]
fluid_nnodes4 = IdNodesS4.shape[0]
fluid_nnodes5 = IdNodes5.shape[0]

fluid_ndof1     = fluid_nnodes1
fluid_ndof2     = fluid_nnodes2*6

if rank==0:
    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)
    print ("Number of elements in porous:",fluid_nelem2)

    print ("Number of nodes in air:",fluid_nnodes1)
    print ("Number of nodes in porous:",fluid_nnodes2)
    print ("Number of nodes at interface:",fluid_nnodes3)

# renumbering air
old = scipy.unique(fluid_elements1)
new = list(range(1,len(scipy.unique(fluid_elements1))+1))
new_nodes=fluid_nodes[scipy.unique(fluid_elements1)-1,:]

dico1 = dict(zip(old,new))
new_elements=scipy.zeros((fluid_nelem1,4),dtype=int)
for e in range(fluid_nelem1):
    for i in range(4):
        new_elements[e,i]=dico1[fluid_elements1[e][i]]

fluid_elements1 = new_elements
fluid_nodes1    = new_nodes

new_elements=scipy.zeros((fluid_nelem5,4),dtype=int)
for e in range(fluid_nelem5):
    for i in range(4):
        new_elements[e,i]=dico1[fluid_elements5[e][i]]

fluid_elements5 = new_elements

# renumbering porous
old = scipy.unique(fluid_elements2)
new = list(range(1,len(scipy.unique(fluid_elements2))+1))
new_nodes=fluid_nodes[scipy.unique(fluid_elements2)-1,:]

dico2 = dict(zip(old,new))
new_elements=scipy.zeros((fluid_nelem2,4),dtype=int)
for e in range(fluid_nelem2):
    for i in range(4):
        new_elements[e,i]=dico2[fluid_elements2[e][i]]

fluid_elements2 = new_elements
fluid_nodes2    = new_nodes

# renumbering porous external surfaces
for i in range(len(IdNodesS4)):
    IdNodesS4[i]=dico2[IdNodesS4[i]]

# Boundary conditions on air cavity
IdNodesFixed_porous_us_x=IdNodesS4
##IdNodesFixed_porous_us_y=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS6]))
IdNodesFixed_porous_us_y=IdNodesS4
IdNodesFixed_porous_us_z=IdNodesS4
IdNodesFixed_porous_uf_x=IdNodesS4
IdNodesFixed_porous_uf_y=IdNodesS4
IdNodesFixed_porous_uf_z=IdNodesS4

Fixed_Dofs_porous = scipy.hstack([(IdNodesFixed_porous_us_x-1)*6,
                                  (IdNodesFixed_porous_us_y-1)*6+1,
                                  (IdNodesFixed_porous_us_z-1)*6+2,
                                  (IdNodesFixed_porous_uf_x-1)*6+3,
                                  (IdNodesFixed_porous_uf_y-1)*6+4,
                                  (IdNodesFixed_porous_uf_z-1)*6+5
                                  ])

# get connectivity at air-porous interface
IdNodesS3_for_1=scipy.zeros(fluid_nnodes3,dtype=int)
IdNodesS3_for_2=scipy.zeros(fluid_nnodes3,dtype=int)
for i in range(fluid_nnodes3):
    IdNodesS3_for_1[i]=dico1[IdNodesS3[i]]
    IdNodesS3_for_2[i]=dico2[IdNodesS3[i]]

InterfaceConnectivity=scipy.zeros((fluid_nelem3,6),dtype=int)
for e in range(fluid_nelem3):
    for i in range(3):
        InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]]
        InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]]

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,4)
    silex_lib_gmsh.WriteResults(results_file+'_porous_material_Mesh2',fluid_nodes2,fluid_elements2,4)
    silex_lib_gmsh.WriteResults(results_file+'_porous_air_interface_Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
    silex_lib_gmsh.WriteResults(results_file+'_porous_fixed_Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
    silex_lib_gmsh.WriteResults(results_file+'_air_controlled_volume_Mesh5',fluid_nodes1,fluid_elements5,4)

##############################################################
# Load structure mesh
##############################################################
struc_nodes        = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh',3)
struc_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',2,6)
struc_boun,tmp     = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',1,7)

struc_nnodes   = struc_nodes.shape[0]
struc_nelem    = struc_elements.shape[0]
struc_ndof     = struc_nnodes*6

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_struc_mesh_surface_6',struc_nodes,struc_elements,2)
    silex_lib_gmsh.WriteResults2(results_file+'_struc_boun_mesh_line_7',struc_nodes,struc_boun,1)

if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)

##################################################################
# compute level set
##################################################################

tic = time.clock()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes1,struc_nodes,struc_elements)

toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)

tic = time.clock()

tangent_nodes,tangent_mesh=silex_lib_xfem_acou_tet4.buildtangentedgemesh(struc_nodes,struc_elements,struc_boun)

LevelSetTangent,tmp = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes1,tangent_nodes,tangent_mesh)

toc = time.clock()
if rank==0:
    print ("time to compute tangent level set:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LS_signed_distance',fluid_nodes1,fluid_elements1,4,[[[LevelSet],'nodal',1,'Level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_tangent_level_set',fluid_nodes1,fluid_elements1,4,[[[LevelSetTangent],'nodal',1,'Tangent level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_distance',fluid_nodes1,fluid_elements1,4,[[[distance],'nodal',1,'Distance']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_tangent_mesh',tangent_nodes,tangent_mesh,2)

##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()

LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSet)
LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]
EnrichedElements,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes1,fluid_elements1[LSEnrichedElements])
EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.clock()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

tic = time.clock()

EdgeEnrichedElements,nbenrelts = silex_lib_xfem_acou_tet4.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes1,fluid_elements1)
EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1

EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSetTangent)
EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])

toc = time.clock()
if rank==0:
    print ("time to find edge enriched elements:",toc-tic)

HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

AllElementsExceptEdgeEnrichedElements=scipy.setdiff1d(range(fluid_nelem1),EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes1,fluid_elements1[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes1,fluid_elements1[EnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements',fluid_nodes1,fluid_elements1[EdgeEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements_in_all_mesh',fluid_nodes1,fluid_elements1[EdgeEnrichedElementsInAllMesh],4)
    silex_lib_gmsh.WriteResults2(results_file+'_heaviside_enriched_elements',fluid_nodes1,fluid_elements1[HeavisideEnrichedElements],4)


##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

##############################################################
# Compute Porous Matrices
##############################################################
porous_material_prop=[E_sol,nu_sol,ro_sol,to_por,po_por,sg_por,lambda_por,lambda_prime_por,ro_fl,visco_fl,pdtl_fl,gamma_fl,p0_fl,ce_fl]

SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),Fixed_Dofs_porous)

##############################################################
# Compute Coupling Porous-air Matrices
##############################################################

#print(silex_lib_porous_tet4_fortran.computecouplingporousair.__doc__)

IIpf,JJpf,Vpf=silex_lib_porous_tet4_fortran.computecouplingporousair(fluid_nodes1,InterfaceConnectivity,po_por)
CPF=scipy.sparse.csc_matrix( (Vpf,(IIpf,JJpf)), shape=(fluid_ndof2,fluid_ndof1) )

#SolvedDof = scipy.hstack([SolvedDofF,SolvedDofP+fluid_ndof1])

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.clock()

#Enrichednodes = scipy.unique(fluid_elements1[HeavisideEnrichedElements])
Enrichednodes = scipy.unique(fluid_elements1[EnrichedElements])

NegativeLSelements,PositiveLSelements,NegativeLStgtElements,PositiveLStgtElements,nbNegLS,nbPosLS,nbNegLSt,nbPosLSt=silex_lib_xfem_acou_tet4.getpositivenegativeelts(fluid_elements1,LevelSet,LevelSetTangent)

NegativeLSelements=NegativeLSelements[list(range(nbNegLS))]
PositiveLSelements=PositiveLSelements[list(range(nbPosLS))]
NegativeLStgtElements=NegativeLStgtElements[list(range(nbNegLSt))]
PositiveLStgtElements=PositiveLStgtElements[list(range(nbPosLSt))]


#IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements1[AllElementsExceptEdge],fluid_nodes1,LevelSet,celerity,rho)
IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements1[NegativeLStgtElements],fluid_nodes1,LevelSet,celerity,rho)
#IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements1,fluid_nodes1,LevelSet,celerity,rho)

KAAheaviside = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof1,fluid_ndof1) )
MAAheaviside = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof1,fluid_ndof1) )
KAFheaviside = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof1,fluid_ndof1) )
MAFheaviside = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofA=Enrichednodes-1

toc = time.clock()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

##################################################################
# Compute Edge enrichment
##################################################################
tic = time.clock()

PartiallyPositiveLStgtElements=scipy.hstack([PositiveLStgtElements,EdgeEnrichedElementsInAllMesh])

#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],LevelSet,LevelSetTangent,celerity,rho)
II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[PartiallyPositiveLStgtElements],LevelSet,LevelSetTangent,celerity,rho)
#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1,LevelSet,LevelSetTangent,celerity,rho)

KAAedge = scipy.sparse.csc_matrix( (vkaa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAAedge = scipy.sparse.csc_matrix( (vmaa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
KAFedge = scipy.sparse.csc_matrix( (vkfa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAFedge = scipy.sparse.csc_matrix( (vmfa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )

toc = time.clock()
if rank==0:
    print ("time to compute edge enrichment:",toc-tic)

KAA=KAAheaviside+KAAedge
MAA=MAAheaviside+MAAedge
KAF=KAFheaviside+KAFedge
MAF=MAFheaviside+MAFedge

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_NegativeLSelements',fluid_nodes1,fluid_elements1[NegativeLSelements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PositiveLSelements',fluid_nodes1,fluid_elements1[PositiveLSelements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_NegativeLStgtElements',fluid_nodes1,fluid_elements1[NegativeLStgtElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PositiveLStgtElements',fluid_nodes1,fluid_elements1[PositiveLStgtElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_PartiallyPositiveLStgtElements',fluid_nodes1,fluid_elements1[PartiallyPositiveLStgtElements],4)

##################################################################
# Construct the whole system
##################################################################

# To impose the load on the fluid:
# fluid node number 1
UF = scipy.zeros(fluid_ndof1,dtype=float)
UF[9-1]=3.1250E-05

SolvedDof = scipy.hstack([SolvedDofF,SolvedDofA+fluid_ndof1,SolvedDofP+2*fluid_ndof1])

##################################################################
# Compute eigen modes of the fluid
##################################################################
tic = time.clock()

eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(KFF[SolvedDofF,:][:,SolvedDofF],nb_mode_F,MFF[SolvedDofF,:][:,SolvedDofF],sigma=0,which='LM')

freq_eigv_F=list(scipy.sqrt(eigen_values_F)/(2*scipy.pi))
eigen_vector_F_list=[]
for i in range(nb_mode_F):
    tmp=eigen_vectors_F[:,i].real
    eigen_vector_F_list.append(tmp[SolvedDofF])

if rank==0:
    print ("fluid eigen frequencies : ",freq_eigv_F)
 
if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes1,fluid_elements1,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

toc = time.clock()
if rank==0:
    print ("time for computing the fluid modes:",toc-tic)

##################################################################
# Compute Psi_FA for the fluid: Psi_FA = - KFF^{-1} * KFA
##################################################################
tic = time.clock()

if rank==0:
    print ("Compute PSI_FA")

omega_cst=1.0*2.0*scipy.pi
MySolve = scipy.sparse.linalg.factorized( KFF-omega_cst*omega_cst*MFF ) # Makes LU decomposition.

if rank==0:
    print("LU decomposition has been made")

if rank==0:
    Psi_FA=scipy.zeros((len(SolvedDofF),len(SolvedDofA)))

i=0
j=0
while j<len(SolvedDofA):
    j=i+rank
    if j<len(SolvedDofA):
        One_dof=[SolvedDofA[j]]
        KFA_i_column=-KAF[SolvedDofF,:][:,One_dof]
        tmp=scipy.zeros(len(SolvedDofF))
        tmp[SolvedDofF]=KFA_i_column.todense()
        Xi=MySolve( tmp )
        if rank!=0:
            comm.send([Xi,j], dest=0, tag=11)
        if rank==0:
            for k in range(nproc):
                if k!=0:
                    if i+k<len(SolvedDofA):
                        [Xi,j]=comm.recv(source=k, tag=11)
                        Psi_FA[SolvedDofF,j]=Xi
                else:
                    Psi_FA[SolvedDofF,j]=Xi
        i=i+nproc

if rank==0:
    print("End of PSI_FA computing, send to other proc.")

if rank==0:
    for i in range(nproc):
        if i!=0:
            comm.send(Psi_FA, dest=i, tag=11)
if rank!=0:
    Psi_FA=comm.recv(source=0, tag=11)

Psi_FA=scipy.sparse.csc_matrix(Psi_FA)
toc = time.clock()

if rank==0:
    print ("time to compute PSI_FA:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_Psi_FA',fluid_nodes1,fluid_elements1,4,[[[Psi_FA[:,40].todense(),Psi_FA[:,50].todense(),Psi_FA[:,80].todense()],'nodal',1,'PSI FA pressure']])


##################################################################
# Construct the whole system
##################################################################

# Fluid part
tic = time.clock()

VK_diag_mm = eigen_values_F
VM_diag_mm = eigen_values_F/eigen_values_F
IIDmm = list(range(nb_mode_F))
JJDmm = list(range(nb_mode_F))

PhiFm=eigen_vectors_F

K_diag_mm= scipy.sparse.csc_matrix( (VK_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
M_diag_mm= scipy.sparse.csc_matrix( (VM_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
Khat_AA = KAA[SolvedDofA,:][:,SolvedDofA]+Psi_FA.T*KAF[SolvedDofF,:][:,SolvedDofA]
toc = time.clock()
if rank==0:
    print ("time to compute Khat_hat:",toc-tic)

tic = time.clock()

Mstar_FA = MFF[SolvedDofF,:][:,SolvedDofF]*Psi_FA+MAF[SolvedDofF,:][:,SolvedDofA]
toc = time.clock()
if rank==0:
    print ("time to compute Mstar:",toc-tic)

tic = time.clock()

Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.sparse.csc_matrix((Psi_FA.T).todense()*Mstar_FA.todense())+MAF[SolvedDofA,:][:,SolvedDofF]*Psi_FA

toc = time.clock()
if rank==0:
    print ("time to compute Mhat:",toc-tic)

CmP=PhiFm.T*CPF[SolvedDofP,:][:,SolvedDofF].T
CAP=Psi_FA.T*CPF[SolvedDofP,:][:,SolvedDofF].T

eigen_vectors_F=scipy.sparse.csc_matrix(eigen_vectors_F)
Mhat_mA = eigen_vectors_F.T*Mstar_FA


Freduced_F=scipy.array(scipy.dot(PhiFm.T,UF[SolvedDofF]) , dtype='c16')

##############################################################
# FRF computation
##############################################################

Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    press_save=[]

    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("proc number",rank,"frequency=",freq)

        tic = time.clock()
        IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)
        KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
        MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )

        K=scipy.sparse.construct.bmat( [ [K_diag_mm,None,-CmP],
                                         [None,Khat_AA,-CAP],
                                         [None,None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )

        M=scipy.sparse.construct.bmat( [ [M_diag_mm,Mhat_mA,None],
                                         [Mhat_mA.T,Mhat_AA,None],
                                         [CmP.T,CAP.T,MPP[SolvedDofP,:][:,SolvedDofP]] ] )

##        K=scipy.sparse.construct.bmat( [ [K_diag_mm,None,-CmP],
##                                         [None,Khat_AA,None],
##                                         [None,None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
##
##        M=scipy.sparse.construct.bmat( [ [M_diag_mm,Mhat_mA,None],
##                                         [Mhat_mA.T,Mhat_AA,None],
##                                         [CmP.T,None,MPP[SolvedDofP,:][:,SolvedDofP]] ] )

        
        F=omega**2*Freduced_F
        F  = scipy.append(F,scipy.zeros((len(SolvedDofA)+len(SolvedDofP)) , dtype='c16'))

        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype='c16')  , F , comm=mycomm )

        alpha_m    = sol[list(range(nb_mode_F))]
        P_A        = sol[list(range(nb_mode_F,nb_mode_F+len(SolvedDofA),1))]
        press      = eigen_vectors_F*alpha_m+Psi_FA*P_A
        enrichment = scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(nb_mode_F,nb_mode_F+len(SolvedDofA)))]
        CorrectedPressure=scipy.array(press)
        CorrectedPressure[SolvedDofA]=(CorrectedPressure[SolvedDofA].T+scipy.array(enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T))
        #frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,CorrectedPressure))
        frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes1,press.todense(),enrichment,LevelSet,LevelSetTangent))

        if rank==0:
            press_save.append(CorrectedPressure.real)
        
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes1,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

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
        print('Last eigenfrequency in fluid basis: ',freq_eigv_F[-1])




