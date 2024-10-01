###########################################################
## AIR CAVITY
## IMPEDANCE PAROI
## THIN FLEXIBLE STRUCTURE
## XFEM
###########################################################

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
import silex_lib_dkt_fortran as silex_lib_dkt

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
# mpirun -np 20  python3 Main_xfem_fluid_porous_flex_struc_no_reduction_8.py
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

mesh_file='geom/cavity12_with_porous_air'
results_file='results/cavity12_with_impedance_air_flexible_structure'

flag_write_gmsh_results=1

freq_ini     = 10.0
freq_end     = 200.0
nb_freq_step_per_proc=800

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=343.0 # ok
rho=1.21 # ok

# shell structure
material_Struc=[]
material_Struc.append(75000.0e6) # E Young
material_Struc.append(0.33) # nu
material_Struc.append(5.0e-3) # thickness
material_Struc.append(2700.0) # rho


# impedance paroi : article Walid-JFD : CMAME 2008
d_imp_paroi = 50.0 # Pa.s/m
k_imp_paroi = 5.0e6 # Pa/m


##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
fluid_elements_S3,IdNodesS3 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,3)# air-porous interface surface : impedance paroi
#fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
#fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem1   = fluid_elements1.shape[0]
#fluid_nelem2   = fluid_elements2.shape[0]
fluid_nelem3   = fluid_elements_S3.shape[0]
#fluid_nelem4   = fluid_elements_S4.shape[0]
fluid_nelem5   = fluid_elements5.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
#fluid_nnodes2 = IdNodes2.shape[0]
fluid_nnodes3 = IdNodesS3.shape[0]
#fluid_nnodes4 = IdNodesS4.shape[0]
fluid_nnodes5 = IdNodes5.shape[0]

fluid_ndof1     = fluid_nnodes1

if rank==0:
    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)

    print ("Number of nodes in air:",fluid_nnodes1)
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

### renumbering porous
##old = scipy.unique(fluid_elements2)
##new = list(range(1,len(scipy.unique(fluid_elements2))+1))
##new_nodes=fluid_nodes[scipy.unique(fluid_elements2)-1,:]
##
##dico2 = dict(zip(old,new))
##new_elements=scipy.zeros((fluid_nelem2,4),dtype=int)
##for e in range(fluid_nelem2):
##    for i in range(4):
##        new_elements[e,i]=dico2[fluid_elements2[e][i]]
##
##fluid_elements2 = new_elements
##fluid_nodes2    = new_nodes

### renumbering porous external surfaces
##for i in range(len(IdNodesS4)):
##    IdNodesS4[i]=dico2[IdNodesS4[i]]

### Boundary conditions on air cavity
##IdNodesFixed_porous_us_x=IdNodesS4
####IdNodesFixed_porous_us_y=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS6]))
##IdNodesFixed_porous_us_y=IdNodesS4
##IdNodesFixed_porous_us_z=IdNodesS4
##IdNodesFixed_porous_uf_x=IdNodesS4
##IdNodesFixed_porous_uf_y=IdNodesS4
##IdNodesFixed_porous_uf_z=IdNodesS4
##
##Fixed_Dofs_porous = scipy.hstack([(IdNodesFixed_porous_us_x-1)*6,
##                                  (IdNodesFixed_porous_us_y-1)*6+1,
##                                  (IdNodesFixed_porous_us_z-1)*6+2,
##                                  (IdNodesFixed_porous_uf_x-1)*6+3,
##                                  (IdNodesFixed_porous_uf_y-1)*6+4,
##                                  (IdNodesFixed_porous_uf_z-1)*6+5
##                                  ])

# get connectivity at air- impedance paroi
IdNodesS3_for_1=scipy.zeros(fluid_nnodes3,dtype=int)
#IdNodesS3_for_2=scipy.zeros(fluid_nnodes3,dtype=int)
for i in range(fluid_nnodes3):
    IdNodesS3_for_1[i]=dico1[IdNodesS3[i]] # for the air mesh
    #IdNodesS3_for_2[i]=dico2[IdNodesS3[i]] # for the porous mesh

InterfaceConnectivity=scipy.zeros((fluid_nelem3,3),dtype=int)
for e in range(fluid_nelem3):
    for i in range(3):
        InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]] # for the air mesh
##        InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]] # for the porous mesh

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,4)
    #silex_lib_gmsh.WriteResults(results_file+'_porous_material_Mesh2',fluid_nodes2,fluid_elements2,4)
    silex_lib_gmsh.WriteResults(results_file+'_porous_air_interface_Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
    #silex_lib_gmsh.WriteResults(results_file+'_porous_fixed_Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
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

FS=scipy.zeros(struc_ndof)
#IddofLoadStructure=[(IdNodeLoadStructure-1)*6+2]
#FS[IddofLoadStructure]=1.0

##############################################################
# Compute structure matrices
##############################################################
tic = time.clock()

IIks,JJks,Vks,Vms=silex_lib_dkt.stiffnessmatrix(struc_nodes,struc_elements,material_Struc)

KSS = scipy.sparse.csc_matrix( (Vks,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )
MSS = scipy.sparse.csc_matrix( (Vms,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )

toc = time.clock()
if rank==0:
    print ("time for computing structure:",toc-tic)


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

##if rank==0:
##    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes1,fluid_elements1[EnrichedElements],4)


##################################################################
# Compute coupling STRUCTURE / AIR terms
##################################################################
tic = time.clock()

IIc1,JJc1,Vc1=silex_lib_xfem_acou_tet4.computexfemcoupling1(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements)
IIc2,JJc2,Vc2=silex_lib_xfem_acou_tet4.computexfemcoupling2(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements,LevelSet)

CSA=0.5*scipy.sparse.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(struc_ndof,fluid_ndof1) )+0.5*scipy.sparse.csc_matrix( (Vc2,(IIc2,JJc2)), shape=(struc_ndof,fluid_ndof1) )

toc = time.clock()
if rank==0:
    print ("time to compute coupling matrices:",toc-tic)


##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

##############################################################
# Compute impedance matrices
##############################################################
print(silex_lib_porous_tet4_fortran.stiffnessimpedanceparoi.__doc__)
IIimp,JJimp,Vimp=silex_lib_porous_tet4_fortran.stiffnessimpedanceparoi(fluid_nodes1,InterfaceConnectivity)

CII=scipy.sparse.csc_matrix( (Vimp,(IIimp,JJimp)), shape=(fluid_ndof1,fluid_ndof1) )


##############################################################
# Compute Coupling Porous-air Matrices
##############################################################

#print(silex_lib_porous_tet4_fortran.computecouplingporousair.__doc__)

##IIpf,JJpf,Vpf=silex_lib_porous_tet4_fortran.computecouplingporousair(fluid_nodes1,InterfaceConnectivity,po_por)
##CPF=scipy.sparse.csc_matrix( (Vpf,(IIpf,JJpf)), shape=(fluid_ndof2,fluid_ndof1) )

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
UF = scipy.zeros(2*fluid_ndof1+struc_ndof,dtype=float)
UF[9-1]=3.1250E-05

SolvedDof = scipy.hstack([SolvedDofF,SolvedDofA+fluid_ndof1,SolvedDofS+2*fluid_ndof1])

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

        #IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)
        #KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
        #MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
        Kimp=-(omega**2/(k_imp_paroi-1j*omega*d_imp_paroi))*CII

        K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF]+Kimp[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA],None],
                                         [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA],None],
                                         [None,       -CSA[SolvedDofS,:][:,SolvedDofA],KSS[SolvedDofS,:][:,SolvedDofS]]
                                         ] )
        
        M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA],None],
                                         [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA],CSA[SolvedDofS,:][:,SolvedDofA].T],
                                         [None,        None,                                  MSS[SolvedDofS,:][:,SolvedDofS]]] )

        F=scipy.array(omega**2*UF[SolvedDof] , dtype='c16')
        
        sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16') , F , comm=mycomm )

        press1 = scipy.zeros((fluid_ndof1),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]

        enrichment=scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=press1
        CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,CorrectedPressure))

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)
            
    frfsave=[frequencies,frf]
    if rank!=0:
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
            if i==0:
                data = frfsave
            else:
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
        print('nb of total dofs: ',len(SolvedDofF)+len(SolvedDofA)+len(SolvedDofS),)


