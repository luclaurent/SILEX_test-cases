###########################################################
## AIR CAVITY
## IMPEDANCE PAROI
## THIN FLEXIBLE STRUCTURE
## XFEM
## CB REDUCTION FOR AIR
## MODAL REDUCTION FOR STRUCTURE
###########################################################

# good results / good management of reduced dofs on the air-porous interface

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

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

from SILEXlib import silex_lib_tet4_acou as xfemT4acou
from SILEXlib import silex_lib_gmsh as lgmsh
from SILEXlib import silex_lib_dkt_fortran as dkt

from SILEXlib import silex_lib_tet4_porous as porousT4

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
# mpirun -np 20 python3 Main_xfem_fluid_flex_struc_CB_reduction_impedance_paroi_12.py
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

mesh_file='geom/cavity'
results_file='results/cavity_with_impedance_air_flexible_structure_CB_reduction'

# fichier : results_file='results/cavity12_with_impedance_air_flexible_structure_CB_reduction'
# nb_mode_F = 400
# nb_mode_S = 70
# freq_ini     = 10.0
# freq_end     = 200.0
# nb_freq_step_per_proc=800


flag_write_gmsh_results=1

nb_mode_F = 100
nb_mode_S = 10
freq_ini     = 10.0
freq_end     = 200.0
nb_freq_step_per_proc=200

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

# structure damping
modal_damping_S=0.0

# impedance paroi : article Walid-JFD : CMAME 2008
d_imp_paroi = 50.0 # Pa.s/m
k_imp_paroi = 5.0e6 # Pa/m


##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
fluid_elements_S3,IdNodesS3 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,3)# air-porous interface surface
#fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
#fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

#for i in range(len(IdNodesS3)):
#    print(fluid_nodes[IdNodesS3[i]-1,:])

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
old = np.unique(fluid_elements1)
new = list(range(1,len(np.unique(fluid_elements1))+1))
new_nodes=fluid_nodes[np.unique(fluid_elements1)-1,:]

dico1 = dict(zip(old,new))
new_elements=np.zeros((fluid_nelem1,4),dtype=int)
for e in range(fluid_nelem1):
    for i in range(4):
        new_elements[e,i]=dico1[fluid_elements1[e][i]]

fluid_elements1 = new_elements
fluid_nodes1    = new_nodes

new_elements=np.zeros((fluid_nelem5,4),dtype=int)
for e in range(fluid_nelem5):
    for i in range(4):
        new_elements[e,i]=dico1[fluid_elements5[e][i]]

fluid_elements5 = new_elements

### renumbering porous
##old = np.unique(fluid_elements2)
##new = list(range(1,len(np.unique(fluid_elements2))+1))
##new_nodes=fluid_nodes[np.unique(fluid_elements2)-1,:]
##
##dico2 = dict(zip(old,new))
##new_elements=np.zeros((fluid_nelem2,4),dtype=int)
##for e in range(fluid_nelem2):
##    for i in range(4):
##        new_elements[e,i]=dico2[fluid_elements2[e][i]]
##
##fluid_elements2 = new_elements
##fluid_nodes2    = new_nodes
##
### renumbering porous external surfaces
##for i in range(len(IdNodesS4)):
##    IdNodesS4[i]=dico2[IdNodesS4[i]]

### Boundary conditions on air cavity
##IdNodesFixed_porous_us_x=IdNodesS4
####IdNodesFixed_porous_us_y=np.unique(np.hstack([IdNodesS4,IdNodesS6]))
##IdNodesFixed_porous_us_y=IdNodesS4
##IdNodesFixed_porous_us_z=IdNodesS4
##IdNodesFixed_porous_uf_x=IdNodesS4
##IdNodesFixed_porous_uf_y=IdNodesS4
##IdNodesFixed_porous_uf_z=IdNodesS4
##
##Fixed_Dofs_porous = np.hstack([(IdNodesFixed_porous_us_x-1)*6,
##                                  (IdNodesFixed_porous_us_y-1)*6+1,
##                                  (IdNodesFixed_porous_us_z-1)*6+2,
##                                  (IdNodesFixed_porous_uf_x-1)*6+3,
##                                  (IdNodesFixed_porous_uf_y-1)*6+4,
##                                  (IdNodesFixed_porous_uf_z-1)*6+5
##                                  ])

# get connectivity at air-porous interface
IdNodesS3_for_1=np.zeros(fluid_nnodes3,dtype=int)
#IdNodesS3_for_2=np.zeros(fluid_nnodes3,dtype=int)
for i in range(fluid_nnodes3):
    IdNodesS3_for_1[i]=dico1[IdNodesS3[i]] # for the air mesh
    #IdNodesS3_for_2[i]=dico2[IdNodesS3[i]] # for the porous mesh

InterfaceConnectivity=np.zeros((fluid_nelem3,3),dtype=int)
for e in range(fluid_nelem3):
    for i in range(3):
        InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]] # for the air mesh
        #InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]] # for the porous mesh

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,4)
    #silex_lib_gmsh.WriteResults(results_file+'_porous_material_Mesh2',fluid_nodes2,fluid_elements2,4)
    silex_lib_gmsh.WriteResults(results_file+'_porous_air_interface_Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
    #silex_lib_gmsh.WriteResults(results_file+'_porous_fixed_Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
    silex_lib_gmsh.WriteResults(results_file+'_air_controlled_volume_Mesh5',fluid_nodes1,fluid_elements5,4)

    #silex_lib_gmsh.WriteResults(results_file+'_air_porous_interface1',fluid_nodes1,InterfaceConnectivity[:,range(3)],2)
    #silex_lib_gmsh.WriteResults(results_file+'_air_porous_interface2',fluid_nodes2,InterfaceConnectivity[:,range(3,6,1)],2)

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

FixedStrucDof=np.hstack([FixedStrucDofUx,FixedStrucDofUy,FixedStrucDofUz,FixedStrucDofRx,FixedStrucDofRy,FixedStrucDofRz])

SolvedDofS=np.setdiff1d(range(struc_ndof),FixedStrucDof)

FS=np.zeros(struc_ndof)
#IddofLoadStructure=[(IdNodeLoadStructure-1)*6+2]
#FS[IddofLoadStructure]=1.0

##############################################################
# Compute structure matrices
##############################################################
tic = time.process_time()

IIks,JJks,Vks,Vms=silex_lib_dkt.stiffnessmatrix(struc_nodes,struc_elements,material_Struc)

KSS = scipy.sparse.csc_matrix( (Vks,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )
MSS = scipy.sparse.csc_matrix( (Vms,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )

toc = time.process_time()
if rank==0:
    print ("time for computing structure:",toc-tic)

##################################################################
# Compute eigen modes of the structure
##################################################################
tic = time.process_time()

eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],nb_mode_S,MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*np.pi))

toc = time.process_time()
if rank==0:
    print ("time for computing the structure modal basis:",toc-tic)

eigen_vector_S_list=[]
for i in range(nb_mode_S):
    Q=np.zeros(struc_ndof)
    Q[SolvedDofS]=eigen_vectors_S[:,i]
    disp=np.zeros((struc_nnodes,3))
    disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
    disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
    disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
    eigen_vector_S_list.append(disp)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_structure_modes',struc_nodes,struc_elements,2,[[eigen_vector_S_list,'nodal',3,'modes']])

##################################################################
# compute level set
##################################################################

tic = time.process_time()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes1,struc_nodes,struc_elements)

toc = time.process_time()
if rank==0:
    print ("time to compute level set:",toc-tic)

tic = time.process_time()

tangent_nodes,tangent_mesh=silex_lib_xfem_acou_tet4.buildtangentedgemesh(struc_nodes,struc_elements,struc_boun)

LevelSetTangent,tmp = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes1,tangent_nodes,tangent_mesh)

toc = time.process_time()
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
tic = time.process_time()

LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSet)
LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]
EnrichedElements,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes1,fluid_elements1[LSEnrichedElements])
EnrichedElements=np.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.process_time()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

tic = time.process_time()

EdgeEnrichedElements,nbenrelts = silex_lib_xfem_acou_tet4.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes1,fluid_elements1)
EdgeEnrichedElements=np.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1

EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSetTangent)
EdgeEnrichedElementsInAllMesh=np.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])

toc = time.process_time()
if rank==0:
    print ("time to find edge enriched elements:",toc-tic)

HeavisideEnrichedElements=np.setdiff1d(EnrichedElements,EdgeEnrichedElements)

AllElementsExceptEdgeEnrichedElements=np.setdiff1d(range(fluid_nelem1),EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes1,fluid_elements1[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes1,fluid_elements1[EnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements',fluid_nodes1,fluid_elements1[EdgeEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements_in_all_mesh',fluid_nodes1,fluid_elements1[EdgeEnrichedElementsInAllMesh],4)
    silex_lib_gmsh.WriteResults2(results_file+'_heaviside_enriched_elements',fluid_nodes1,fluid_elements1[HeavisideEnrichedElements],4)

##################################################################
# Compute coupling STRUCTURE / AIR terms
##################################################################
tic = time.process_time()

IIc1,JJc1,Vc1=silex_lib_xfem_acou_tet4.computexfemcoupling1(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements)
IIc2,JJc2,Vc2=silex_lib_xfem_acou_tet4.computexfemcoupling2(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements,LevelSet)

CSA=0.5*scipy.sparse.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(struc_ndof,fluid_ndof1) )+0.5*scipy.sparse.csc_matrix( (Vc2,(IIc2,JJc2)), shape=(struc_ndof,fluid_ndof1) )

toc = time.process_time()
if rank==0:
    print ("time to compute coupling matrices:",toc-tic)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

SolvedDofB=np.hstack([IdNodesS3_for_1-1,9-1]) # 9 : node number where acoustic source is imposed
SolvedDofI=np.setdiff1d(SolvedDofF,SolvedDofB)

##############################################################
# Compute impedance matrices
##############################################################
#print(silex_lib_porous_tet4_fortran.stiffnessimpedanceparoi.__doc__)
IIimp,JJimp,Vimp=silex_lib_porous_tet4_fortran.stiffnessimpedanceparoi(fluid_nodes1,InterfaceConnectivity)

CII=scipy.sparse.csc_matrix( (Vimp,(IIimp,JJimp)), shape=(fluid_ndof1,fluid_ndof1) )


##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.process_time()

#Enrichednodes = np.unique(fluid_elements1[HeavisideEnrichedElements])
Enrichednodes = np.unique(fluid_elements1[EnrichedElements])

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

toc = time.process_time()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

##################################################################
# Compute Edge enrichment
##################################################################
tic = time.process_time()

PartiallyPositiveLStgtElements=np.hstack([PositiveLStgtElements,EdgeEnrichedElementsInAllMesh])

#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],LevelSet,LevelSetTangent,celerity,rho)
II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[PartiallyPositiveLStgtElements],LevelSet,LevelSetTangent,celerity,rho)
#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1,LevelSet,LevelSetTangent,celerity,rho)

KAAedge = scipy.sparse.csc_matrix( (vkaa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAAedge = scipy.sparse.csc_matrix( (vmaa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
KAFedge = scipy.sparse.csc_matrix( (vkfa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAFedge = scipy.sparse.csc_matrix( (vmfa,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )

toc = time.process_time()
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
# Compute eigen modes of the fluid: internal dof I
##################################################################
tic = time.process_time()

if rank==0:
    print ("Compute fluid modes")
eigen_values_I,eigen_vectors_I= scipy.sparse.linalg.eigsh(KFF[SolvedDofI,:][:,SolvedDofI],nb_mode_F,MFF[SolvedDofI,:][:,SolvedDofI],sigma=0,which='LM')

freq_eigv_I=list(scipy.sqrt(eigen_values_I)/(2*np.pi))
eigen_vector_F_list=[]
for i in range(nb_mode_F):
    tmp=np.zeros((fluid_ndof1) , dtype='float')
    tmp[SolvedDofI]=eigen_vectors_I[:,i].real
    eigen_vector_F_list.append(tmp)

if rank==0:
    print ("LAST fluid eigen frequencies : ",freq_eigv_I[-1])
 
if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes1,fluid_elements1,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

toc = time.process_time()
if rank==0:
    print ("time for computing the fluid modes:",toc-tic)


##################################################################
# Compute Psi_IA for the fluid: Psi_IA = - KII^{-1} * KIA
##################################################################
tic = time.process_time()

if rank==0:
    print ("Compute PSI_IA")

omega_cst=1.0*2.0*np.pi
MySolve = scipy.sparse.linalg.factorized( KFF[SolvedDofI,:][:,SolvedDofI]-(omega_cst**2)*MFF[SolvedDofI,:][:,SolvedDofI] ) # Makes LU decomposition.

if rank==0:
    print("LU decomposition has been made")

if rank==0:
    Psi_IA=np.zeros((len(SolvedDofI),len(SolvedDofA)))

i=0
j=0
while j<len(SolvedDofA):
    j=i+rank
    if j<len(SolvedDofA):
        One_dof=[SolvedDofA[j]]
        KIA_i_column=-KAF[SolvedDofI,:][:,One_dof]
        #tmp=np.zeros(len(SolvedDofF))
        #tmp[SolvedDofI]=KIA_i_column.todense()
        Xi=MySolve( KIA_i_column.todense() )
        if rank!=0:
            comm.send([Xi,j], dest=0, tag=11)
        if rank==0:
            for k in range(nproc):
                if k!=0:
                    if i+k<len(SolvedDofA):
                        [Xi,j]=comm.recv(source=k, tag=11)
                        Psi_IA[:,j]=np.array(Xi)[:,0]
                else:
                    Psi_IA[:,j]=np.array(Xi)[:,0]
        i=i+nproc

if rank==0:
    print("End of PSI_IA computing, send to other proc.")

if rank==0:
    for i in range(nproc):
        if i!=0:
            comm.send(Psi_IA, dest=i, tag=11)
if rank!=0:
    Psi_IA=comm.recv(source=0, tag=11)

Psi_IA=scipy.sparse.csc_matrix(Psi_IA)
toc = time.process_time()

if rank==0:
    print ("time to compute PSI_IA:",toc-tic)

##if rank==0:
##    eigen_vector_I_list=[]
##    for i in range(len(SolvedDofA)):
##        tmp=np.zeros((fluid_ndof1) , dtype='float')
##        tmp[SolvedDofI]=Psi_IA[:,i].todense()
##        eigen_vector_I_list.append(tmp)
##
##if (flag_write_gmsh_results==1) and (rank==0):
##    silex_lib_gmsh.WriteResults2(results_file+'_Psi_IA',fluid_nodes1,fluid_elements1,4,[[eigen_vector_I_list,'nodal',1,'PSI IA pressure']])


##################################################################
# Compute Psi_IB for the fluid: Psi_IB = - KII^{-1} * KIB
##################################################################
tic = time.process_time()

if rank==0:
    print ("Compute PSI_IB")

if rank==0:
    Psi_IB=np.zeros((len(SolvedDofI),len(SolvedDofB)))

i=0
j=0
while j<len(SolvedDofB):
    j=i+rank
    if j<len(SolvedDofB):
        One_dof=[SolvedDofB[j]]
        KIB_i_column=-KFF[SolvedDofI,:][:,One_dof]
        Xi=MySolve( KIB_i_column.todense() )
        if rank!=0:
            comm.send([Xi,j], dest=0, tag=11)
        if rank==0:
            for k in range(nproc):
                if k!=0:
                    if i+k<len(SolvedDofB):
                        [Xi,j]=comm.recv(source=k, tag=11)
                        Psi_IB[:,j]=np.array(Xi)[:,0]
                else:
                    Psi_IB[:,j]=np.array(Xi)[:,0]
        i=i+nproc

if rank==0:
    print("End of PSI_IB computing, send to other proc.")

if rank==0:
    for i in range(nproc):
        if i!=0:
            comm.send(Psi_IB, dest=i, tag=11)
if rank!=0:
    Psi_IB=comm.recv(source=0, tag=11)

Psi_IB=scipy.sparse.csc_matrix(Psi_IB)
toc = time.process_time()

if rank==0:
    print ("time to compute PSI_FB:",toc-tic)

##if rank==0:
##    eigen_vector_I_list=[]
##    for i in range(len(SolvedDofB)):
##        tmp=np.zeros((fluid_ndof1) , dtype='float')
##        tmp[SolvedDofI]=Psi_IB[:,i].todense()
##        eigen_vector_I_list.append(tmp)
##
##if (flag_write_gmsh_results==1) and (rank==0):
##    silex_lib_gmsh.WriteResults2(results_file+'_Psi_IB',fluid_nodes1,fluid_elements1,4,[[eigen_vector_I_list,'nodal',1,'PSI IB pressure']])


##################################################################
# Construct structure projection
##################################################################

VK_diag_nn = eigen_values_S
VM_diag_nn = eigen_values_S/eigen_values_S
IIDnn = list(range(nb_mode_S))
JJDnn = list(range(nb_mode_S))

Knn= scipy.sparse.csc_matrix( (VK_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )
Mnn= scipy.sparse.csc_matrix( (VM_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )

# Projection of air-structure coupling on structure modal basis
CnA = eigen_vectors_S.T*CSA[SolvedDofS,:][:,SolvedDofA]


##################################################################
# Compute structure damping matrix
##################################################################
VDnn = 2.0*modal_damping_S*scipy.sqrt(eigen_values_S)
IIDnn = list(range(nb_mode_F+len(SolvedDofB)+len(SolvedDofA),nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S))
JJDnn = list(range(nb_mode_F+len(SolvedDofB)+len(SolvedDofA),nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S))

totaldofs=nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S

D = scipy.sparse.coo_matrix( (VDnn,(IIDnn,JJDnn)), shape=(totaldofs,totaldofs) )



##################################################################
# Construct the whole system
##################################################################

# Fluid part
tic = time.process_time()

VK_diag_mm = eigen_values_I
VM_diag_mm = eigen_values_I/eigen_values_I
IIDmm = list(range(nb_mode_F))
JJDmm = list(range(nb_mode_F))


#PhiFm=eigen_vectors_F

K_diag_mm= scipy.sparse.csc_matrix( (VK_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
M_diag_mm= scipy.sparse.csc_matrix( (VM_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
Khat_AA = KAA[SolvedDofA,:][:,SolvedDofA]+Psi_IA.T*KAF[SolvedDofI,:][:,SolvedDofA]
Khat_BB = KFF[SolvedDofB,:][:,SolvedDofB]+Psi_IB.T*KFF[SolvedDofI,:][:,SolvedDofB]

Khat_BA = KFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA

toc = time.process_time()
if rank==0:
    print ("time to compute Khat_hat:",toc-tic)

tic = time.process_time()

Mstar_IA = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IA+MAF[SolvedDofI,:][:,SolvedDofA]
Mstar_IB = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IB+MFF[SolvedDofI,:][:,SolvedDofB]

toc = time.process_time()
if rank==0:
    print ("time to compute Mstar:",toc-tic)

tic = time.process_time()

Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.sparse.csc_matrix((Psi_IA.T).todense()*Mstar_IA.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IB.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

#Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+Psi_IA.T*Mstar_IA+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
#Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+Psi_IB.T*Mstar_IB+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

#Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.dot((Psi_IA.T).todense(),Mstar_IA.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
#Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+scipy.dot((Psi_IB.T).todense(),Mstar_IB.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

toc = time.process_time()
if rank==0:
    print ("time to compute Mhat:",toc-tic)

#CmP=PhiFm.T*CPF[SolvedDofP,:][:,SolvedDofF].T
#CAP=Psi_FA.T*CPF[SolvedDofP,:][:,SolvedDofF].T

eigen_vectors_I=scipy.sparse.csc_matrix(eigen_vectors_I)

Mhat_mA = eigen_vectors_I.T*Mstar_IA
Mhat_mB = eigen_vectors_I.T*Mstar_IB

Mhat_BA = scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IA.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA


##################################################################
# Build Second member
##################################################################

# To impose the load on the fluid:
# fluid node number 1
UF = np.zeros(fluid_ndof1,dtype=float)
UF[9-1]=3.1250E-05

#SolvedDof = np.hstack([SolvedDofF,SolvedDofA+fluid_ndof1,SolvedDofP+2*fluid_ndof1])
Freduced_F=np.hstack([np.zeros(nb_mode_F,dtype=float),UF[SolvedDofB]])

##############################################################
# FRF computation
##############################################################


Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    if rank==0:
        print('nb of total dofs: ',nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S,)

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*np.pi*freq

        print ("proc number",rank,"frequency=",freq)

        tic = time.process_time()

        Kimp=-(omega**2/(k_imp_paroi-1j*omega*d_imp_paroi))*CII

        K=scipy.sparse.construct.bmat( [ [K_diag_mm,None,       None,                           None],
                                         [None,     Khat_BB+Kimp[SolvedDofB,:][:,SolvedDofB],    Khat_BA,                           None],
                                         [None,     Khat_BA.T,  Khat_AA,                           None],
                                         [None,     None,       -CnA,Knn]
                                         ]
                                       )
        
        M=scipy.sparse.construct.bmat( [ [M_diag_mm,    Mhat_mB,    Mhat_mA,                          None],
                                         [Mhat_mB.T,    Mhat_BB,    Mhat_BA,                               None],
                                         [Mhat_mA.T,    Mhat_BA.T,  Mhat_AA,                               CnA.T],
                                         [None,         None,       None,                              Mnn]
                                         ]
                                       )
        
        F  = scipy.append(omega**2*Freduced_F,np.zeros((len(SolvedDofA)+nb_mode_S) , dtype='c16'))

        sol = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M+omega*D*1j,dtype='c16')  , F , comm=mycomm )

        alpha_m    = sol[list(range(nb_mode_F))]

        P_B        = sol[list(range(nb_mode_F,nb_mode_F+len(SolvedDofB),1))]
        P_A        = sol[list(range(nb_mode_F+len(SolvedDofB),nb_mode_F+len(SolvedDofB)+len(SolvedDofA),1))]
        
        P_I        = eigen_vectors_I*alpha_m+Psi_IB*P_B+Psi_IA*P_A
        
        press      = np.zeros(fluid_ndof1,dtype=complex)
        press[SolvedDofI] = P_I
        press[SolvedDofB] = P_B
        enrichment = np.zeros(fluid_ndof1,dtype=complex)
        enrichment[SolvedDofA]= P_A
        CorrectedPressure=np.array(press)
        CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA].T+np.array(enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,CorrectedPressure))
        #frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes1,press,enrichment,LevelSet,LevelSetTangent))

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)
            alpha_n=sol[list(range(nb_mode_F+len(SolvedDofB)+len(SolvedDofA),nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S,1))].real
            Q=np.zeros((struc_ndof),dtype=float)
            Q[SolvedDofS]=scipy.dot(eigen_vectors_S,alpha_n)
            disp=np.zeros((struc_nnodes,3))
            disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
            disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
            disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
            disp_save.append(disp)
            

    frfsave=[frequencies,frf]
    if rank!=0:
        comm.send(frfsave, dest=0, tag=11)

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes1,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_struct_frf',struc_nodes,struc_elements,2,[[disp_save,'nodal',3,'displacement']])

    # Save the FRF problem
    Allfrequencies=np.zeros(nb_freq_step)
    Allfrf=np.zeros(nb_freq_step)
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
        Allfrfsave=[np.array(list(Allfrequencies)),np.array(list(Allfrf))]
        f=open(results_file+'_results.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()
        print('Last eigenfrequency in fluid basis: ',freq_eigv_I[-1])
        print('Last eigenfrequency in structure basis: ',freq_eigv_S[-1])
        print('nb of total dofs: ',nb_mode_F+len(SolvedDofB)+len(SolvedDofA)+nb_mode_S,)




