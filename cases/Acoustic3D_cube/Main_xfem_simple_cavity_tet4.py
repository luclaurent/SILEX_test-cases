###########################################################
## AIR CAVITY
## POROUS MATERIAL
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
sys.path.append('../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_gmsh
import silex_lib_dkt_fortran as silex_lib_dkt

import silex_lib_porous_tet4_fortran
#import silex_lib_tet4_fortran
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
# mpirun -np 20 python3 Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
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

mesh_file='simple_cavity_tet4'
results_file='simple_cavity_tet4'

freq_ini     = 150.0
freq_end     = 500.0
nb_freq_step = 500

flag_write_gmsh_results=1

# air
celerity=343.0 # ok
rho=1.21 # ok

# shell structure
material_Struc=[]
material_Struc.append(75000.0e6) # E Young
material_Struc.append(0.33) # nu
material_Struc.append(5.0e-3) # thickness
material_Struc.append(2700.0) # rho

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity 
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume

fluid_nnodes   = fluid_nodes.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
fluid_nelem1 = fluid_elements1.shape[0]
fluid_ndof1     = fluid_nnodes1

if rank==0:
    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)


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

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,4)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

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
EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.process_time()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

tic = time.process_time()
EdgeEnrichedElements,nbenrelts = silex_lib_xfem_acou_tet4.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes1,fluid_elements1)
EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1
EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSetTangent)
EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])
toc = time.process_time()
if rank==0:
    print ("time to find edge enriched elements:",toc-tic)

HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes1,fluid_elements1[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes1,fluid_elements1[EnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements',fluid_nodes1,fluid_elements1[EdgeEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements_in_all_mesh',fluid_nodes1,fluid_elements1[EdgeEnrichedElementsInAllMesh],4)
    silex_lib_gmsh.WriteResults2(results_file+'_heaviside_enriched_elements',fluid_nodes1,fluid_elements1[HeavisideEnrichedElements],4)

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.process_time()

Enrichednodes = scipy.unique(fluid_elements1[EnrichedElements])

IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements1[EnrichedElements],fluid_nodes1,LevelSet,celerity,rho)

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

#PartiallyPositiveLStgtElements=scipy.hstack([PositiveLStgtElements,EdgeEnrichedElementsInAllMesh])
#print(silex_lib_xfem_acou_tet4.computeedgeenrichment2.__doc__)
#II,JJ,vkaa,vmaa,vkfa,vmfa,NodesCutEdgeMesh,NbNodesCutEdgeMesh,EltCutEdgeMEsh,NbEltCutEdgeMEsh= silex_lib_xfem_acou_tet4.computeedgeenrichment2(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],LevelSet,LevelSetTangent,celerity,rho)
#print(silex_lib_xfem_acou_tet4.cutedgeglobalmesh.__doc__)
#NodesCutEdgeMesh,NbNodesCutEdgeMesh,EltCutEdgeMEsh,NbEltCutEdgeMEsh,LSnewnodes,LStgtnewnodes= silex_lib_xfem_acou_tet4.cutedgeglobalmesh(fluid_nodes1,fluid_elements1[EdgeEnrichedElements[0:1]],LevelSet,LevelSetTangent)#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[PartiallyPositiveLStgtElements],LevelSet,LevelSetTangent,celerity,rho)
NodesCutEdgeMesh,NbNodesCutEdgeMesh,EltCutEdgeMEsh,NbEltCutEdgeMEsh,LSnewnodes,LStgtnewnodes= silex_lib_xfem_acou_tet4.cutedgeglobalmesh(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],LevelSet,LevelSetTangent)#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1[PartiallyPositiveLStgtElements],LevelSet,LevelSetTangent,celerity,rho)
NodesCutEdgeMesh=NodesCutEdgeMesh[range(NbNodesCutEdgeMesh)]
EltCutEdgeMEsh=EltCutEdgeMEsh[range(NbEltCutEdgeMEsh)]
LSnewnodes=LSnewnodes[range(NbNodesCutEdgeMesh)]
LStgtnewnodes=LStgtnewnodes[range(NbNodesCutEdgeMesh)]
#silex_lib_gmsh.WriteResults(results_file+'_CutEdgeMEsh',NodesCutEdgeMesh,EltCutEdgeMEsh,4)
silex_lib_gmsh.WriteResults2(results_file+'_CutEdgeMEsh_LS',NodesCutEdgeMesh,EltCutEdgeMEsh,4,[[[LSnewnodes],'nodal',1,'Level set'],[[LStgtnewnodes],'nodal',1,'Tangent level set']])
#silex_lib_gmsh.WriteResults2(results_file+'_CutEdgeMEsh_LS_tangent_level_set',NodesCutEdgeMesh,EltCutEdgeMEsh,4,[[[LStgtnewnodes],'nodal',1,'Tangent level set']])
#STOP
#II,JJ,vkaa,vmaa,vkfa,vmfa= silex_lib_xfem_acou_tet4.computeedgeenrichment(NodesCutEdgeMesh,EltCutEdgeMEsh,LevelSet,LevelSetTangent,celerity,rho)
#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment(fluid_nodes1,fluid_elements1,LevelSet,LevelSetTangent,celerity,rho)

# calcul du volume pour tester (en activant les sorties shell du volume :) !!)
#Ik,Jk,Vk=silex_lib_tet4_fortran.stiffnessmatrix(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],[1,1])

II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment2(fluid_nodes1,fluid_elements1[EdgeEnrichedElements],LevelSet,LevelSetTangent,celerity,rho)
#II,JJ,vkaa,vmaa,vkfa,vmfa = silex_lib_xfem_acou_tet4.computeedgeenrichment2(fluid_nodes1,fluid_elements1[EdgeEnrichedElements[0:1]],LevelSet,LevelSetTangent,celerity,rho)

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

##################################################################
# Construct the whole system
#################################################################

K=scipy.sparse.construct.bmat( [
            [KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
            [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]] )

M=scipy.sparse.construct.bmat( [
            [MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
            [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]] )

##################################################################
# Build Second member
##################################################################

# To impose the load on the fluid:
# fluid node number 1
UF = scipy.zeros(2*fluid_ndof1,dtype=float)
UF[1-1]=3.1250E-05

SolvedDof = scipy.hstack([SolvedDofF,SolvedDofA+fluid_ndof1])

##############################################################
# FRF computation
##############################################################

Flag_frf_analysis=1
frequencies=[]
frf=[]
frfgradient=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    if rank==0:
        print('nb of total dofs: ',len(SolvedDofF)+len(SolvedDofA))

    press_save=[]
    disp_save=[]
    
    #for i in range(nb_freq_step):
    for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

        #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq
        #print('omega = ',omega)

        print ("proc number",rank,"frequency=",freq)

        tic = time.process_time()        
        
        F=scipy.array(omega**2*UF[SolvedDof] , dtype='float')

        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float')  , F , comm=mycomm  )
        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=float) , scipy.array(F , dtype=float) )

        press1 = scipy.zeros((fluid_ndof1),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        enrichment=scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=scipy.zeros((fluid_ndof1),dtype=complex)
        CorrectedPressure[SolvedDofA]=press1[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
        #frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
        frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,press1,enrichment,LevelSet,LevelSet*0-1.0))
        #frf.append(scipy.dot(scipy.dot(M,sol),sol))

        #print(silex_lib_xfem_acou_tet4.makexfemposfile.__doc__)
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,LevelSet,press1.real,enrichment.real,'press_plus.pos')
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,-LevelSet,press1.real,-enrichment.real,'press_moins.pos')

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)        

    frfsave=[scipy.array(frequencies),scipy.array(frf)]

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

    f=open(results_file+'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()



