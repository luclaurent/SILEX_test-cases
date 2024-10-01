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

import silex_lib_xfem_acou_tet10
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

mesh_file='simple_cavity_tet10'
results_file='simple_cavity_tet10'

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

tic = time.clock()

fluid_nodes1    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',11,1) # air, cavity 
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',11,5) # air, ONLY controlled volume

fluid_nnodes   = fluid_nodes1.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
fluid_nelem1 = fluid_elements1.shape[0]
fluid_ndof1     = fluid_nnodes1

if rank==0:
    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes1,fluid_elements1,11)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet10.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

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
    silex_lib_gmsh.WriteResults2(results_file+'_LS_signed_distance',fluid_nodes1,fluid_elements1,11,[[[LevelSet],'nodal',1,'Level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_tangent_level_set',fluid_nodes1,fluid_elements1,11,[[[LevelSetTangent],'nodal',1,'Tangent level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_distance',fluid_nodes1,fluid_elements1,11,[[[distance],'nodal',1,'Distance']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_tangent_mesh',tangent_nodes,tangent_mesh,2)


##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()
# we take only the 4 corners nodes : it should be ok if tetra has straight lines
LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1[:,list(range(4))],LevelSet) 

LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]
EnrichedElements,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes1,fluid_elements1[LSEnrichedElements,:][:,list(range(4))])
EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.clock()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

tic = time.clock()
EdgeEnrichedElements,nbenrelts = silex_lib_xfem_acou_tet4.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes1,fluid_elements1[:,list(range(4))])
EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1
EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1[:,list(range(4))],LevelSetTangent)
EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])
toc = time.clock()
if rank==0:
    print ("time to find edge enriched elements:",toc-tic)

HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes1,fluid_elements1[LSEnrichedElements],11)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes1,fluid_elements1[EnrichedElements],11)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements',fluid_nodes1,fluid_elements1[EdgeEnrichedElements],11)
    silex_lib_gmsh.WriteResults2(results_file+'_edge_enriched_elements_in_all_mesh',fluid_nodes1,fluid_elements1[EdgeEnrichedElementsInAllMesh],11)
    silex_lib_gmsh.WriteResults2(results_file+'_heaviside_enriched_elements',fluid_nodes1,fluid_elements1[HeavisideEnrichedElements],11)

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.clock()

Enrichednodes = scipy.unique(fluid_elements1[EnrichedElements])
print(silex_lib_xfem_acou_tet10.globalxfemacousticmatrices.__doc__)
II,JJ,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet10.globalxfemacousticmatrices(fluid_elements1,fluid_nodes1,LevelSetTangent,LevelSet,celerity,rho)

KAA = scipy.sparse.csc_matrix( (Vaak,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAA= scipy.sparse.csc_matrix( (Vaam,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
KAF = scipy.sparse.csc_matrix( (Vafk,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAF = scipy.sparse.csc_matrix( (Vafm,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofA=Enrichednodes-1

toc = time.clock()
if rank==0:
    print ("time to compute enrichment:",toc-tic)

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

        tic = time.clock()        
        
        F=scipy.array(omega**2*UF[SolvedDof] , dtype='float')

        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float')  , F , comm=mycomm  )
        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=float) , scipy.array(F , dtype=float) )

        press1 = scipy.zeros((fluid_ndof1),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        enrichment=scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=scipy.zeros((fluid_ndof1),dtype=complex)
        CorrectedPressure[SolvedDofA]=press1[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
        frf.append(silex_lib_xfem_acou_tet10.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,press1))
        #frf.append(scipy.dot(scipy.dot(M,sol),sol))

        #print(silex_lib_xfem_acou_tet4.makexfemposfile.__doc__)
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,LevelSet,press1.real,enrichment.real,'press_plus.pos')
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,-LevelSet,press1.real,-enrichment.real,'press_moins.pos')

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)        

    frfsave=[scipy.array(frequencies),scipy.array(frf)]

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes1,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

    f=open(results_file+'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()



