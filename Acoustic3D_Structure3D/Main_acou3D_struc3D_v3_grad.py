###########################################################
## AIR CAVITY
## 3D RIGID STRUCTURE
## XFEM
###########################################################

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
#import silex_lib_dkt_fortran as silex_lib_dkt

#import silex_lib_tet4_fortran_test as silex_lib_tet4

#import silex_lib_porous_tet4_fortran

# manage mumps librairie
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
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7-v4.py
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

# gmsh -3 -format msh2 cavity_acou3D_struc_3D_v3_air.geo 
# gmsh -3 -format msh2 cavity_acou3D_struc_3D_v3_struc.geo 
 
mesh_file='geom/cavity_acou3D_struc_3D_v3'
results_file='results/cavity_acou3D_struc_3D_v3'

flag_write_gmsh_results=1

freq_ini     = 10.0
freq_end     = 120.0
nb_freq_step_per_proc=2

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=343.0 # ok
rho=1.21 # ok

##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_air.msh',3)
fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_air.msh',4,1) # air, cavity + controlled volume
fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_air.msh',4,5) # air, ONLY controlled volume

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem1   = fluid_elements1.shape[0]
fluid_nelem5   = fluid_elements5.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
fluid_nnodes5 = IdNodes5.shape[0]

fluid_ndof     = fluid_nnodes

if rank==0:
    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)
    print ("Number of nodes in air:",fluid_nnodes1)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_air_cavity_Mesh1',fluid_nodes,fluid_elements1,4)
    silex_lib_gmsh.WriteResults(results_file+'_air_controlled_volume_Mesh5',fluid_nodes,fluid_elements5,4)

##############################################################
# Load structure mesh
##############################################################
struc_nodes        = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh',3)
struc_elements,Idnodes_S_air_interface     = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',2,2)

struc_nnodes   = struc_nodes.shape[0]
struc_nelem    = struc_elements.shape[0]

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_struc_surface',struc_nodes,struc_elements,2)

if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)

##################################################################
# compute level set
##################################################################

tic = time.clock()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)

toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LS_signed_distance',fluid_nodes,fluid_elements1,4,[[[LevelSet],'nodal',1,'Level set']])
    silex_lib_gmsh.WriteResults2(results_file+'_LS_distance',fluid_nodes,fluid_elements1,4,[[[distance],'nodal',1,'Distance']])
    silex_lib_gmsh.WriteResults2(results_file+'_struc_air_interface',struc_nodes,struc_elements,2)

##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()

LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1,LevelSet)
LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]
silex_lib_gmsh.WriteResults2(results_file+'_LS_enriched_elements',fluid_nodes,fluid_elements1[LSEnrichedElements],4)
##EnrichedElements=LSEnrichedElements#[EnrichedElements-1]
LSEnrichednodes=scipy.unique(fluid_elements1[LSEnrichedElements])

tmp=[]
for i in LSEnrichednodes:
    for j in range(4):
        tmpp=scipy.where(fluid_elements1[:,j]==i)[0]
        for k in range(len(tmpp)):
            tmp.append(tmpp[k])
##    tmp.append(scipy.where(fluid_elements1[:,1]==i))
##    tmp.append(scipy.where(fluid_elements1[:,2]==i))
##    tmp.append(scipy.where(fluid_elements1[:,3]==i))

tmp=scipy.unique(scipy.array(tmp))
##tmp1,elttest0,tmp2=scipy.intersect1d(fluid_elements1[:,0],LSEnrichednodes,return_indices=True)
#silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements_test0',fluid_nodes,fluid_elements1[tmp],4)
#[75804, 97252, 97253,34973, 93135, 93137, 93248,83787, 93136,93525]
EnrichedElements0,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements1[tmp])
EnrichedElements0=scipy.unique(EnrichedElements0[list(range(NbEnrichedElements))])
EnrichedElements0=EnrichedElements0-1
EnrichedElements=tmp[EnrichedElements0]
toc = time.clock()
if rank==0:
    print ("time to find enriched elements:",toc-tic)

tic = time.clock()

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes,fluid_elements1[EnrichedElements],4)

LS_moins_enriched = scipy.setdiff1d(LSEnrichedElements,EnrichedElements)
enriched_moins_LS = scipy.setdiff1d(EnrichedElements,LSEnrichedElements)
silex_lib_gmsh.WriteResults2(results_file+'_LS_moins_enriched',fluid_nodes,fluid_elements1[LS_moins_enriched],4)
silex_lib_gmsh.WriteResults2(results_file+'_enriched_moins_LS',fluid_nodes,fluid_elements1[enriched_moins_LS],4)
##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=list(range(fluid_ndof))


##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.clock()

Enrichednodes = scipy.unique(fluid_elements1[EnrichedElements])

IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements1,fluid_nodes,LevelSet,celerity,rho)

KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofA=Enrichednodes-1

toc = time.clock()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

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
UF = scipy.zeros(2*fluid_ndof,dtype=float)
UF[9-1]=3.1250E-05

SolvedDof = scipy.hstack([SolvedDofF,SolvedDofA+fluid_ndof])

#################################################################
# Compute gradients with respect to parameters
##################################################################

print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)

LevelSet_gradient=-scipy.ones(fluid_nnodes)

IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_xfem_acou_tet4.globalacousticgradientmatrices(fluid_elements1[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient,LevelSet)
dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )


##############################################################
# FRF computation
##############################################################

Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    if rank==0:
        print('nb of total dofs: ',len(SolvedDofF)+len(SolvedDofA))

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("proc number",rank,"frequency=",freq)

        tic = time.clock()        
        
        F=scipy.array(omega**2*UF[SolvedDof] , dtype='c16')

        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , F )

        press1 = scipy.zeros((fluid_ndof),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        enrichment=scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=press1
        CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
        #frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
        frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,press1,enrichment,LevelSet,LevelSet*0-1.0))

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)
        
    frfsave=[scipy.array(frequencies),scipy.array(frf)]

    #comm.send(frfsave, dest=0, tag=11)

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

    # Save the FRF problem
    #Allfrequencies=scipy.zeros(nb_freq_step)
    #Allfrf=scipy.zeros(nb_freq_step)
    #k=0
    #if rank==0:
    #   for i in range(nproc):
    #        data = comm.recv(source=i, tag=11)
    #       for j in range(len(data[0])):
    #            Allfrequencies[k]=data[0][j]
    #            Allfrf[k]=data[1][j]
    #           k=k+1

        #Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        #Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
    f=open(results_file+'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()




