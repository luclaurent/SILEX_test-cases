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
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io
import getopt
from pathlib import Path
from loguru import logger

import pylab as pl
import pickle


import mumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
acousticsFEM = silex_lib_fem.LinearAcousticsTET10()
acousticsXFEM = silex_lib_xfem.LinearAcousticsTET10()
structureFEM = silex_lib_fem.DKT()


from mpi4py import MPI

comm = MPI.COMM_WORLD

nproc = comm.Get_size()
rank = comm.Get_rank()
log_format = (
    "<cyan> R{extra[rank]}</cyan> |"
    "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> | "
    "<level>{message}</level>"
)

logger.remove()
logger.configure(extra={"rank": 0})  # Default values
logger.add(
    sys.stdout,
    level="DEBUG",
    format=log_format,
    colorize=True,
    backtrace=True,
    diagnose=True,
)
logger = logger.bind(rank=rank)

# mpirun -np 2 python Main_xfem.py
logger.info("START")


def mpiInfo():
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    return nproc, rank, comm


class comm_mumps_one_proc:
    rank = 0

    def py2f(self):
        return 0


mycomm = comm_mumps_one_proc()

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
    logger.info("time at the beginning of the computation: {}".format(time.ctime()))

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

##############################################################
# Datas
##############################################################

mesh_file='geom/simple_cavity_tet10'
results_file='results/simple_cavity_tet10'
cwd = Path(__file__).resolve().parent

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

# start reading mesh
mesh = msh.mshReader(cwd / (mesh_file + "_fluid.msh"))

fluid_nodes1 = mesh.getNodes()
fluid_elements1 = mesh.getElements(tag=1)["TET10"]  # air, cavity
fluid_elements5 = mesh.getElements(tag=5)["TET10"]  # air, control volum

IdNodes1 = np.unique(fluid_elements1.flatten())
IdNodes5 = np.unique(fluid_elements5.flatten())

fluid_nnodes = fluid_nodes1.shape[0]
fluid_nnodes1 = IdNodes1.shape[0]
fluid_nelem1 = fluid_elements1.shape[0]
fluid_ndof1 = fluid_nnodes1

if rank==0:
    logger.info("Number of nodes:",fluid_nnodes)
    logger.info("Number of elements in air:",fluid_nelem1)

if (flag_write_gmsh_results==1) and (rank==0):
    msh2.mshWriter(
        cwd / (results_file + "_air_cavity_Mesh1.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1},
    )
    msh2.mshWriter(
        cwd / (results_file + "Mesh_control_volume.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements5},
    )

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf, JJf, Vffk, Vffm = acousticsFEM.getMatrices(
    fluid_nodes1, fluid_elements1, [celerity, rho]
)
KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofF=list(range(fluid_ndof1))

##############################################################
# Load structure mesh
##############################################################
# start reading mesh
mesh = msh.mshReader(cwd / (mesh_file + "_struc.msh"))

struc_nodes = mesh.getNodes()[:, 0:3]
struc_elements = mesh.getElements(tag=6)["TRI3"]
struc_boun = mesh.getElements(tag=7)["LIN2"]

struc_nnodes = struc_nodes.shape[0]
struc_nelem = struc_elements.shape[0]
struc_ndof = struc_nnodes * 6

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_struc_mesh_surface_6.msh"),
        struc_nodes,
        {"type": "TRI3", "connectivity": struc_elements},
    )
    msh2.mshWriter(
        cwd / (results_file + "_struc_boun_mesh_line_7.msh"),
        struc_nodes,
        {"type": "LIN2", "connectivity": struc_boun},
    )

if rank==0:
    logger.info("nnodes for structure=",struc_nnodes)
    logger.info("nelem for structure=",struc_nelem)

##################################################################
# compute level set
##################################################################

tic = time.process_time()

LevelSet, distance = acousticsXFEM.getLevelSet(
    fluid_nodes1, struc_nodes, struc_elements
)

toc = time.process_time()
if rank==0:
    logger.info("time to compute level set: {}".format(toc-tic))

tic = time.process_time()
tangent_nodes, tangent_mesh = acousticsXFEM.getTangentMesh(
    struc_nodes, struc_elements, struc_boun
)
LevelSetTangent, _ = acousticsXFEM.getLevelSet(
    fluid_nodes1, tangent_nodes, tangent_mesh
)
toc = time.process_time()
if rank==0:
    logger.info("time to compute tangent level set: {}".format(toc-tic))

if (flag_write_gmsh_results==1) and (rank==0):
    msh2.mshWriter(
        cwd / (results_file + "_level_sets.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1},
        fields=[
            {"name": "LevelSet", "type": "nodal", "data": LevelSet},
            {"name": "TangentLevelSet", "type": "nodal", "data": LevelSetTangent},
            {"name": "distance", "type": "nodal", "data": distance},
        ],
        append=True,
    )
    msh2.mshWriter(
        cwd / (results_file + "_LS_tangent_mesh.msh"),
        tangent_nodes,
        {"type": "TRI6", "connectivity": tangent_mesh},
    )


##################################################################
# Get enriched nodes and elements
##################################################################
## Caution many TET4 procedures are used
tic = time.process_time()
# we take only the 4 corners nodes : it should be ok if tetra has straight lines
LSEnrichedElements = acousticsXFEM.getEnrichedElements(
    fluidElements=fluid_elements1, levelset=LevelSet
)
# LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements1[:,list(range(4))],LevelSet) 

EnrichedElements = acousticsXFEM.getSurfEnrichedElements(
    fluid_nodes1, 
    fluid_elements1[LSEnrichedElements], 
    struc_nodes, 
    struc_elements
)
EnrichedElements = LSEnrichedElements[EnrichedElements]


toc = time.process_time()
if rank==0:
    logger.info("time to find surface enriched elements: {}".format(toc-tic))

tic = time.process_time()
EdgeEnrichedElements = acousticsXFEM.getEdgeEnrichedElements(
    fluid_nodes1,
    fluid_elements1,
    struc_nodes,
    struc_boun
)
#
EdgeEnrichedElementsInAllMesh = acousticsXFEM.getEnrichedElements(
    fluidElements=fluid_elements1,
    levelset=LevelSetTangent
    )

toc = time.process_time()
if rank==0:
    logger.info("time to find edge enriched elements: {}".format(toc-tic))

HeavisideEnrichedElements=np.setdiff1d(EnrichedElements,EdgeEnrichedElements)

if (flag_write_gmsh_results==1) and (rank==0):
    msh2.mshWriter(
        cwd / (results_file + "_LSenriched_elements.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1[LSEnrichedElements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_enriched_elements.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1[EnrichedElements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_edge_enriched_elements.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1[EdgeEnrichedElements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_edge_enriched_elements_in_all_mesh.msh"),
        fluid_nodes1,
        {
            "type": "TET10",
            "connectivity": fluid_elements1[EdgeEnrichedElementsInAllMesh],
        },
    )
    msh2.mshWriter(
        cwd / (results_file + "_heaviside_enriched_elements.msh"),
        fluid_nodes1,
        {"type": "TET10", "connectivity": fluid_elements1[HeavisideEnrichedElements]},
    )

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.process_time()

Enrichednodes = np.unique(fluid_elements1[EnrichedElements])

II, JJ, Vaak, Vaam, Vafk, Vafm = acousticsXFEM.getMatrices(
    fluid_nodes1, fluid_elements1[EnrichedElements], LevelSet, LevelSetTangent, [celerity, rho]
)


KAA = scipy.sparse.csc_matrix( (Vaak,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAA = scipy.sparse.csc_matrix( (Vaam,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
KAF = scipy.sparse.csc_matrix( (Vafk,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )
MAF = scipy.sparse.csc_matrix( (Vafm,(II,JJ)), shape=(fluid_ndof1,fluid_ndof1) )

SolvedDofA=Enrichednodes-1

toc = time.process_time()
if rank==0:
    logger.info("time to compute enrichment: {}".format(toc-tic))

##################################################################
# Construct the whole system
#################################################################

K=scipy.sparse.bmat( [
            [KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
            [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]] )

M=scipy.sparse.bmat( [
            [MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
            [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]] )

##################################################################
# Build Second member
##################################################################

# To impose the load on the fluid:
# fluid node number 1
UF = np.zeros(2*fluid_ndof1,dtype=float)
UF[1-1]=3.1250E-05

SolvedDof = np.hstack([SolvedDofF,SolvedDofA+fluid_ndof1])

##############################################################
# FRF computation
##############################################################

Flag_frf_analysis=1
frequencies=[]
frf=[]
frfgradient=[]

if (Flag_frf_analysis==1):
    logger.info("Proc. {} / time at the beginning of the FRF: {}".format(rank,time.ctime()))

    if rank==0:
        logger.info('nb of total dofs: ',len(SolvedDofF)+len(SolvedDofA))

    press_save=[]
    disp_save=[]
    
    #for i in range(nb_freq_step):
    for freq in np.linspace(freq_ini,freq_end,nb_freq_step):

        #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*np.pi*freq
        #print('omega = ',omega)

        logger.info("proc number {} - frequency= {}".format(rank,freq))

        tic = time.process_time()        
        
        F=np.array(omega**2*UF[SolvedDof] , dtype='float')

        #sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float')  , F , comm=mycomm  )
        sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=float) , np.array(F , dtype=float) )

        press1 = np.zeros((fluid_ndof1),dtype=complex)
        press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        enrichment=np.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
        CorrectedPressure=np.zeros((fluid_ndof1),dtype=complex)
        CorrectedPressure[SolvedDofA]=press1[SolvedDofA]+enrichment[SolvedDofA]*np.sign(LevelSet[SolvedDofA])
        frf.append(
            acousticsFEM.getQuadraticPressure(
                fluid_nodes1,
                fluid_elements5,
                press1
            )
        )
            
            
        #frf.append(scipy.dot(scipy.dot(M,sol),sol))

        #print(silex_lib_xfem_acou_tet4.makexfemposfile.__doc__)
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,LevelSet,press1.real,enrichment.real,'press_plus.pos')
        #silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,-LevelSet,press1.real,-enrichment.real,'press_moins.pos')

        if (flag_write_gmsh_results==1) and (rank==0):
            press_save.append(CorrectedPressure.real)        

    frfsave=[np.array(frequencies),np.array(frf)]

    logger.info("Proc. {} / time at the end of the FRF: {}".format(rank, time.ctime()))

    if (flag_write_gmsh_results==1) and (rank==0):
        msh2.mshWriter(
            cwd / (results_file + str(rank) + "_results_fluid_frf.msh"),
            fluid_nodes1,
            {"type": "TET10", "connectivity": fluid_elements1},
            fields=[
                {
                    "data": press_save,
                    "type": "nodal",
                    "nbsteps": nb_freq_step,
                    "name": "Pressure",
                }
            ],
        )
        

    with open(cwd / (results_file + "_results.frf"), "wb") as f:
        pickle.dump(frfsave, f)



