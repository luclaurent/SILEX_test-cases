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


import pymumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
acousticsFEM = silex_lib_fem.LinearAcousticsTET4()
acousticsXFEM = silex_lib_xfem.LinearAcousticsTET4()
# structureFEM = silex_lib_fem.DKT()


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
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7-v4.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

if rank == 0:
    logger.info("time at the beginning of the computation: {}".format(time.ctime()))

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

mesh_file = "geom/cavity_acou3D_struc_3D_v3"
results_file = "results/cavity_acou3D_struc_3D_v3"
cwd = Path(__file__).resolve().parent

flag_write_gmsh_results = 1

freq_ini = 90.0
freq_end = 120.0
nb_freq_step_per_proc = 2

nb_freq_step = nb_freq_step_per_proc * nproc
deltafreq = (freq_end - freq_ini) / (nb_freq_step - 1)

# air
celerity = 343.0  # ok
rho = 1.21  # ok

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

mesh = msh.mshReader(cwd / (mesh_file + "_air.msh"))

fluid_nodes = mesh.getNodes()
fluid_elements1 = mesh.getElements(tag=1)["TET4"]  # air, cavity
fluid_elements5 = mesh.getElements(tag=5)["TET4"]  # air, control volum

IdNodes1 = np.unique(fluid_elements1.flatten())
IdNodes5 = np.unique(fluid_elements5.flatten())

fluid_nnodes = fluid_nodes.shape[0]
fluid_nelem1 = fluid_elements1.shape[0]
fluid_nelem5 = fluid_elements5.shape[0]

fluid_nnodes1 = IdNodes1.shape[0]
fluid_nnodes5 = IdNodes5.shape[0]

fluid_ndof = fluid_nnodes

if rank == 0:
    logger.info("Number of nodes: {}".format(fluid_nnodes))
    logger.info("Number of elements in air: {}".format(fluid_nelem1))
    logger.info("Number of nodes in air: {}".format(fluid_nnodes1))

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_air_cavity_Mesh1.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1},
    )
    msh2.mshWriter(
        cwd / (results_file + "_air_controlled_volume_Mesh5.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements5},
    )

##############################################################
# Load structure mesh
##############################################################

# start reading mesh
mesh = msh.mshReader(cwd / (mesh_file + "_struc.msh"))

struc_nodes = mesh.getNodes()
struc_elements = mesh.getElements(tag=2)["TRI3"]

Idnodes_S_air_interface = np.unique(struc_elements.flatten())

struc_nnodes = struc_nodes.shape[0]
struc_nelem = struc_elements.shape[0]
struc_ndof = struc_nnodes * 6

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_struc_surface.msh"),
        struc_nodes,
        {"type": "TRI3", "connectivity": struc_elements},
    )

if rank == 0:
    logger.info("nnodes for structure= {}".format(struc_nnodes))
    logger.info("nelem for structure= {}".format(struc_nelem))


##################################################################
# compute level set and its gradient according to a parameter
##################################################################

tic = time.process_time()

# LS from a structure mesh
LevelSet_from_Mesh, distance = acousticsXFEM.getLevelSet(
    fluid_nodes, struc_nodes, struc_elements
)

# LS from a simple analytic shape (sphere)
lx3 = 2.0  # Xc
ly3 = 2.0  # Yc
R = 1.0  # sphere radius

LevelSet = (
    np.sqrt(
        (fluid_nodes[:, 0] - lx3) ** 2
        + (fluid_nodes[:, 1] - lx3) ** 2
        + (fluid_nodes[:, 2] - 0.0) ** 2
    )
    - R
)
# Compute LS gradient according to Xc
LevelSet_gradient = (lx3 - fluid_nodes[:, 0]) / (
    np.sqrt(
        (fluid_nodes[:, 0] - lx3) ** 2
        + (fluid_nodes[:, 1] - lx3) ** 2
        + (fluid_nodes[:, 2] - 0.0) ** 2
    )
)


toc = time.process_time()
if rank == 0:
    logger.info("time to compute level set: {}".format(toc - tic))

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_level_sets.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1},
        fields=[
            {"name": "LevelSet (from analytics", "type": "nodal", "data": LevelSet},
            {
                "name": "TangentLevelSet (from analytics)",
                "type": "nodal",
                "data": LevelSet_gradient,
            },
            {
                "name": "LevelSet (from mesh)",
                "type": "nodal",
                "data": LevelSet_from_Mesh,
            },
            {"name": "distance (from mesh)", "type": "nodal", "data": distance},
        ],
        append=True,
    )
    msh2.mshWriter(
        cwd / (results_file + "_struc_air_interface.msh"),
        struc_nodes,
        {"type": "TRI3", "connectivity": struc_elements},
    )

##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.process_time()

LSEnrichedElements = acousticsXFEM.getEnrichedElements(
    fluidElements=fluid_elements1, levelset=LevelSet
)

msh2.mshWriter(
    cwd / (results_file + "_LSenriched_elements.msh"),
    fluid_nodes,
    {"type": "TET4", "connectivity": fluid_elements1[LSEnrichedElements]},
)

##EnrichedElements=LSEnrichedElements#[EnrichedElements-1]
LSEnrichednodes = np.unique(fluid_elements1[LSEnrichedElements])

tmp = []
for i in LSEnrichednodes:
    for j in range(4):
        tmpp = np.where(fluid_elements1[:, j] == i)[0]
        for k in range(len(tmpp)):
            tmp.append(tmpp[k])
##    tmp.append(np.where(fluid_elements1[:,1]==i))
##    tmp.append(np.where(fluid_elements1[:,2]==i))
##    tmp.append(np.where(fluid_elements1[:,3]==i))

tmp = np.unique(np.array(tmp))
##tmp1,elttest0,tmp2=scipy.intersect1d(fluid_elements1[:,0],LSEnrichednodes,return_indices=True)
# silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements_test0',fluid_nodes,fluid_elements1[tmp],4)
# [75804, 97252, 97253,34973, 93135, 93137, 93248,83787, 93136,93525]
EnrichedElements0 = acousticsXFEM.getSurfEnrichedElements(
    fluid_nodes, fluid_elements1[tmp], struc_nodes, struc_elements
)
EnrichedElements = tmp[EnrichedElements0]

toc = time.process_time()
if rank == 0:
    logger.info("time to find enriched elements: {}".format(toc - tic))

tic = time.process_time()


LS_moins_enriched = np.setdiff1d(LSEnrichedElements, EnrichedElements)
enriched_moins_LS = np.setdiff1d(EnrichedElements, LSEnrichedElements)

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_enriched_elements.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1[EnrichedElements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_LS_moins_enriched.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1[LS_moins_enriched]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_enriched_moins_LS.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1[enriched_moins_LS]},
    )

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf, JJf, Vffk, Vffm = acousticsFEM.getMatrices(
    fluid_nodes, fluid_elements1, [celerity, rho]
)

KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

SolvedDofF = list(range(fluid_ndof))


##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.process_time()

Enrichednodes = np.unique(fluid_elements1[EnrichedElements])

IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = acousticsXFEM.getMatrices(
    fluid_nodes, fluid_elements1[EnrichedElements], LevelSet, [celerity, rho]
)

KAA = scipy.sparse.csc_matrix((Vaak, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
MAA = scipy.sparse.csc_matrix((Vaam, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
KAF = scipy.sparse.csc_matrix((Vafk, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))
MAF = scipy.sparse.csc_matrix((Vafm, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))

SolvedDofA = Enrichednodes - 1

toc = time.process_time()
if rank == 0:
    logger.info("time to compute Heaviside enrichment: {}".format(toc - tic))

##################################################################
# Construct the whole system
#################################################################

K = scipy.sparse.bmat(
    [
        [KFF[SolvedDofF, :][:, SolvedDofF], KAF[SolvedDofF, :][:, SolvedDofA]],
        [KAF[SolvedDofA, :][:, SolvedDofF], KAA[SolvedDofA, :][:, SolvedDofA]],
    ]
)

M = scipy.sparse.bmat(
    [
        [MFF[SolvedDofF, :][:, SolvedDofF], MAF[SolvedDofF, :][:, SolvedDofA]],
        [MAF[SolvedDofA, :][:, SolvedDofF], MAA[SolvedDofA, :][:, SolvedDofA]],
    ]
)

##################################################################
# Build Second member
##################################################################

# To impose the load on the fluid:
# fluid node number 1
UF = np.zeros(2 * fluid_ndof, dtype=float)
UF[9 - 1] = 3.1250e-05

SolvedDof = np.hstack([SolvedDofF, SolvedDofA + fluid_ndof])

#################################################################
# Compute gradients with respect to parameters
##################################################################


# logger.info(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
IIf, JJf, Vfak_gradient, Vfam_gradient = acousticsXFEM.getGradientMatrices(
    nodes=fluid_nodes,
    elements=fluid_elements1,
    levelset=LevelSet,
    levelsetGradient=LevelSet_gradient,
    material=[celerity, rho],
)

dKFA_dtheta = scipy.sparse.csc_matrix(
    (Vfak_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
)
dMFA_dtheta = scipy.sparse.csc_matrix(
    (Vfam_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
)

dK = scipy.sparse.bmat(
    [
        [None, dKFA_dtheta[SolvedDofF, :][:, SolvedDofA]],
        [dKFA_dtheta[SolvedDofA, :][:, SolvedDofF], None],
    ]
)
dM = scipy.sparse.bmat(
    [
        [None, dMFA_dtheta[SolvedDofF, :][:, SolvedDofA]],
        [dMFA_dtheta[SolvedDofA, :][:, SolvedDofF], None],
    ]
)

##############################################################
# FRF computation
##############################################################

Flag_frf_analysis = 1
frequencies = []
frf = []
frfgradient = []

if Flag_frf_analysis == 1:
    logger.info(
        "Proc. {} / time at the beginning of the FRF: {}".format(rank, time.ctime())
    )

    if rank == 0:
        logger.info("nb of total dofs: ", len(SolvedDofF) + len(SolvedDofA))

    press_save = []
    disp_save = []
    dpress_save = []

    for i in range(nb_freq_step_per_proc):
        freq = freq_ini + i * nproc * deltafreq + rank * deltafreq
        frequencies.append(freq)
        omega = 2 * np.pi * freq

        logger.info("proc number {} - frequency={}".format(rank, freq))

        tic = time.process_time()

        dtype = complex

        F = np.array(omega**2 * UF[SolvedDof], dtype=dtype)

        sol = scipy.sparse.linalg.spsolve(K - (omega**2) * M, F)

        # sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , F )

        press1 = np.zeros((fluid_ndof), dtype=dtype)
        press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
        enrichment = np.zeros((fluid_nnodes), dtype=dtype)
        enrichment[SolvedDofA] = sol[
            list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
        ]
        CorrectedPressure = np.zeros((fluid_ndof), dtype=dtype)
        CorrectedPressure[SolvedDofA] = press1[SolvedDofA] + enrichment[
            SolvedDofA
        ] * np.sign(LevelSet[SolvedDofA])
        # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
        frf.append(
            acousticsXFEM.getQuadraticPressure(
                fluid_nodes,
                fluid_elements5,
                press1,
                enrichment,
                LevelSet,
                LevelSet * 0 - 1.0,
            )
        )

        # frf.append(scipy.dot(scipy.dot(M,sol),sol))

        # logger.info(silex_lib_xfem_acou_tet4.makexfemposfile.__doc__)
        # silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,LevelSet,press1.real,enrichment.real,'press_plus.pos')
        # silex_lib_xfem_acou_tet4.makexfemposfile(fluid_nodes,fluid_elements1,-LevelSet,press1.real,-enrichment.real,'press_moins.pos')

        if (flag_write_gmsh_results == 1) and (rank == 0):
            press_save.append(CorrectedPressure.real)

        tmp = -(dK - (omega**2) * dM) * sol
        Dsol_Dtheta = sol = scipy.sparse.linalg.spsolve(K - (omega**2) * M, tmp)
        # Dsol_Dtheta = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )

        Dpress_Dtheta = np.zeros(fluid_ndof, dtype=dtype)
        Dpress_Dtheta[SolvedDofF] = Dsol_Dtheta[list(range(len(SolvedDofF)))]
        Denrichment_Dtheta = np.zeros(fluid_ndof, dtype=dtype)
        Denrichment_Dtheta[SolvedDofA] = Dsol_Dtheta[
            list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
        ]
        # logger.info(silex_lib_xfem_acou_tet4.computegradientcomplexquadratiquepressure.__doc__)

        #####################
        #####################
        # store gradients
        frfgradient.append(
            acousticsXFEM.getGradientQuadraticPressure(
                nodes=fluid_nodes,
                elements=fluid_elements5,
                pressureField=press1,
                gradientPressureField=Dpress_Dtheta,
                levelset=LevelSet,
            )
        )
        dpress_save.append(Dpress_Dtheta.copy())

    frfsave = [np.array(frequencies), np.array(frf)]

    # comm.send(frfsave, dest=0, tag=11)

    logger.info("Proc. {} / time at the end of the FRF: {}".format(rank, time.ctime()))

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + str(rank) + "_results_fluid_frf.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1},
            fields=[
                {
                    "data": press_save,
                    "type": "nodal",
                    "nbsteps": nb_freq_step,
                    "name": "Pressure",
                },
                {
                    "data": dpress_save,
                    "type": "nodal",
                    "nbsteps": nb_freq_step,
                    "name": "Pressure gradient",
                },
            ],
            append=True,
        )

    # Save the FRF problem
    # Allfrequencies=np.zeros(nb_freq_step)
    # Allfrf=np.zeros(nb_freq_step)
    # k=0
    # if rank==0:
    #   for i in range(nproc):
    #        data = comm.recv(source=i, tag=11)
    #       for j in range(len(data[0])):
    #            Allfrequencies[k]=data[0][j]
    #            Allfrf[k]=data[1][j]
    #           k=k+1

    # Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
    # Allfrfsave=[np.array(list(Allfrequencies)),np.array(list(Allfrf))]
    with open(cwd / (results_file + "_results.frf"), "wb") as f:
        pickle.dump(frfsave, f)
