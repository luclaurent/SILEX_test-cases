import string
import time
from pathlib import Path
import numpy as np
from loguru import logger
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle


# import mumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
objFEM = silex_lib_fem.LinearAcousticsTRI3()
objXFEM = silex_lib_xfem.LinearAcousticsTRI3()

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
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file = "results/xfem"
results_file = "results/xfem_1"

##############################################################
# Material, Boundary conditions
##############################################################

# air
celerity = 340.0
rho = 1.2

freq_ini = 100.0

flag_write_gmsh_results = 1

flag_edge_enrichment = 0
# flag_edge_enrichment=1

# freq_comparaison = 210.0


##############################################################
# Load fluid mesh
##############################################################
tic = time.process_time()

cwd = Path(__file__).resolve().parent

# start reading mesh
mesh = msh.mshReader(cwd / (mesh_file + "_fluid.msh"))

fluid_nodes = mesh.getNodes()[:, 0:2]
fluid_elements = mesh.getElements(tag=1)["TRI3"]
Idnodes = np.unique(fluid_elements.flatten())

fluid_nnodes = fluid_nodes.shape[0]
fluid_nelem = fluid_elements.shape[0]
fluid_ndof = len(np.unique(fluid_elements.flatten()))

fluid_elements_boun = mesh.getElements(tag=2)["LIN2"]
IdnodeS2 = np.unique(fluid_elements_boun.flatten())

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_fluid_mesh.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements},
    )
    msh2.mshWriter(
        cwd / (results_file + "_fluid_boundary.msh"),
        fluid_nodes,
        {"type": "LIN2", "connectivity": fluid_elements_boun},
    )
logger.info("nnodes for fluid= {}".format(fluid_nnodes))
logger.info("nelem for fluid=", fluid_nelem)


##################################################################
# compute level set
##################################################################

tic = time.process_time()

x_pos_struc = 0.6200
h_struc = 0.65

LevelSet = fluid_nodes[:, 0] - x_pos_struc
LevelSetTangent = fluid_nodes[:, 1] - h_struc

# level set gradient with respect to parameters
LevelSet_gradient = -np.ones(fluid_nnodes)

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_level_sets.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements},
        fields=[
            {
                "data": LevelSet,
                "type": "nodal",
                "dim": 1,
                "name": "levelset",
            },
            {
                "data": LevelSetTangent,
                "type": "nodal",
                "dim": 1,
                "name": "levelset_tangent",
            },
        ],
        append=True,
    )
toc = time.process_time()
logger.info("time to compute level set: {}".format(toc - tic))


##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.process_time()

struc_nodes = np.array(
    [[x_pos_struc, 0.0], [x_pos_struc, h_struc / 2.0], [x_pos_struc, h_struc]]
)
struc_elements = np.array([[1, 2], [2, 3]])
struc_boun = np.array([3])

msh2.mshWriter(
    cwd / (results_file + "_struc_mesh.msh"),
    struc_nodes,
    {"type": "LIN2", "connectivity": struc_elements},
)

EnrichedElements = objXFEM.getEnrichedElements(
    fluid_nodes, fluid_elements, struc_nodes, struc_elements
)
toc = time.process_time()
logger.info("time to find surface enriched elements: {}".format(toc - tic))

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_enriched_elements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[EnrichedElements]},
    )
tic = time.process_time()

EdgeEnrichedElements = objXFEM.getEdgeEnrichedElements(
    fluid_nodes, fluid_elements, struc_nodes, struc_boun
)

toc = time.process_time()
logger.info("time to find edge enriched elements: {}".format(toc - tic))

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_edge_enriched_elements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[EdgeEnrichedElements]},
    )
##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf, JJf, Vffk, Vffm = objFEM.getMatrices(fluid_nodes, fluid_elements, [celerity, rho])

KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

SolvedDofF = np.setdiff1d(list(range(fluid_ndof)), IdnodeS2 - 1)
# SolvedDofF=range(fluid_ndof)

toc = time.process_time()
logger.info("time to compute fluid matrices: {}".format(toc - tic))

##################################################################
# Compute enrichment: Heaviside + Edge
##################################################################
tic = time.process_time()


NegativeLSelements, PositiveLSelements, NegativeLStgtElements, PositiveLStgtElements = (
    objXFEM.getLocationElementsLS(fluid_elements, LevelSet, LevelSetTangent)
)

EdgeEnrichedElementsInAllMesh = objXFEM.getEnrichedElements(
    fluidElements=fluid_elements, levelset=LevelSetTangent
)


IdElementTip = objXFEM.getElementContainingPoint(
    fluid_elements, fluid_nodes, [0.6, 0.65]
)
# .getelementcontainingpoint(fluid_elements,fluid_nodes,[0.6,0.65])

if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_NegativeLSelements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[NegativeLSelements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_PositiveLSelements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[PositiveLSelements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_NegativeLStgtElements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[NegativeLStgtElements]},
    )
    msh2.mshWriter(
        cwd / (results_file + "_PositiveLStgtElements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[PositiveLStgtElements]},
    )


IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = objXFEM.getMatrices(
    fluid_nodes,
    fluid_elements,
    LevelSet,
    LevelSetTangent,
    [celerity, rho],
    flag_edge_enrichment,
)

KAA = scipy.sparse.csc_matrix((Vaak, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
MAA = scipy.sparse.csc_matrix((Vaam, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
KAF = scipy.sparse.csc_matrix((Vafk, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))
MAF = scipy.sparse.csc_matrix((Vafm, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))

toc = time.process_time()
logger.info("time to compute Heaviside enrichment: {}".format(toc - tic))

# Enrichednodes = np.unique(fluid_elements[scipy.hstack(([HeavisideEnrichedElements,EdgeEnrichedElements]))])
# Enrichednodes = np.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements,EdgeEnrichedElementsInAllMesh]))])
# Enrichednodes = np.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements]))])
# Enrichednodes = np.unique(fluid_elements[scipy.hstack(([NegativeLStgtElements]))])
Enrichednodes = np.unique(fluid_elements[EnrichedElements])
# Enrichednodes = np.unique(fluid_elements)
SolvedDofA = Enrichednodes - 1

msh2.mshWriter(
    cwd / (results_file + "_EnrichedElements.msh"),
    fluid_nodes,
    {"type": "TRI3", "connectivity": fluid_elements[EnrichedElements.flatten()]},
)
#################################################################
# Construct the whole system
##################################################################

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

#################################################################
# Compute gradients with respect to parameters
##################################################################

IIf, JJf, Vfak_gradient, Vfam_gradient, _ = objXFEM.getGradientMatrices(
    fluid_nodes,
    fluid_elements[EnrichedElements],
    levelset=LevelSet,
    levelsetTg=None,
    levelsetGradient=LevelSet_gradient,
    material=[celerity, rho],
)


dMFA_dtheta = scipy.sparse.csc_matrix(
    (Vfam_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
)
dKFA_dtheta = scipy.sparse.csc_matrix(
    (Vfak_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
)

with open(cwd / (results_file + "_KAF_MAF_6200.pck"), "wb") as f:
    pickle.dump(
        [
            dKFA_dtheta[SolvedDofF, :][:, SolvedDofA],
            dMFA_dtheta[SolvedDofF, :][:, SolvedDofA],
            KAF[SolvedDofF, :][:, SolvedDofA],
            MAF[SolvedDofF, :][:, SolvedDofA],
        ],
        f,
    )


##############################################################
# FRF computation of the FSI problem
##############################################################

FF = np.zeros(fluid_ndof)

freq = freq_ini
omega = 2 * np.pi * freq

FF[SolvedDofF] = -(
    KFF[SolvedDofF, :][:, IdnodeS2 - 1]
    - (omega * omega) * MFF[SolvedDofF, :][:, IdnodeS2 - 1]
) * (np.zeros((len(IdnodeS2))) + 1.0)
FA = np.zeros(fluid_ndof)
F = FF[SolvedDofF]
F = np.concatenate((F, FA[SolvedDofA]))
F = np.array(F)

sol = scipy.sparse.linalg.spsolve(K - (omega * omega) * M, F)

press = np.zeros(fluid_ndof)
press[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
enrichment = np.zeros(fluid_nnodes)
enrichment[SolvedDofA] = sol[
    list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
]
CorrectedPressure = press
CorrectedPressure[np.ix_(SolvedDofA)] = CorrectedPressure[SolvedDofA] + enrichment[
    SolvedDofA
] * np.sign(LevelSet[SolvedDofA])
# quadratic_pressure=silex_lib_tri3_acou.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure)
quadratic_pressure = objXFEM.getQuadraticPressure(
    fluid_nodes,
    fluid_elements,
    press,
    enrichment,
    LevelSet,
    LevelSetTangent,
    flag_edge_enrichment,
)


if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
        cwd / (results_file + "_results_fluid_frf.msh"),
        fluid_nodes,
        [{"type": "TRI3", "connectivity": fluid_elements}],
        fields={
            "data": CorrectedPressure,
            "type": "nodal",
            "dim": 1,
            "name": "pressure",
        },
        append=True,
    )
