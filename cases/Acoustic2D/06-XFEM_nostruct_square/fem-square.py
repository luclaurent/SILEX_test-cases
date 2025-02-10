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
from SILEXlib import silex_lib_fem
from SILEXlib import MeshField

# load classes
objFEM = silex_lib_fem.LinearAcousticsTRI3()

from mpi4py import MPI

comm = MPI.COMM_WORLD

nproc = comm.Get_size()
rank = comm.Get_rank()

# mpirun -np 2 python Main_xfem.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file = "xfem"
results_file = "xfem_1"

##############################################################
# Material, Boundary conditions
##############################################################

# air
celerity = 340.0
rho = 1.2

freq_ini = 500.0

flag_write_gmsh_results = 1

flag_edge_enrichment = 0
# flag_edge_enrichment=1

# freq_comparaison = 210.0

dirichlet = False


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
        cwd / (results_file + "_fem_fluid_boundary.msh"),
        fluid_nodes,
        {"type": "LIN2", "connectivity": fluid_elements_boun},
    )
logger.info("nnodes for fluid= {}".format(fluid_nnodes))
logger.info("nelem for fluid=", fluid_nelem)


##############################################################
# Compute Standard Fluid Matrices
##############################################################
tic = time.process_time()

IIf, JJf, Vffk, Vffm = objFEM.getMatrices(fluid_nodes, fluid_elements, [celerity, rho])

KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

if dirichlet:
    SolvedDofF = np.setdiff1d(list(range(fluid_ndof)), IdnodeS2 - 1)
else:
    SolvedDofF=range(fluid_ndof)

toc = time.process_time()
logger.info("time to compute fluid matrices: {}".format(toc - tic))

#################################################################
# Construct the whole system
##################################################################
K = KFF[SolvedDofF, :][:, SolvedDofF]
M = MFF[SolvedDofF, :][:, SolvedDofF]

##############################################################
# FRF computation of the FSI problem
##############################################################

FF = np.zeros(fluid_ndof)
press = np.zeros(fluid_ndof)

freq = freq_ini
omega = 2 * np.pi * freq

if dirichlet:
    press[IdnodeS2-1] = np.ones(len(IdnodeS2))
    FF[SolvedDofF]=-(KFF[np.ix_(SolvedDofF,IdnodeS2-1)]-
                 (omega*omega)*MFF[np.ix_(SolvedDofF,IdnodeS2-1)])@press[IdnodeS2-1]
else:
    FF[IdnodeS2 - 1] = 1e-3

F = FF[SolvedDofF]

sol = scipy.sparse.linalg.spsolve(K - (omega * omega) * M, F)

press[SolvedDofF] = sol[list(range(len(SolvedDofF)))]

# quadratic_pressure=silex_lib_tri3_acou.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure)
quadratic_pressure = objFEM.getQuadraticPressure(
    fluid_nodes,
    fluid_elements,
    press
)


if (flag_write_gmsh_results == 1) and (rank == 0):
    msh2.mshWriter(
    cwd / (results_file + "_results_fem_fluid_frf_raw.msh"),
    fluid_nodes,
    [{"type": "TRI3", "connectivity": fluid_elements}],
    fields={
        "data": press,
        "type": "nodal",
        "dim": 1,
        "name": "pressure",
    },
    append=True,
    )

print('')