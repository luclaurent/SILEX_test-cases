import string
import time
from pathlib import Path
import numpy as np
from loguru import logger
import scipy
import scipy.sparse
import scipy.sparse.linalg
import gmsh

import pylab as pl
import pickle

# import mumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem,silex_lib_tri3
from SILEXlib import MeshField

# load classes
objFEM = silex_lib_fem.LinearAcousticsTRI3()
objXFEM = silex_lib_xfem.LinearAcousticsTRI3()

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

freq_ini = 10.0
freq_end = 1000.0
nbfreq = 1000

flag_write_gmsh_results = 1

flag_edge_enrichment = 0
flag_edge_enrichment=1

dirichlet = False

# freq_comparaison = 210.0


# remesh
cwd = Path(__file__).resolve().parent
if rank == 0:
    gmsh.initialize()
    gmsh.open(str(cwd / (mesh_file + "_fluid.geo")))
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.refine()
    gmsh.model.mesh.refine()
    gmsh.model.mesh.refine()
    gmsh.model.mesh.refine()
    gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    # gmsh.model.mesh.refine()
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)  
    gmsh.write(str(cwd / (mesh_file + "_fluid.msh")))

comm.Barrier()

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

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

LevelSet = fluid_nodes[:, 0] - 0.75
LevelSetTangent = fluid_nodes[:, 1] - 1.5

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

struc_nodes = np.array([[0.75, 0.0], [0.75, 0.3], [0.75, 1.5]])
struc_elements = np.array([[1, 2], [2, 3]])
struc_boun = np.array([3])

if rank ==0:
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

if dirichlet:
    SolvedDofF = np.setdiff1d(list(range(fluid_ndof)), IdnodeS2 - 1)
else:
    SolvedDofF = range(fluid_ndof)

toc = time.process_time()
logger.info("time to compute fluid matrices: {}".format(toc - tic))

##################################################################
# Compute enrichment: Heaviside + Edge
##################################################################
tic = time.process_time()

# HeavisideEnrichedElements=np.setdiff1d(EnrichedElements,EdgeEnrichedElements)

# Enrichednodes = np.unique(fluid_elements[HeavisideEnrichedElements])
# Enrichednodes = np.unique(fluid_elements[EnrichedElements])

# print xvibacoufo.getpositivenegativeelts.__doc__

NegativeLSelements, PositiveLSelements, NegativeLStgtElements, PositiveLStgtElements = (
    objXFEM.getLocationElementsLS(fluid_elements, LevelSet, LevelSetTangent)
)

# NegativeLSelements=NegativeLSelements[list(range(nbNegLS))]
# PositiveLSelements=PositiveLSelements[list(range(nbPosLS))]
# NegativeLStgtElements=NegativeLStgtElements[list(range(nbNegLSt))]
# PositiveLStgtElements=PositiveLStgtElements[list(range(nbPosLSt))]

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

# Enrichednodes = np.unique(fluid_elements[np.hstack(([HeavisideEnrichedElements,EdgeEnrichedElements]))])
# Enrichednodes = np.unique(fluid_elements[np.hstack(([EnrichedElements,PositiveLStgtElements,EdgeEnrichedElementsInAllMesh]))])
# Enrichednodes = np.unique(fluid_elements[np.hstack(([EnrichedElements,PositiveLStgtElements]))])
# Enrichednodes = np.unique(fluid_elements[np.hstack(([NegativeLStgtElements]))])
Enrichednodes = np.unique(fluid_elements[EnrichedElements])
# Enrichednodes = np.unique(fluid_elements)
SolvedDofA = Enrichednodes - 1

if rank ==0:
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

##############################################################
# FRF computation of the FSI problem
##############################################################

FF = np.zeros(fluid_ndof)
solvFRF = list()
solvPRESS = list()
solvENRICH = list()

freqlist = np.linspace(freq_ini, freq_end, nbfreq)

chunks = [[] for _ in range(nproc)]
for i, chunk in enumerate(freqlist):
    chunks[i % nproc].append(chunk)
# Divide the data among processes
freqlistMPI = comm.scatter(chunks, root=0)

for freq in freqlistMPI:
    omega = 2 * np.pi * freq

    if dirichlet:    
        press[IdnodeS2-1] = np.ones(len(IdnodeS2))
        FF[SolvedDofF] = -(
        KFF[np.ix_(SolvedDofF,IdnodeS2-1)]
        - (omega * omega) * MFF[np.ix_(SolvedDofF,IdnodeS2-1)]
        ) * press[IdnodeS2-1]
    else:
        UF = 1e-3
        forceSurf = np.array([UF,0,UF,0 ])*omega**2
        ptsbnd = [fluid_nodes[IdnodeS2-1,0].min(),
                fluid_nodes[IdnodeS2-1,1].min(),
                fluid_nodes[IdnodeS2-1,0].max(),
                fluid_nodes[IdnodeS2-1,1].max()]
        FF= silex_lib_tri3.forceonline(fluid_nodes,fluid_elements_boun,forceSurf,ptsbnd)
        FF = FF[::2]


    FA = np.zeros(fluid_ndof)
    F = FF[SolvedDofF]
    F = np.concatenate((F, FA[SolvedDofA]))

    sol = scipy.sparse.linalg.spsolve(K - (omega * omega) * M, F)

    press = np.zeros(fluid_ndof)
    press[IdnodeS2-1] = np.ones(len(IdnodeS2))
    press[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
    enrichment = np.zeros(fluid_nnodes)
    enrichment[SolvedDofA] = sol[
        list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
    ]
    CorrectedPressure = press.copy()
    CorrectedPressure[np.ix_(SolvedDofA)] = CorrectedPressure[SolvedDofA] + \
                                            enrichment[SolvedDofA] * np.sign(LevelSet[SolvedDofA])
    solvPRESS.append(press)
    solvENRICH.append(enrichment)
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
    solvFRF.append(quadratic_pressure)

comm.Barrier()
# Send the results back to the master processes
newData = comm.gather(solvFRF,root=0)
newDataPRESS = comm.gather(solvPRESS,root=0)
newDataENRICH = comm.gather(solvENRICH,root=0)

comm.Barrier()



if rank ==0:
    freqListFinal = np.concatenate(chunks)
    IXsort = np.argsort(freqListFinal)
    freqListFinal = freqListFinal[IXsort]
    quadpressure = np.concatenate(newData)
    quadpressure = quadpressure[IXsort] 
    solpressure = np.concatenate(newDataPRESS)
    solpressure = solpressure[IXsort]
    solenrich = np.concatenate(newDataENRICH)
    solenrich = solenrich[IXsort]

    data ={"freq": freqListFinal, "mag": quadpressure, "nbnodes": fluid_nnodes, "nelem": fluid_nelem}
    with open("solvFRF"+str(fluid_nnodes).zfill(6)+".pkl", "wb") as f:
        pickle.dump(data, f)


    # prepare final data
    objMesh = MeshField.MeshField(fluid_nodes, fluid_elements, LevelSet, LevelSetTangent)
    objMesh.addField(np.array(solpressure).transpose(),np.array(solenrich).transpose())
    # extract final data
    data = objMesh.getData(['nodes','mesh','levelset','levelset_tangent','fields'])
    msh2.mshWriter(
            cwd / (results_file + "_results_fluid_frf_along.msh"),
            data['nodes'],
            [{"type": "TRI3", "connectivity": data['TRI3']},
            {"type": "QUA4", "connectivity": data['QUA4']}],
            fields={
                "data": data['fields'],
                "type": "nodal",
                "nbsteps": nbfreq,
                "name": "pressure",
            },
            append=True,
        )

    # prepare final data
    objMesh = MeshField.MeshField(fluid_nodes, fluid_elements, LevelSet, LevelSetTangent)
    objMesh.addField(press,enrichment)
    # extract final data
    data = objMesh.getData(['nodes','mesh','levelset','levelset_tangent','fields'])


    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_results_fluid_frf.msh"),
            data['nodes'],
            [{"type": "TRI3", "connectivity": data['TRI3']},
            {"type": "QUA4", "connectivity": data['QUA4']}],
            fields={
                "data": data['fields'],
                "type": "nodal",
                "dim": 1,
                "name": "pressure",
            },
            append=True,
        )
        msh2.mshWriter(
        cwd / (results_file + "_results_fluid_frf_raw.msh"),
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

    print('')