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


def computeFreqPerProc(nbStep, nbProc, freqInit, freqEnd):
    # compute integer number of freq per proc and remaining steps
    nbFreqProc = nbStep // nbProc
    nbFreqProcRemain = nbStep % nbProc
    # compute frequencies steps
    varCase = 1
    if nbFreqProcRemain == 0:
        varCase = 0
    listFreq = np.zeros((nbFreqProc + varCase, nbProc))
    listAllFreq = np.linspace(freqInit, freqEnd, nbStep)
    # logger.info(np.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc + varCase):
            if itC * nbProc + itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1

    # logger.info(listFreq)
    return listFreq


###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py

# mpirun -np 2 python Main_xfem_2.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################


def RunPb(
    freqMin,
    freqMax,
    nbStep,
    nbProc,
    rank,
    comm,
    positionStructX,
    positionStructY,
    radiusStruct,
    angleStruct,
):
    # parallepipedic cavity with plane structure
    mesh_file = "geom/xfem3_square"
    results_file_ini = "results_square/xfem_3_"
    cwd = Path(__file__).resolve().parent

    listFreqPerProc = computeFreqPerProc(nbStep, nbProc, freqMin, freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity = 340.0
    rho = 1.2
    fluid_damping = 1 + 0.01j

    nproc = comm.Get_size()
    rank = comm.Get_rank()

    flag_write_gmsh_results = 1

    flag_edge_enrichment = 0

    h_struc = 2.0

    ValpStructX = positionStructX * 10000.0
    ValpStructY = positionStructY * 10000.0
    ValpStructR = radiusStruct * 10000.0
    ValpStructT = angleStruct * 10000.0

    file_extension = (
        str(ValpStructX)[0:5]
        + "_"
        + str(ValpStructY)[0:5]
        + "_"
        + str(ValpStructR)[0:4]
        + "_"
        + str(ValpStructT)[0:4]
    )
    results_file = results_file_ini + file_extension
    logger.info(results_file)

    # center and radius of half-circle
    x_pos_struc = positionStructX
    y_pos_struc = positionStructY
    radius_hcircle = radiusStruct
    logger.info(x_pos_struc)
    logger.info(y_pos_struc)
    logger.info(radius_hcircle)
    logger.info(angleStruct)
    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

    # start reading mesh
    mesh = msh.mshReader(cwd / (mesh_file + "_fluid.msh"))

    fluid_nodes = mesh.getNodes()[:, 0:2]
    fluid_elements = mesh.getElements(tag=1)["TRI3"]
    fluid_elements5 = mesh.getElements(tag=5)["TRI3"]
    Idnodes = np.unique(fluid_elements.flatten())
    Idnodes5 = np.unique(fluid_elements5.flatten())

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
            cwd / (results_file + "_control_volume_fluid_mesh.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements5},
        )
        msh2.mshWriter(
            cwd / (results_file + "_fluid_boundary.msh"),
            fluid_nodes,
            {"type": "LIN2", "connectivity": fluid_elements_boun},
        )
    logger.info("nnodes for fluid= {}".format(fluid_nnodes))
    logger.info("nelem for fluid= {}".format(fluid_nelem))
    logger.info("nelem for control volume= {}".format(fluid_elements5.shape[0]))

    ##################################################################
    # compute level set
    ##################################################################

    tic = time.process_time()

    # create coordinates of nodes of the structure (half circle)
    nbNodesHC = 50
    nbNodesSC = 25
    nbNodesWall = 1
    #
    angleU = angleStruct * np.pi / 180
    thicknessU = 1.0
    # points
    xc2 = x_pos_struc - radius_hcircle * np.cos(angleU)
    yc2 = y_pos_struc - radius_hcircle * np.sin(angleU)
    xc4 = x_pos_struc + radius_hcircle * np.cos(angleU)
    yc4 = y_pos_struc + radius_hcircle * np.sin(angleU)
    # parameter for circles
    thetaHC = np.linspace(-angleU - np.pi / 2, -angleU + np.pi / 2, nbNodesHC)
    thetaSC = np.linspace(-angleU - np.pi / 2, -angleU + np.pi / 2, nbNodesSC)
    # inner large circle
    xNodesIHC = x_pos_struc - (radius_hcircle - thicknessU / 2.0) * np.sin(thetaHC)
    yNodesIHC = y_pos_struc - (radius_hcircle - thicknessU / 2.0) * np.cos(thetaHC)
    # outer larger circle
    xNodesOHC = x_pos_struc - (radius_hcircle + thicknessU / 2.0) * np.sin(
        thetaHC[::-1]
    )
    yNodesOHC = y_pos_struc - (radius_hcircle + thicknessU / 2.0) * np.cos(
        thetaHC[::-1]
    )
    # small circle 2
    xNodesC2 = xc2 - thicknessU / 2 * np.sin(thetaSC[::-1] + np.pi)
    yNodesC2 = yc2 - thicknessU / 2 * np.cos(thetaSC[::-1] + np.pi)
    # small circle 4
    xNodesC4 = xc4 + thicknessU / 2 * np.sin(thetaSC[::-1])
    yNodesC4 = yc4 + thicknessU / 2 * np.cos(thetaSC[::-1])

    struc_nodes = np.vstack(
        [
            np.hstack([xNodesIHC, xNodesC2[1:], xNodesOHC[1:], xNodesC4[1:]]),
            np.hstack([yNodesIHC, yNodesC2[1:], yNodesOHC[1:], yNodesC4[1:]]),
        ]
    ).transpose()
    logger.info(struc_nodes)

    nbNodesAllStruct = 2 * nbNodesHC + 2 * nbNodesSC - 4
    lCA = np.linspace(1, nbNodesAllStruct - 1, nbNodesAllStruct - 1)
    lCB = np.linspace(2, nbNodesAllStruct, nbNodesAllStruct - 1)

    struc_elements = np.vstack([lCA, lCB]).transpose()
    struc_elements = np.vstack([struc_elements, [1, nbNodesAllStruct]])
    logger.info(struc_elements)

    LevelSet = np.zeros(fluid_nnodes)
    LevelSet_gradient_X = np.zeros(fluid_nnodes)
    LevelSet_gradient_Y = np.zeros(fluid_nnodes)
    LevelSet_gradient_R = np.zeros(fluid_nnodes)
    LevelSet_gradient_T = np.zeros(fluid_nnodes)
    IndicZone = np.zeros(fluid_nnodes)
    # level set 3 cases
    for itN in range(fluid_nnodes):
        xNcurr = fluid_nodes[itN, 0]
        yNcurr = fluid_nodes[itN, 1]
        # check areas
        if angleU != np.pi / 2 or angleU != 3 * np.pi / 2:
            zoneA = yNcurr - np.tan(angleU) * (xNcurr - x_pos_struc) - y_pos_struc < 0
        else:
            zoneA = xNcurr - x_pos_struc < 0

        zoneAA = (
            np.sqrt((xNcurr - x_pos_struc) ** 2 + (yNcurr - y_pos_struc) ** 2)
            - radius_hcircle
            >= 0
        )
        zoneAB = (
            np.sqrt((xNcurr - x_pos_struc) ** 2 + (yNcurr - y_pos_struc) ** 2)
            - radius_hcircle
            < 0
        )
        if angleU != np.pi and angleU != 0.0 and angleU != 2 * np.pi:
            zoneB = (
                yNcurr
                - np.tan(angleU + np.pi / 2) * (xNcurr - x_pos_struc)
                - y_pos_struc
                > 0
            )
            zoneC = (
                yNcurr
                - np.tan(angleU + np.pi / 2) * (xNcurr - x_pos_struc)
                - y_pos_struc
                < 0
            )
        else:
            zoneB = xNcurr - x_pos_struc > 0
            zoneC = xNcurr - x_pos_struc < 0

        if zoneA:
            if zoneAA:
                radiusOC = radius_hcircle + thicknessU / 2.0
                LevelSet[itN] = (
                    np.sqrt((xNcurr - x_pos_struc) ** 2 + (yNcurr - y_pos_struc) ** 2)
                    - radiusOC
                )
                LevelSet_gradient_X[itN] = -(xNcurr - x_pos_struc) / (
                    LevelSet[itN] + radiusOC
                )
                LevelSet_gradient_Y[itN] = -(yNcurr - y_pos_struc) / (
                    LevelSet[itN] + radiusOC
                )
                LevelSet_gradient_R[itN] = -1
                LevelSet_gradient_T[itN] = 0
                IndicZone[itN] = 1
            elif zoneAB:
                radiusIC = radius_hcircle - thicknessU / 2.0
                LevelSet[itN] = -(
                    np.sqrt((xNcurr - x_pos_struc) ** 2 + (yNcurr - y_pos_struc) ** 2)
                    - radiusIC
                )
                LevelSet_gradient_X[itN] = (xNcurr - x_pos_struc) / (
                    -LevelSet[itN] + radiusIC
                )
                LevelSet_gradient_Y[itN] = (yNcurr - y_pos_struc) / (
                    -LevelSet[itN] + radiusIC
                )
                LevelSet_gradient_R[itN] = 1
                LevelSet_gradient_T[itN] = 0
                IndicZone[itN] = 2
        elif zoneC:
            radiusC2 = thicknessU / 2.0
            LevelSet[itN] = (
                np.sqrt((xNcurr - xc2) ** 2 + (yNcurr - yc2) ** 2) - radiusC2
            )
            LevelSet_gradient_X[itN] = -(xNcurr - xc2) / (LevelSet[itN] + radiusC2)
            LevelSet_gradient_Y[itN] = -(yNcurr - yc2) / (LevelSet[itN] + radiusC2)
            LevelSet_gradient_R[itN] = 0
            LevelSet_gradient_T[itN] = (
                -radius_hcircle * np.sin(angleU) * LevelSet_gradient_X[itN]
                + radius_hcircle * np.cos(angleU) * LevelSet_gradient_Y[itN]
            )
            IndicZone[itN] = 3
        elif zoneB:
            radiusC4 = thicknessU / 2.0
            LevelSet[itN] = (
                np.sqrt((xNcurr - xc4) ** 2 + (yNcurr - yc4) ** 2) - radiusC4
            )
            LevelSet_gradient_X[itN] = -(xNcurr - xc4) / (LevelSet[itN] + radiusC4)
            LevelSet_gradient_Y[itN] = -(yNcurr - yc4) / (LevelSet[itN] + radiusC4)
            LevelSet_gradient_R[itN] = 0
            LevelSet_gradient_T[itN] = (
                radius_hcircle * np.sin(angleU) * LevelSet_gradient_X[itN]
                - radius_hcircle * np.cos(angleU) * LevelSet_gradient_Y[itN]
            )
            IndicZone[itN] = 4

    LevelSetTangent = fluid_nodes[:, 1] - max(fluid_nodes[:, 1])

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_level_sets.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements},
            fields=[
                {
                    "data": LevelSetTangent,
                    "type": "nodal",
                    "dim": 1,
                    "name": "Tangent level set",
                },
                {"data": LevelSet, "type": "nodal", "dim": 1, "name": "Level set"},
                {
                    "data": LevelSet_gradient_X,
                    "type": "nodal",
                    "dim": 1,
                    "name": "Level set Grad X",
                },
                {
                    "data": LevelSet_gradient_Y,
                    "type": "nodal",
                    "dim": 1,
                    "name": "Level set Grad Y",
                },
                {
                    "data": LevelSet_gradient_R,
                    "type": "nodal",
                    "dim": 1,
                    "name": "Level set Grad R",
                },
                {
                    "data": LevelSet_gradient_T,
                    "type": "nodal",
                    "dim": 1,
                    "name": "Level set Grad T",
                },
                {"data": IndicZone, "type": "nodal", "dim": 1, "name": "Indic Zone"},
            ],
            append=True,
        )

    toc = time.process_time()
    logger.info("time to compute level set:", toc - tic)

    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.process_time()

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
    logger.info("time to find surface enriched elements:", toc - tic)

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
    logger.info("time to find edge enriched elements:", toc - tic)

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

    IIf, JJf, Vffk, Vffm = objFEM.getMatrices(
        fluid_nodes, fluid_elements, [celerity, rho]
    )

    KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = np.setdiff1d(list(range(fluid_ndof)), IdnodeS2 - 1)

    toc = time.process_time()
    logger.info("time to compute fluid matrices:", toc - tic)

    ##################################################################
    # Compute enrichment: Heaviside + Edge
    ##################################################################
    tic = time.process_time()

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
    logger.info("time to compute Heaviside enrichment:", toc - tic)

    Enrichednodes = np.unique(fluid_elements[EnrichedElements])

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
            [
                fluid_damping * KFF[SolvedDofF, :][:, SolvedDofF],
                fluid_damping * KAF[SolvedDofF, :][:, SolvedDofA],
            ],
            [
                fluid_damping * KAF[SolvedDofA, :][:, SolvedDofF],
                fluid_damping * KAA[SolvedDofA, :][:, SolvedDofA],
            ],
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

    # gradient wrt Xc
    IIf_X, JJf_X, Vfak_gradient_X, Vfam_gradient_X, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_X,
        material=[celerity, rho],
    )
    dMFA_dtheta_X = scipy.sparse.csc_matrix(
        (Vfam_gradient_X, (IIf_X, JJf_X)), shape=(fluid_ndof, fluid_ndof)
    )
    dKFA_dtheta_X = scipy.sparse.csc_matrix(
        (Vfak_gradient_X, (IIf_X, JJf_X)), shape=(fluid_ndof, fluid_ndof)
    )

    # gradient wrt Yc
    IIf_Y, JJf_Y, Vfak_gradient_Y, Vfam_gradient_Y, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_Y,
        material=[celerity, rho],
    )
    dMFA_dtheta_Y = scipy.sparse.csc_matrix(
        (Vfam_gradient_Y, (IIf_Y, JJf_Y)), shape=(fluid_ndof, fluid_ndof)
    )
    dKFA_dtheta_Y = scipy.sparse.csc_matrix(
        (Vfak_gradient_Y, (IIf_Y, JJf_Y)), shape=(fluid_ndof, fluid_ndof)
    )

    # gradient wrt R
    IIf_R, JJf_R, Vfak_gradient_R, Vfam_gradient_R, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_R,
        material=[celerity, rho],
    )
    dMFA_dtheta_R = scipy.sparse.csc_matrix(
        (Vfam_gradient_R, (IIf_R, JJf_R)), shape=(fluid_ndof, fluid_ndof)
    )
    dKFA_dtheta_R = scipy.sparse.csc_matrix(
        (Vfak_gradient_R, (IIf_R, JJf_R)), shape=(fluid_ndof, fluid_ndof)
    )

    # gradient wrt R
    IIf_T, JJf_T, Vfak_gradient_T, Vfam_gradient_T, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_T,
        material=[celerity, rho],
    )
    dMFA_dtheta_T = scipy.sparse.csc_matrix(
        (Vfam_gradient_T, (IIf_T, JJf_T)), shape=(fluid_ndof, fluid_ndof)
    )
    dKFA_dtheta_T = scipy.sparse.csc_matrix(
        (Vfak_gradient_T, (IIf_T, JJf_T)), shape=(fluid_ndof, fluid_ndof)
    )

    # export data
    with open(cwd / (results_file + "_KAF_MAF.pck"), "wb") as f:
        pickle.dump(
            [
                dMFA_dtheta_X[SolvedDofF, :][:, SolvedDofA],
                dKFA_dtheta_X[SolvedDofF, :][:, SolvedDofA],
                KAF[SolvedDofF, :][:, SolvedDofA],
                MAF[SolvedDofF, :][:, SolvedDofA],
            ],
            f,
        )

    dK_X = scipy.sparse.bmat(
        [
            [None, fluid_damping * dKFA_dtheta_X[SolvedDofF, :][:, SolvedDofA]],
            [fluid_damping * dKFA_dtheta_X[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dM_X = scipy.sparse.bmat(
        [
            [None, dMFA_dtheta_X[SolvedDofF, :][:, SolvedDofA]],
            [dMFA_dtheta_X[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dK_Y = scipy.sparse.bmat(
        [
            [None, fluid_damping * dKFA_dtheta_Y[SolvedDofF, :][:, SolvedDofA]],
            [fluid_damping * dKFA_dtheta_Y[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dM_Y = scipy.sparse.bmat(
        [
            [None, dMFA_dtheta_Y[SolvedDofF, :][:, SolvedDofA]],
            [dMFA_dtheta_Y[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dK_R = scipy.sparse.bmat(
        [
            [None, fluid_damping * dKFA_dtheta_R[SolvedDofF, :][:, SolvedDofA]],
            [fluid_damping * dKFA_dtheta_R[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dM_R = scipy.sparse.bmat(
        [
            [None, dMFA_dtheta_R[SolvedDofF, :][:, SolvedDofA]],
            [dMFA_dtheta_R[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dK_T = scipy.sparse.bmat(
        [
            [None, fluid_damping * dKFA_dtheta_T[SolvedDofF, :][:, SolvedDofA]],
            [fluid_damping * dKFA_dtheta_T[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    dM_T = scipy.sparse.bmat(
        [
            [None, dMFA_dtheta_T[SolvedDofF, :][:, SolvedDofA]],
            [dMFA_dtheta_T[SolvedDofA, :][:, SolvedDofF], None],
        ]
    )

    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis = 1
    FF = np.zeros(fluid_ndof)
    frequencies = []
    frf = []
    frfgradient_X = []
    frfgradient_Y = []
    frfgradient_R = []
    frfgradient_T = []

    if Flag_frf_analysis == 1:
        logger.info("time at the beginning of the FRF: {}".format(time.ctime()))

        press_save = []
        dpress_save_X = []
        dpress_save_Y = []
        dpress_save_R = []
        dpress_save_T = []
        disp_save = []

        # extract frequencies for the associated processors
        freqCompute = listFreqPerProc[:, rank]
        freqCompute = freqCompute[freqCompute > 0]

        for freq in freqCompute:
            # freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2 * np.pi * freq
            logger.info("proc number", rank, "frequency=", freq)

            FF[SolvedDofF] = -(
                KFF[SolvedDofF, :][:, IdnodeS2 - 1]
                - (omega**2) * MFF[SolvedDofF, :][:, IdnodeS2 - 1]
            ) * (np.ones((len(IdnodeS2))))
            FA = np.zeros(fluid_ndof)
            F = FF[SolvedDofF]
            F = np.concatenate((F, FA[SolvedDofA]))
            # F  = scipy.sparse.csc_matrix(F)

            ##
            pbFreq = K - (omega**2) * M

            sol = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), F
            )

            # stop
            tmp_X = -(dK_X - (omega**2) * dM_X) * sol
            tmp_Y = -(dK_Y - (omega**2) * dM_Y) * sol
            tmp_R = -(dK_R - (omega**2) * dM_R) * sol
            tmp_T = -(dK_T - (omega**2) * dM_T) * sol

            # Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
            Dsol_Dtheta_X = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), tmp_X
            )
            Dsol_Dtheta_Y = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), tmp_Y
            )
            Dsol_Dtheta_R = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), tmp_R
            )
            Dsol_Dtheta_T = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), tmp_T
            )

            press = np.zeros(fluid_ndof, dtype=complex)
            press[IdnodeS2 - 1] = np.ones(len(IdnodeS2))
            press[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            enrichment = np.zeros(fluid_nnodes, dtype=complex)
            enrichment[SolvedDofA] = sol[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            CorrectedPressure = press
            # logger.info(SolvedDofA)
            CorrectedPressure[SolvedDofA] = CorrectedPressure[SolvedDofA] + enrichment[
                SolvedDofA
            ] * np.sign(LevelSet[SolvedDofA])
            #####################
            #####################
            frf.append(
                objXFEM.getQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    CorrectedPressure,
                    0.0 * enrichment,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )
            #####################
            #####################
            press_save.append(CorrectedPressure.copy())
            #####################
            #####################
            Dpress_Dtheta = np.zeros([fluid_ndof, 4], dtype=complex)
            Dpress_Dtheta[SolvedDofF, 0] = Dsol_Dtheta_X[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF, 1] = Dsol_Dtheta_Y[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF, 2] = Dsol_Dtheta_R[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF, 3] = Dsol_Dtheta_T[list(range(len(SolvedDofF)))]
            #####################
            #####################
            Denrichment_Dtheta = np.zeros([fluid_ndof, 4], dtype=complex)
            Denrichment_Dtheta[SolvedDofA, 0] = Dsol_Dtheta_X[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            Denrichment_Dtheta[SolvedDofA, 1] = Dsol_Dtheta_Y[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            Denrichment_Dtheta[SolvedDofA, 2] = Dsol_Dtheta_R[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            Denrichment_Dtheta[SolvedDofA, 3] = Dsol_Dtheta_T[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            #####################
            #####################
            DCorrectedPressure_Dtheta = np.array(Dpress_Dtheta)
            DCorrectedPressure_Dtheta[SolvedDofA, 0] = DCorrectedPressure_Dtheta[
                SolvedDofA, 0
            ].T + np.array(
                Denrichment_Dtheta[SolvedDofA, 0] * np.sign(LevelSet[SolvedDofA]).T
            )
            DCorrectedPressure_Dtheta[SolvedDofA, 1] = DCorrectedPressure_Dtheta[
                SolvedDofA, 1
            ].T + np.array(
                Denrichment_Dtheta[SolvedDofA, 1] * np.sign(LevelSet[SolvedDofA]).T
            )
            DCorrectedPressure_Dtheta[SolvedDofA, 2] = DCorrectedPressure_Dtheta[
                SolvedDofA, 2
            ].T + np.array(
                Denrichment_Dtheta[SolvedDofA, 2] * np.sign(LevelSet[SolvedDofA]).T
            )
            DCorrectedPressure_Dtheta[SolvedDofA, 3] = DCorrectedPressure_Dtheta[
                SolvedDofA, 3
            ].T + np.array(
                Denrichment_Dtheta[SolvedDofA, 3] * np.sign(LevelSet[SolvedDofA]).T
            )
            #####################
            #####################
            frfgradient_X.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    CorrectedPressure,
                    enrichment * 0.0,
                    DCorrectedPressure_Dtheta[:, 0],
                    Denrichment_Dtheta[:, 0] * 0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )

            frfgradient_Y.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    CorrectedPressure,
                    enrichment * 0.0,
                    DCorrectedPressure_Dtheta[:, 1],
                    Denrichment_Dtheta[:, 1] * 0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )

            frfgradient_R.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    CorrectedPressure,
                    enrichment * 0.0,
                    DCorrectedPressure_Dtheta[:, 2],
                    Denrichment_Dtheta[:, 2] * 0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )

            frfgradient_T.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    CorrectedPressure,
                    enrichment * 0.0,
                    DCorrectedPressure_Dtheta[:, 3],
                    Denrichment_Dtheta[:, 3] * 0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )

            #####################
            #####################
            dpress_save_X.append(DCorrectedPressure_Dtheta[:, 0].copy())
            dpress_save_Y.append(DCorrectedPressure_Dtheta[:, 1].copy())
            dpress_save_R.append(DCorrectedPressure_Dtheta[:, 2].copy())
            dpress_save_T.append(DCorrectedPressure_Dtheta[:, 3].copy())

        #####################
        #####################
        logger.info("time at the end of the FRF: {}".format(time.ctime()))
        frfsave = [
            frequencies,
            frf,
            frfgradient_X,
            frfgradient_Y,
            frfgradient_R,
            frfgradient_T,
        ]
        if rank != 0:
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results == 1) and (rank == 0):
            msh2.mshWriter(
                cwd / (results_file + "_results_fluid_frf.msh"),
                fluid_nodes,
                [{"type": "TRI3", "connectivity": fluid_elements}],
                fields=[
                    {
                        "data": np.real(press_save),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure (real)",
                    },
                    {
                        "data": np.imag(press_save),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure (imaginary)",
                    },
                    {
                        "data": np.absolute(press_save),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure (norm)",
                    },
                    {
                        "data": np.real(dpress_save_X),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient X (real)",
                    },
                    {
                        "data": np.imag(dpress_save_X),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient X (imaginary)",
                    },
                    {
                        "data": np.absolute(dpress_save_X),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient X (norm)",
                    },
                    {
                        "data": np.real(dpress_save_Y),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient Y (real)",
                    },
                    {
                        "data": np.imag(dpress_save_Y),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient Y (imaginary)",
                    },
                    {
                        "data": np.absolute(dpress_save_Y),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient Y (norm)",
                    },
                    {
                        "data": np.real(dpress_save_R),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient R (real)",
                    },
                    {
                        "data": np.imag(dpress_save_R),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient R (imaginary)",
                    },
                    {
                        "data": np.absolute(dpress_save_R),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient R (norm)",
                    },
                    {
                        "data": np.real(dpress_save_T),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient T (real)",
                    },
                    {
                        "data": np.imag(dpress_save_T),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient T (imaginary)",
                    },
                    {
                        "data": np.absolute(dpress_save_T),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient T (norm)",
                    },
                ],
                append=True,
            )

        #####################
        #####################
        # save the FRF problem
        Allfrequencies = np.zeros(nbStep)
        Allfrf = np.zeros(nbStep)
        Allfrfgradient = np.zeros([nbStep, 4])
        k = 0
        if rank == 0:
            for i in range(nproc):
                if i == 0:
                    data = frfsave
                    logger.info(data)
                else:
                    logger.info(i)
                    data = comm.recv(source=i, tag=11)
                    # data=data_buffer

                for j in range(len(data[0])):
                    Allfrequencies[k] = data[0][j]
                    Allfrf[k] = data[1][j]
                    Allfrfgradient[k, 0] = data[2][j]
                    Allfrfgradient[k, 1] = data[3][j]
                    Allfrfgradient[k, 2] = data[4][j]
                    Allfrfgradient[k, 3] = data[5][j]
                    k = k + 1
            #####################
            #####################
            (
                Allfrequencies,
                Allfrf,
                Allfrfgradient_X,
                Allfrfgradient_Y,
                Allfrfgradient_R,
                Allfrfgradient_T,
            ) = zip(
                *sorted(
                    zip(
                        Allfrequencies,
                        Allfrf,
                        Allfrfgradient[:, 0],
                        Allfrfgradient[:, 1],
                        Allfrfgradient[:, 2],
                        Allfrfgradient[:, 3],
                    )
                )
            )
            Allfrfsave = [
                np.array(list(Allfrequencies)),
                np.array(list(Allfrf)),
                np.array(list(Allfrfgradient_X)),
                np.array(list(Allfrfgradient_Y)),
                np.array(list(Allfrfgradient_R)),
                np.array(list(Allfrfgradient_T)),
            ]  # ,press_save,dpress_save]
            with open(cwd / (results_file + "_results.frf"), "wb") as f:
                pickle.dump(Allfrfsave, f)
            #####################
            #####################
            # save on mat file
            scipy.io.savemat(
                cwd / (results_file + "_results.mat"), mdict={"AllFRF": Allfrfsave}
            )
            scipy.io.savemat(
                cwd / (results_file_ini + "results.mat"), mdict={"AllFRF": Allfrfsave}
            )
            #####################
            #####################


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
# function for dealing with options
def manageOpt(argv, dV):
    # load default values
    freqMin = dV.freqMin
    freqMax = dV.freqMax
    nbStep = dV.nbStep
    posStructX = dV.posStructX
    posStructY = dV.posStructY
    radiusStruct = dV.radiusStruct
    angleStruct = dV.angleStruct

    # load info from MPI
    nbProc, rank, comm = mpiInfo()
    # load options
    opts, args = getopt.getopt(argv, "p:s:F:f:hp")
    for opt, arg in opts:
        if opt == "-s":
            nbStep = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-p":
            tmp = np.array(arg.split(","), dtype=scipy.float32)
            posStructX = tmp[0]
            if len(tmp) > 1:
                posStructY = tmp[1]
            if len(tmp) > 2:
                radiusStruct = tmp[2]
            if len(tmp) > 3:
                angleStruct = tmp[3]
        elif opt == "-h":
            usage()
            sys.exit()
    # print chosen parameters
    print("Number of processors: ", nbProc)
    print("Number of frequency steps: ", nbStep)
    print("Maximum frequency: ", freqMax)
    print("Minimum frequency: ", freqMin)
    print("Structure position X: ", posStructX)
    print("Structure position Y: ", posStructY)
    print("Structure radius: ", radiusStruct)
    print("Structure angle: ", angleStruct)
    print("\n\n")

    # run computation
    RunPb(
        freqMin,
        freqMax,
        nbStep,
        nbProc,
        rank,
        comm,
        posStructX,
        posStructY,
        radiusStruct,
        angleStruct,
    )


# usage definition
def usage():
    dV = defaultV
    logger.info("Usage: ", sys.argv[0], "-psFfh [+arg]")
    logger.info("\t -p : number of processors (default value ", dV.nbProc, ")")
    logger.info(
        "\t -s : number of steps in the frequency range (default value ", dV.nbStep, ")"
    )
    logger.info("\t -F : maximum frequency (default value ", dV.freqMax, ")")
    logger.info("\t -f : minimum frequency (default value ", dV.freqMin, ")")


# default values
class defaultV:
    freqMin = 35.0
    freqMax = 80.0
    nbStep = 22
    posStructX = 1.5
    posStructY = 1
    radiusStruct = 0.75
    angleStruct = 0


### Run autonomous
if __name__ == "__main__":
    # run with options
    dV = defaultV
    manageOpt(sys.argv[1:], dV)
