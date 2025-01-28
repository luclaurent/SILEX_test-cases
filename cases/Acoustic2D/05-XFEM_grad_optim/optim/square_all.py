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

import buildStruct as lvlB

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


def RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm, paraVal, caseDefine):
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

    # number of parameters
    nbPara = len(paraVal)
    # prepare save file
    file_extension = "{:.{}E}".format(paraVal[0], 2)
    if nbPara > 1:
        for i in range(1, nbPara):
            file_extension = file_extension + "_" + "{:.{}E}".format(paraVal[i], 2)

    results_file = results_file_ini + file_extension
    logger.info(results_file)

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

    # build level set
    (
        struc_nodes,
        struc_elements,
        LevelSet,
        LevelSetGradient,
        NamePara,
        LevelSetTangent,
        IndicZone,
    ) = lvlB.buildStruct(paraVal, fluid_nodes, caseDefine)

    # number of parameters
    nbPara = len(NamePara)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        # build fields to be exported
        fields = list()
        fields.append(
            {"data": LevelSet, "type": "nodal", "dim": 1, "name": "Level set"}
        )
        for itP, iN in enumerate(NamePara):
            fields.append(
                {
                    "data": LevelSetGradient[itP],
                    "type": "nodal",
                    "dim": 1,
                    "name": "Level Set Grad " + iN,
                }
            )

        fields.append(
            {"data": IndicZone, "type": "nodal", "dim": 1, "name": "Indic Zones"}
        )
        msh2.mshWriter(
            cwd / (results_file + "_level_sets.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements},
            fields=fields,
            append=True,
        )

    toc = time.process_time()
    logger.info("time to compute level set: {}".format(toc - tic))

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

    IIf, JJf, Vffk, Vffm = objFEM.getMatrices(
        fluid_nodes, fluid_elements, [celerity, rho]
    )

    KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = np.setdiff1d(list(range(fluid_ndof)), IdnodeS2 - 1)

    toc = time.process_time()
    logger.info("time to compute fluid matrices: {}".format(toc - tic))

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
    logger.info("time to compute Heaviside enrichment: {}".format(toc - tic))

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

    # build gradients matrices
    dK = list()
    dM = list()
    for itP in range(0, nbPara):
        logger.info(" Build Matrix for parameter " + NamePara[itP])
        IIf, JJf, Vfak_gradient, Vfam_gradient, _ = objXFEM.getGradientMatrices(
            fluid_nodes,
            fluid_elements[EnrichedElements],
            levelset=LevelSet,
            levelsetTg=None,
            levelsetGradient=LevelSetGradient[itP],
            material=[celerity, rho],
        )
        dMFA_dtheta = scipy.sparse.csc_matrix(
            (Vfam_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
        )
        dKFA_dtheta = scipy.sparse.csc_matrix(
            (Vfak_gradient, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof)
        )
        # export data
        with open(cwd / (results_file + "_KAF_MAF.pck"), "wb") as f:
            pickle.dump(
                [
                    dKFA_dtheta[SolvedDofF, :][:, SolvedDofA],
                    dMFA_dtheta[SolvedDofF, :][:, SolvedDofA],
                    KAF[SolvedDofF, :][:, SolvedDofA],
                    MAF[SolvedDofF, :][:, SolvedDofA],
                ],
                f,
            )
        # build full stiffness and mass gradient matrices
        dK.append(
            scipy.sparse.bmat(
                [
                    [None, fluid_damping * dKFA_dtheta[SolvedDofF, :][:, SolvedDofA]],
                    [fluid_damping * dKFA_dtheta[SolvedDofA, :][:, SolvedDofF], None],
                ]
            )
        )

        dM.append(
            scipy.sparse.bmat(
                [
                    [None, dMFA_dtheta[SolvedDofF, :][:, SolvedDofA]],
                    [dMFA_dtheta[SolvedDofA, :][:, SolvedDofF], None],
                ]
            )
        )

    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis = 1
    FF = np.zeros(fluid_ndof)
    frequencies = []
    frf = []
    frfgradient = list()
    for it in range(0, nbPara):
        frfgradient.append([])

    if Flag_frf_analysis == 1:
        logger.info("time at the beginning of the FRF: {}".format(time.ctime()))

        press_save = []
        dpress_save = list()
        for it in range(0, nbPara):
            dpress_save.append([])

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
            ## solve direct problem
            sol = scipy.sparse.linalg.spsolve(
                scipy.sparse.csc_matrix(pbFreq, dtype=complex), F
            )
            # store the pressure field
            press = np.zeros(fluid_ndof, dtype=complex)
            press[IdnodeS2 - 1] = np.ones(len(IdnodeS2))
            press[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            enrichment = np.zeros(fluid_nnodes, dtype=complex)
            enrichment[SolvedDofA] = sol[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            # compute the corrected pressure field (via enrichment)
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
                    press,
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
            ######################
            #####################
            Dpress_Dtheta = np.zeros([fluid_ndof, nbPara], dtype=complex)
            DCorrectedPressure_Dtheta = np.array(Dpress_Dtheta)
            Denrichment_Dtheta = np.zeros([fluid_ndof, nbPara], dtype=complex)
            ## compute gradients
            for itP in range(0, nbPara):
                tmpG = -(dK[itP] - (omega**2) * dM[itP]) * sol
                # Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
                Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve(
                    scipy.sparse.csc_matrix(pbFreq, dtype=complex), tmpG
                )
                #####################
                #####################
                Dpress_Dtheta[SolvedDofF, itP] = Dsol_Dtheta_RAW[
                    list(range(len(SolvedDofF)))
                ]
                #####################
                #####################
                Denrichment_Dtheta[SolvedDofA, itP] = Dsol_Dtheta_RAW[
                    list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
                ]
                #####################
                #####################
                # compute the corrected gradient pressure field (via enrichment)
                DCorrectedPressure_Dtheta[SolvedDofA, itP] = DCorrectedPressure_Dtheta[
                    SolvedDofA, itP
                ].T + np.array(
                    Denrichment_Dtheta[SolvedDofA, itP]
                    * np.sign(LevelSet[SolvedDofA]).T
                )
                #####################
                #####################
                # store gradients
                frfgradient[itP].append(
                    objXFEM.getGradientQuadraticPressure(
                        fluid_nodes,
                        fluid_elements5,
                        CorrectedPressure,
                        enrichment * 0.0,
                        DCorrectedPressure_Dtheta[:, itP],
                        Denrichment_Dtheta[:, itP] * 0.0,
                        LevelSet,
                        LevelSetTangent,
                        flag_edge_enrichment,
                    )
                )
                #####################
                #####################
                dpress_save[itP].append(DCorrectedPressure_Dtheta[:, itP].copy())

        #####################
        #####################
        logger.info("time at the end of the FRF: {}".format(time.ctime()))
        frfsave = [frequencies, frf, frfgradient]
        if rank != 0:
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW = list()
            # prepare pressure field
            dataW.append(
                {
                    "data": np.real(press_save),
                    "type": "nodal",
                    "nbsteps": len(freqCompute),
                    "name": "pressure (real)",
                }
            )
            dataW.append(
                {
                    "data": np.imag(press_save),
                    "type": "nodal",
                    "nbsteps": len(freqCompute),
                    "name": "pressure (imaginary)",
                }
            )
            dataW.append(
                {
                    "data": np.absolute(press_save),
                    "type": "nodal",
                    "nbsteps": len(freqCompute),
                    "name": "pressure (norm)",
                }
            )
            # prepare gradient pressure field
            itG = 0
            for itP in NamePara:
                dataW.append(
                    {
                        "data": np.real(dpress_save[itG]),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient " + itP + " (real)",
                    }
                )
                dataW.append(
                    {
                        "data": np.imag(dpress_save[itG]),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient " + itP + " (imaginary)",
                    }
                )
                dataW.append(
                    {
                        "data": np.absolute(dpress_save[itG]),
                        "type": "nodal",
                        "nbsteps": len(freqCompute),
                        "name": "pressure gradient " + itP + " (norm)",
                    }
                )
                itG = itG + 1
            # write fields
            msh2.mshWriter(
                cwd / (results_file + "_results_fluid_frf.msh"),
                fluid_nodes,
                [{"type": "TRI3", "connectivity": fluid_elements}],
                fields=dataW,
                append=True,
            )

        #####################
        #####################
        # save the FRF problem
        Allfrequencies = np.zeros(nbStep)
        Allfrf = np.zeros(nbStep)
        Allfrfgradient = np.zeros([nbStep, nbPara])
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
                    for itP in range(0, nbPara):
                        Allfrfgradient[k, itP] = data[2][itP][j]
                    k = k + 1
            #####################
            #####################
            Allfrequencies, Allfrf, Allfrfgradient = zip(
                *sorted(zip(Allfrequencies, Allfrf, Allfrfgradient))
            )
            Allfrfsave = list()
            Allfrfsave.append(np.array(list(Allfrequencies)))
            Allfrfsave.append(np.array(list(Allfrf)))
            for itP in range(0, nbPara):
                Allfrfsave.append(np.array(list(Allfrfgradient[itP])))

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
# function for dealing with options
def manageOpt(argv, dV):
    # load default values
    freqMin = dV.freqMin
    freqMax = dV.freqMax
    nbStep = dV.nbStep
    paraVal = np.array(dV.paraVal)
    caseDefine = dV.caseDef

    # load info from MPI
    nbProc, rank, comm = mpiInfo()
    # load options
    opts, args = getopt.getopt(argv, "p:s:F:f:hp:c:")
    for opt, arg in opts:
        if opt == "-s":
            nbStep = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-p":
            tmp = np.array(arg.split(","), dtype=scipy.float32)
            paraVal = tmp
        elif opt == "-c":
            caseDefine = str(arg)
        elif opt == "-h":
            usage()
            sys.exit()
    # print chosen parameters
    print("Number of processors: ", nbProc)
    print("Number of frequency steps: ", nbStep)
    print("Maximum frequency: ", freqMax)
    print("Minimum frequency: ", freqMin)
    print("Case: ", caseDefine)
    it = 0
    for itP in paraVal:
        print("Parameter num " + str(it) + ": " + str(itP))
        it = it + 1
    print("\n\n")

    # run computation
    RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm, paraVal, caseDefine)


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
    paraVal = [1.5, 1, 0.75, 0]
    caseDef = "thick_u"


### Run autonomous
if __name__ == "__main__":
    # run with options
    dV = defaultV
    manageOpt(sys.argv[1:], dV)
