###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# Libraries
###########################################################
import getopt
import string
import time
import numpy
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io
import gmsh
import mumps

import pickle
from pathlib import Path
from loguru import logger

cwd = Path(__file__).resolve().parent

import sys
import os
from shutil import copyfile

from SILEXlib import silex_lib_xfem, silex_lib_fem
from SILEXlib import MeshField
from meshRW import msh2, msh

objFEM = silex_lib_fem.LinearAcousticsTET4()
objXFEM = silex_lib_xfem.LinearAcousticsTET4()
objDKT = silex_lib_fem.DKT()

from mpi4py import MPI


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


# distribution of frequencies per processor
def computeFreqPerProc(nbStep, nbProc, freqInit, freqEnd):
    # compute integer number of freq per proc and remaining steps
    nbFreqProc = nbStep // nbProc
    nbFreqProcRemain = nbStep % nbProc
    # compute frequencies steps
    varCase = 1
    if nbFreqProcRemain == 0:
        varCase = 0
    listFreq = numpy.zeros((nbFreqProc + varCase, nbProc))
    listAllFreq = numpy.linspace(freqInit, freqEnd, nbStep)
    # logger.info(numpy.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc + varCase):
            if itC * nbProc + itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1

    # logger.info(listFreq)
    return listFreq


# load structure mesh fo values of parameters
def buildStructMesh(fileOrig, destFile, paraVal):
    # copy original file to the used one
    copyfile(fileOrig + ".geo", destFile + ".geo")
    # change value of parameters in the new file
    for key, value in enumerate(paraVal):
        oldText = "<val##" + str(key) + ">"
        newText = "%g" % value
        # logger.info(oldText)
        # logger.info(newText)
        cmdSed = "sed -i 's/" + oldText + "/" + newText + "/g' " + destFile + ".geo"
        # logger.info(cmdSed)
        os.system(cmdSed)

    # run gmsh to build the mesh
    # os.system('gmsh -3 -format msh2 '+destFile+'.geo')


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

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################


def RunPb(
    freqMin, freqMax, nbStep, nbrefine, nbProc, rank, comm, saveResults=1
):  # , caseDefine):
    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##    Start SILEX vibro-acoustics computation   ##")
    logger.info("##################################################")
    logger.info("##################################################")

    # load 3D geometry
    fluid_mesh_file = "geom/rectangle_cavity"
    struct_mesh_file = "geom/simple_struct"
    results_file_ini = "results/rectangle_cavity"

    listFreqPerProc = computeFreqPerProc(nbStep, nbProc, freqMin, freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity = 340.0
    rho = 1.2
    fluid_damping = 1  # + 0.01j
    # shell structure
    material_Struc = []
    material_Struc.append(70000.0e6)  # E Young
    material_Struc.append(0.27)  # nu
    material_Struc.append(4.0e-3)  # thickness
    material_Struc.append(2700.0)  # rho

    dtypecustom = float
    if numpy.iscomplex(fluid_damping):
        dtypecustom = complex

    nproc = comm.Get_size()
    rank = comm.Get_rank()

    flag_write_gmsh_results = saveResults

    flag_edge_enrichment = 0

    results_file = results_file_ini
    logger.info(results_file)

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

    if rank == 0:
        gmsh.initialize()
        gmsh.open(str(cwd / (fluid_mesh_file + ".geo")))
        gmsh.model.mesh.generate(3)
        for _ in range(nbrefine):
            gmsh.model.mesh.refine()
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(cwd / (fluid_mesh_file + ".msh")))

    comm.Barrier()

    # start reading mesh
    mesh = msh.mshReader(cwd / (fluid_mesh_file + ".msh"))
    fluid_nodes = mesh.getNodes()[:, 0:3]
    fluid_elements1 = mesh.getElements(tag=1)["TET4"]  # air, cavity + controlled volume
    fluid_elements5 = mesh.getElements(tag=5)["TET4"]  # air, ONLY controlled volume
    IdNodes1 = numpy.unique(fluid_elements1)
    IdNodes5 = numpy.unique(fluid_elements5)

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
    mesh = msh.mshReader(cwd / (struct_mesh_file + ".msh"))
    struc_nodes = mesh.getNodes()[:, 0:3]
    struc_elements = mesh.getElements(tag=2)["TRI3"]

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

    # Find the fixed dofs and the free dofs of the structure
    tmpZa = scipy.sparse.find(struc_nodes[:, 2] == 0.0)  # z=0
    tmpZb = scipy.sparse.find(struc_nodes[:, 2] == 0.4)
    tmpYa = scipy.sparse.find(struc_nodes[:, 1] == 0.0)
    tmpYb = scipy.sparse.find(struc_nodes[:, 1] == 0.6)
    FixedStrucNodes = (
        numpy.unique(numpy.concatenate([tmpZa[0], tmpZb[0], tmpYa[0], tmpYb[0]])) + 1
    )

    FixedStrucDofUx = (FixedStrucNodes - 1) * 6
    FixedStrucDofUy = (FixedStrucNodes - 1) * 6 + 1
    FixedStrucDofUz = (FixedStrucNodes - 1) * 6 + 2
    # FixedStrucDofRx=(FixedStrucNodes-1)*6+3
    # FixedStrucDofRy=(FixedStrucNodes-1)*6+4
    # FixedStrucDofRz=(FixedStrucNodes-1)*6+5
    FixedStrucDofRx = []
    FixedStrucDofRy = []
    FixedStrucDofRz = []

    FixedStrucDof = numpy.hstack(
        [
            FixedStrucDofUx,
            FixedStrucDofUy,
            FixedStrucDofUz,
            FixedStrucDofRx,
            FixedStrucDofRy,
            FixedStrucDofRz,
        ]
    )

    SolvedDofS = numpy.setdiff1d(range(struc_ndof), FixedStrucDof)

    FS = numpy.zeros(struc_ndof)

    ##############################################################
    # Compute structure matrices
    ##############################################################
    tic = time.process_time()

    IIks, JJks, Vks, Vms = objDKT.getMatrices(
        struc_nodes, struc_elements, material_Struc
    )

    KSS = scipy.sparse.csc_matrix((Vks, (IIks, JJks)), shape=(struc_ndof, struc_ndof))
    MSS = scipy.sparse.csc_matrix((Vms, (IIks, JJks)), shape=(struc_ndof, struc_ndof))

    toc = time.process_time()
    if rank == 0:
        print("time for computing structure: {}".format(toc - tic))

    ##################################################################
    # compute level set
    ##################################################################

    tic = time.process_time()

    # LS from a structure mesh
    LevelSet, distance = objXFEM.getLevelSet(fluid_nodes, struc_nodes, struc_elements)

    toc = time.process_time()
    if rank == 0:
        logger.info("time to compute level set: {}".format(toc - tic))

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_struc_air_interface.msh"),
            struc_nodes,
            {"type": "TRI3", "connectivity": struc_elements},
        )
        # export levelset and levelset gradient
        dataW = list()
        dataW.append({"data": LevelSet, "type": "nodal", "dim": 1, "name": "Level set"})
        msh2.mshWriter(
            cwd / (results_file + "_LS_data.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1},
        )
    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.process_time()

    LSEnrichedElements = objXFEM.getEnrichedElements(
        fluidElements=fluid_elements1, levelset=LevelSet
    )

    LSEnrichednodes = numpy.unique(fluid_elements1[LSEnrichedElements])

    EnrichedElements = LSEnrichedElements
    toc = time.process_time()
    if rank == 0:
        logger.info("time to find enriched elements: {}".format(toc - tic))

    tic = time.process_time()

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_enriched_elements.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1[EnrichedElements]},
        )

        LS_moins_enriched = numpy.setdiff1d(LSEnrichedElements, EnrichedElements)
        enriched_moins_LS = numpy.setdiff1d(EnrichedElements, LSEnrichedElements)
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

    ##################################################################
    # Compute coupling STRUCTURE / AIR terms
    ##################################################################
    tic = time.process_time()

    IIc1, JJc1, Vc1 = objXFEM.getCoupling(
        fluidNodes=fluid_nodes,
        structNodes=struc_nodes,
        fluidElements=fluid_elements1,
        structElements=struc_elements,
        enrichedElements=EnrichedElements,
    )
    IIc2, JJc2, Vc2 = objXFEM.getCoupling(
        fluidNodes=fluid_nodes,
        structNodes=struc_nodes,
        fluidElements=fluid_elements1,
        structElements=struc_elements,
        enrichedElements=EnrichedElements,
        levelset=LevelSet,
    )

    CSA = 0.5 * scipy.sparse.csc_matrix(
        (Vc1, (IIc1, JJc1)), shape=(struc_ndof, fluid_ndof)
    ) + 0.5 * scipy.sparse.csc_matrix(
        (Vc2, (IIc2, JJc2)), shape=(struc_ndof, fluid_ndof)
    )

    toc = time.process_time()
    if rank == 0:
        logger.info("time to compute coupling matrices: {}".format(toc - tic))

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.process_time()

    IIf, JJf, Vffk, Vffm = objFEM.getMatrices(
        fluid_nodes, fluid_elements1, [celerity, rho]
    )

    KFF = scipy.sparse.csc_matrix((Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix((Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = list(range(fluid_ndof))

    ##################################################################
    # Compute Heaviside enrichment
    ##################################################################
    tic = time.process_time()

    Enrichednodes = numpy.unique(fluid_elements1[EnrichedElements])

    IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = objXFEM.getMatrices(
        fluid_nodes,
        fluid_elements1,
        LevelSet,
        [celerity, rho],
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
    IXFF = numpy.ix_(SolvedDofF, SolvedDofF)
    IXSS = numpy.ix_(SolvedDofS, SolvedDofS)
    IXFA = numpy.ix_(SolvedDofF, SolvedDofA)
    IXAF = numpy.ix_(SolvedDofA, SolvedDofF)
    IXAA = numpy.ix_(SolvedDofA, SolvedDofA)
    IXSA = numpy.ix_(SolvedDofS, SolvedDofA)
    IXAS = numpy.ix_(SolvedDofA, SolvedDofS)
    #
    K = scipy.sparse.bmat(
        [
            [KFF[IXFF], KAF[IXFA], None],
            [KAF[IXAF], KAA[IXAA], None],
            [None, -CSA[IXSA], KSS[IXSS]],
        ]
    )

    M = scipy.sparse.bmat(
        [
            [MFF[IXFF], MAF[IXFA], None],
            [MAF[IXAF], MAA[IXAA], CSA[IXSA].T],
            [None, None, MSS[IXSS]],
        ]
    )

    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 1
    UF = numpy.zeros(2 * fluid_ndof, dtype=float)
    UF[1 - 1] = 3.1250e-05

    SolvedDof = numpy.hstack([SolvedDofF, SolvedDofA + fluid_ndof])

    ##############################################################
    # FRF computation
    ##############################################################

    Flag_frf_analysis = 1
    frequencies = []
    frf = []

    if Flag_frf_analysis == 1:
        logger.info(
            "Proc. {} / time at the beginning of the FRF: {}".format(rank, time.ctime())
        )

        if rank == 0:
            logger.info("nb of total dofs: ", len(SolvedDofF) + len(SolvedDofA))

        press_save = []
        enrichment_save = []
        disp_save = []
        uncorrectedpress_save = []

        FA = numpy.zeros(len(SolvedDofA))
        FS = numpy.zeros(len(SolvedDofS))

        # extract frequencies for the associated processors
        freqCompute = listFreqPerProc[:, rank]
        freqCompute = freqCompute[freqCompute > 0]
        it = 0
        itmax = len(freqCompute)
        for freq in freqCompute:
            it = it + 1
            # freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2 * numpy.pi * freq

            logger.info(
                "Freq. step {} - proc number {} - frequency={}".format(it, rank, freq)
            )

            tic = time.process_time()

            FF = numpy.array(omega**2 * UF[SolvedDofF], dtype=dtypecustom)
            F = numpy.concatenate([FF, FA, FS],dtype=dtypecustom)
            if False:
                tic= time.process_time()
                solver = scipy.sparse.linalg.factorized(K - (omega**2) * M)
                sol = solver(numpy.array(F, dtype=dtypecustom))
                logger.info('time to solve the linear system: {}'.format(time.process_time()-tic))
            else:
                tic= time.process_time()
                ctx = mumps.Context()
                ctx.factor(K - (omega**2) * M)
                sol = ctx.solve(F)
                logger.info('time to solve the linear system: {}'.format(time.process_time()-tic))

            ## pressure field without enrichment
            press1 = numpy.zeros((fluid_ndof), dtype=dtypecustom)
            press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            ## enrichment field
            enrichment = numpy.zeros((fluid_nnodes), dtype=dtypecustom)
            enrichment[SolvedDofA] = sol[
                list(range(len(SolvedDofF), len(SolvedDofF) + len(SolvedDofA)))
            ]
            ## correction of the pressure field with enrichment
            CorrectedPressure = press1.copy()
            CorrectedPressure[SolvedDofA] = press1[SolvedDofA] + enrichment[
                SolvedDofA
            ] * numpy.sign(LevelSet[SolvedDofA])
            ## compute and store FRF on the test volume
            # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
            frf.append(
                objXFEM.getQuadraticPressure(
                    fluid_nodes,
                    fluid_elements1,
                    press1,
                    enrichment,
                    LevelSet,
                    LevelSet * 0 - 1.0,
                )
            )

            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(CorrectedPressure)
                enrichment_save.append(enrichment)
                uncorrectedpress_save.append(press1)

            #####################
            #####################
            ######################
            #####################

        frfsave = [frequencies, frf]
        if rank != 0:
            comm.send(frfsave, dest=0, tag=11)

        logger.info(
            "Proc. {} / time at the end of the FRF: {}".format(rank, time.ctime())
        )

        # prepare final data
        objMesh = MeshField.MeshField(
            fluid_nodes, fluid_elements1, LevelSet, LevelSet * 0.0 - 1.0
        )
        objMesh.addField(
            numpy.array(uncorrectedpress_save).transpose(),
            numpy.array(enrichment_save).transpose(),
        )

        # numpy.array(press_save),
        #             numpy.array(enrichment_save))
        # extract final data
        data = objMesh.getData(
            ["nodes", "mesh", "levelset", "levelset_tangent", "fields"]
        )

        if (flag_write_gmsh_results == 1) and (rank == 0):
            msh2.mshWriter(
                cwd / (results_file + "_results_fluid_frf.msh"),
                data["nodes"],
                [
                    {"type": "TET4", "connectivity": data["TET4"]},
                    {"type": "PRI6", "connectivity": data["PRI6"]},
                ],
                fields=[
                    {
                        "data": data["fields"],
                        "type": "nodal",
                        "nbsteps": itmax,
                        "name": "pressure",
                    },
                    {
                        "data": data["LS"],
                        "type": "nodal",
                        "nbsteps": 1,  # itmax,
                        "name": "levelset",
                    },
                    {
                        "data": data["statusnodes"],
                        "type": "nodal",
                        "nbsteps": 1,  # itmax,
                        "name": "statusnodes",
                    },
                ],
                append=True,
            )

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW = list()
            # prepare pressure field
            dataW.append(
                {
                    "data": numpy.real(press_save),
                    "type": "nodal",
                    "nbsteps": itmax,
                    "name": "pressure (real)",
                }
            )
            dataW.append(
                {
                    "data": numpy.imag(press_save),
                    "type": "nodal",
                    "nbsteps": itmax,
                    "name": "pressure (imaginary)",
                }
            )
            dataW.append(
                {
                    "data": numpy.absolute(press_save),
                    "type": "nodal",
                    "nbsteps": itmax,
                    "name": "pressure (norm)",
                }
            )
            logger.info("Write pressure field and gradients in msh file")
            msh2.mshWriter(
                cwd / (results_file + "_results_fluid_frf_raw.msh"),
                fluid_nodes,
                [{"type": "TET4", "connectivity": fluid_elements1}],
                fields=dataW,
                append=True,
            )

            logger.info(">>> Done!!")

            #

        #####################
        #####################
        # save the FRF problem
        Allfrequencies = numpy.zeros(nbStep)
        Allfrf = numpy.zeros(nbStep)
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
                    k = k + 1
            #####################
            #####################zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            IXsort = numpy.argsort(Allfrequencies)
            AllfreqSorted = numpy.zeros(nbStep)
            AllfrfSorted = numpy.zeros(nbStep)
            for itS in range(0, nbStep):
                AllfreqSorted[itS] = Allfrequencies[IXsort[itS]]
                AllfrfSorted[itS] = Allfrf[IXsort[itS]]

            # Allfrequencies, Allfrf,Allfrfgradient = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            Allfrfsave = list()
            Allfrfsave.append(AllfreqSorted)
            Allfrfsave.append(AllfrfSorted)

            data = {
                "freq": AllfreqSorted,
                "mag": AllfrfSorted,
                "nbnodes": fluid_nnodes,
                "nelem": fluid_nelem1,
            }
            with open("solvFRF" + str(fluid_nnodes).zfill(6) + ".pkl", "wb") as f:
                pickle.dump(data, f)
            #####################
            #####################
            # save on mat file
            scipy.io.savemat(results_file + "_results.mat", mdict={"data": data})
            scipy.io.savemat(results_file_ini + "results.mat", mdict={"data": data})
            #####################
            #####################
            return Allfrfsave


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
    nbRefine = dV.nbrefine
    # caseDefine = dV.caseDef

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
            tmp = numpy.array(arg.split(","), dtype=numpy.float32)
            paraVal = tmp
        elif opt == "-c":
            caseDefine = str(arg)
        elif opt == "-g":
            tmp = numpy.array(arg.split(","), dtype=numpy.int)
            gradCompute = tmp
        elif opt == "-h":
            usage()
            sys.exit()
    # print chosen parameters
    logger.info("Number of processors: ", nbProc)
    logger.info("Number of frequency steps: ", nbStep)
    logger.info("Maximum frequency: ", freqMax)
    logger.info("Minimum frequency: ", freqMin)
    # logger.info("Case: ",caseDefine)
    logger.info("\n\n")

    # run computation
    RunPb(freqMin, freqMax, nbStep, nbRefine, nbProc, rank, comm)  # ,caseDefine)


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
    freqMin = 10.0
    freqMax = 1000.0
    nbStep = 500
    nbrefine = 3
    # caseDef= 'thick_u'


### Run autonomous
if __name__ == "__main__":
    # run with options
    dV = defaultV
    manageOpt(sys.argv[1:], dV)
