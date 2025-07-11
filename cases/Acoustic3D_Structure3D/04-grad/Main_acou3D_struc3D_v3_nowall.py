###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
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
import os
from shutil import copyfile
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


# distribution of frequencies per processor
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
    freqMin, freqMax, nbStep, nbProc, rank, comm, saveResults=1
):  # , caseDefine):
    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##    Start SILEX vibro-acoustics computation   ##")
    logger.info("##################################################")
    logger.info("##################################################")

    # load 3D geometry
    orig_mesh_file = "geom/cavity_acou3D_struc_3D_v3_para"
    mesh_file = "geom/cavity_acou3D_struc_3D_v3"
    results_file_ini = "results/cavity_acou3D_struc_3D_v3"
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

    flag_write_gmsh_results = saveResults

    flag_edge_enrichment = 0

    # prepare save file
    results_file = results_file_ini + "_nowall"
    logger.info(results_file)

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
        logger.info("Number of nodes:", fluid_nnodes)
        logger.info("Number of elements in air:", fluid_nelem1)
        logger.info("Number of nodes in air:", fluid_nnodes1)

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
    # Construct the whole system
    #################################################################

    K = scipy.sparse.bmat([[fluid_damping * KFF[SolvedDofF, :][:, SolvedDofF]]])

    M = scipy.sparse.bmat([[MFF[SolvedDofF, :][:, SolvedDofF]]])

    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 1
    UF = np.zeros(2 * fluid_ndof, dtype=float)
    UF[9 - 1] = 3.1250e-05

    SolvedDof = np.hstack([SolvedDofF])

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
            logger.info("nb of total dofs: ", len(SolvedDofF))

        press_save = []
        disp_save = []

        # extract frequencies for the associated processors
        freqCompute = listFreqPerProc[:, rank]
        freqCompute = freqCompute[freqCompute > 0]
        it = 0
        itmax = len(freqCompute)
        for freq in freqCompute:
            it = it + 1
            # freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2 * np.pi * freq

            logger.info("Freq. step {} - proc number {} - frequency={}".format(it, rank, freq))

            tic = time.process_time()

            F = np.array(omega**2 * UF[SolvedDof], dtype="c16")

            if rank >= 0:
                # logger.info(K)
                # logger.info(M)
                # logger.info(omega)
                #  sol = mumps.spsolve(K-(omega**2)*M, F+0.j,comm=mycomm)
                sol = scipy.sparse.linalg.spsolve(
                    K - (omega**2) * M, F + 0.0j)

                # sol = mumps.spsolve(scipy.sparse.coo_matrix( \
                #   K-(omega**2)*M, dtype='complex'), F+0.j,comm=mycomm)
                # sol
                # sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(
                #     K-(omega**2)*M, dtype='c16'), F)

            ## pressure field without enrichment
            press1 = np.zeros((fluid_ndof), dtype=complex)
            press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            ## compute and store FRF on the test volume
            # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
            frf.append(
                acousticsXFEM.getQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    press1,
                    0.0 * press1,
                    np.real(0.0 * press1) + 1.0,
                    np.real(0.0 * press1) * 0 - 1.0,
                )
            )

            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(press1)

        frfsave = [frequencies, frf]
        if rank != 0:
            comm.send(frfsave, dest=0, tag=11)

        logger.info(
            "Proc. {} / time at the end of the FRF: {}".format(rank, time.ctime())
        )

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
            logger.info("Write pressure field and gradients in msh file")
            msh2.mshWriter(
                cwd / (results_file + str(rank) + "_results_fluid_frf.msh"),
                fluid_nodes,
                {"type": "TET4", "connectivity": fluid_elements1},
                fields=dataW,
                append=True,
            )

            logger.info(">>> Done!!")

        #####################
        #####################
        # save the FRF problem
        Allfrequencies = np.zeros(nbStep)
        Allfrf = np.zeros(nbStep)
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
            IXsort = np.argsort(Allfrequencies)
            AllfreqSorted = np.zeros(nbStep)
            AllfrfSorted = np.zeros(nbStep)
            for itS in range(0, nbStep):
                AllfreqSorted[itS] = Allfrequencies[IXsort[itS]]
                AllfrfSorted[itS] = Allfrf[IXsort[itS]]

            Allfrfsave = list()
            Allfrfsave.append(AllfreqSorted)
            Allfrfsave.append(AllfrfSorted)

            with open(cwd / (results_file + "_results.frf"), "wb") as f:
                pickle.dump(frfsave, f)
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
        elif opt == "-h":
            usage()
            sys.exit()
    # print chosen parameters
    print("Number of processors: ", nbProc)
    print("Number of frequency steps: ", nbStep)
    print("Maximum frequency: ", freqMax)
    print("Minimum frequency: ", freqMin)
    # print ("Case: ",caseDefine)
    print("\n\n")

    # run computation
    RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm)  # ,caseDefine)


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
    freqMax = 600.0
    nbStep = 5 # 2000
    nbProc = 1
    # caseDef= 'thick_u'


### Run autonomous
if __name__ == "__main__":
    # run with options
    dV = defaultV
    manageOpt(sys.argv[1:], dV)
