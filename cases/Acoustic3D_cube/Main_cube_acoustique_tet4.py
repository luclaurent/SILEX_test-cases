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
objFEM = silex_lib_fem.LinearAcousticsTET4()

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

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4 python3.4 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3.4 Main_toto.py
#
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='geom/cube_acou_tet4'
results_file='results/cube_acou_tet4'
cwd = Path(__file__).resolve().parent


freq_ini     = 150.0
freq_end     = 1000.0
nb_freq_step = 20

#deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=340.0
rho=1.2
#fluid_damping=(1.0+0.002j)
fluid_damping=1.0

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

mesh = msh.mshReader(cwd / (mesh_file + "_fluid.msh"))

fluid_nodes = mesh.getNodes()[:, 0:3]
fluid_elements = mesh.getElements(tag=1)["TET4"]
fluid_elements5 = mesh.getElements(tag=5)["TET4"]
fluid_elements6 = mesh.getElements(tag=6)["TRI3"]
Idnodes = np.unique(fluid_elements.flatten())
Idnodes5 = np.unique(fluid_elements5.flatten())
Idnodes6 = np.unique(fluid_elements6.flatten())

fluid_nnodes = fluid_nodes.shape[0]
fluid_nelem = fluid_elements.shape[0]
fluid_ndof = len(np.unique(fluid_elements.flatten()))


if rank == 0:
    msh2.mshWriter(
        cwd / (results_file + "Mesh.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements},
    )
    msh2.mshWriter(
        cwd / (results_file + "Mesh_control_volume.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements5},
    )
    msh2.mshWriter(
        cwd / (results_file + "Mesh_surface.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements6},
    )
logger.info("nnodes for fluid= {}".format(fluid_nnodes))
logger.info("nelem for fluid= {}".format(fluid_nelem))
logger.info("nelem for control volume= {}".format(fluid_elements5.shape[0]))
logger.info("nelem for surface= {}".format(fluid_elements6.shape[0]))

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf, JJf, Vffk, Vffm = objFEM.getMatrices(fluid_nodes, fluid_elements, [celerity, rho])

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=list(range(fluid_ndof))

##################################################################
# Construct the whole system
##################################################################

K=KFF[SolvedDofF,:][:,SolvedDofF]
M=MFF[SolvedDofF,:][:,SolvedDofF]

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )

F0=np.zeros((fluid_ndof))
F0[1-1]=1.0

##############################################################
# FRF computation
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    logger.info("Proc. {} / time at the beginning of the FRF: {}".format(rank, time.ctime()))

    press_save=[]
    disp_save=[]

    for freq in np.linspace(freq_ini,freq_end,nb_freq_step):

        frequencies.append(freq)
        omega=2*np.pi*freq

        logger.info("frequency= {}".format(freq))

        F=np.array(omega**2*F0 , dtype='d')
        #F=np.array(omega**2*F0 , dtype='c16')
        #F[SolvedDofF]=-(KFF[SolvedDofF,:][:,1-1]-(omega**2)*MFF[SolvedDofF,:][:,1-1])*(np.zeros((len([1])))+1.0)


        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , np.array(F.todense() , dtype=complex) )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )

        sol = pymumps.spsolve(K-(omega**2)*M , F , comm=mycomm )
        #sol = scipy.sparse.linalg.spsolve(K-(omega**2)*M , F )

        #sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        

        #press = np.zeros((fluid_ndof),dtype=float)
        press = np.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        frf.append(objFEM.getQuadraticPressure(fluid_nodes,fluid_elements5, press))

        if rank==0:
            press_save.append(press.real)

    frfsave=[np.array(frequencies),np.array(frf)]

    if rank==0:
        msh2.mshWriter(
            cwd / (results_file + "_results_fluid_frf.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements},
            fields={"name": "pressure", 
                    "nbsteps": nb_freq_step,
                    "data": press_save, 
                    "type": "nodal"},
            append=True,
        )
    logger.info("Time at the end of the FRF: {}".format(time.ctime()))

    Allfrfsave=[frequencies,frf]
    with open(cwd / (results_file + "_results.frf"), "wb") as f:
        pickle.dump(frfsave, f)

