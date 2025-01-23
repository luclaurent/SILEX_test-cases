import string
import time
from loguru import logger
from pathlib import Path
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen,cg
#from numpy.linalg import solve, norm

#import os
#from scipy.linalg import decomp
import pylab as pl
import pickle

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem

# load classes
objFEM = silex_lib_fem.LinearAcousticsTRI3()
objXFEM = silex_lib_xfem.LinearAcousticsTRI3()

#import xvibacoufo
#import shell_lib
import mumpspy


from mpi4py import MPI
comm = MPI.COMM_WORLD

# To run it:
#mpirun -np 4  python Main_classic.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='results/classic-test1'
results_file='results/classic-test1'

##############################################################
# Material, Boundary conditions
##############################################################

# air
celerity=340.0
rho=1.2

freq_ini     = 210.0
freq_end     = 300.0
nb_freq_step_per_proc=50

nproc=comm.Get_size()
rank = comm.Get_rank()

log_format =( 
    "<cyan> R{extra[rank]}</cyan> |"
    "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> | "
    "<level>{message}</level>"
)

logger.remove()
logger.configure(extra={"rank": 0})  # Default values
logger.add(sys.stdout, level='DEBUG', format=log_format, colorize=True, backtrace=True, diagnose=True)
logger = logger.bind(rank=rank)

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=1

freq_comparaison = 210.0

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

fluid_elements_tip = mesh.getElements(tag=3)["LIN2"]
IdnodeTip=np.unique(fluid_elements_tip.flatten())


FF=np.zeros((fluid_ndof))

logger.info("number of fluid nodes : {}".format(fluid_nnodes))
logger.info("number of fluid elements : {}".format(fluid_nelem))

if (flag_write_gmsh_results==1) and (rank==0):
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
    msh2.mshWriter(
        cwd / (results_file + "_fluid_tip.msh"),
        fluid_nodes,
        {"type": "LIN2", "connectivity": fluid_elements_tip},
    )
   
##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()
#print silex_acou_lib_tri3.globalacousticmatrices.__doc__

IIf, JJf, Vffk, Vffm = objFEM.getMatrices(fluid_nodes, fluid_elements, [celerity, rho])

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=np.setdiff1d(range(fluid_ndof),IdnodeS2-1)
#SolvedDofF=range(fluid_ndof)

toc = time.process_time()
logger.info("time to compute fluid matrices: {}".format(toc - tic))

##############################################################
# FRF computation of the FSI problem
##############################################################

Flag_frf_analysis=1

FF=np.zeros((fluid_ndof))
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    logger.info("time at the beginning of the FRF: {}".format(time.ctime())) #time.ctime()

    press_save=[]
    disp_save=[]
    AllPressTipProcI=[]
    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*np.pi*freq

        FF[np.ix_(SolvedDofF)]=-(KFF[np.ix_(SolvedDofF,IdnodeS2-1)]-(omega*omega)*MFF[np.ix_(SolvedDofF,IdnodeS2-1)])*(np.zeros((len(IdnodeS2)))+1.0)
        logger.info("proc number {} frequency={}".format(rank,freq))


        tic = time.process_time()
        sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(KFF[np.ix_(SolvedDofF,SolvedDofF)]-(omega*omega)*MFF[np.ix_(SolvedDofF,SolvedDofF)])
                                          , np.array(FF[np.ix_(SolvedDofF)], dtype=float))
        logger.info('time to solve the linear system: {}'.format(time.process_time()-tic))
        # tic = time.process_time()
        # solver = mumpspy.MumpsSolver()
        # sol = solver.solve(KFF[np.ix_(SolvedDofF,SolvedDofF)]-(omega*omega)*MFF[np.ix_(SolvedDofF,SolvedDofF)], FF[np.ix_(SolvedDofF)])
        # logger.info('time to solve the linear system: {}'.format(time.process_time()-tic))
        #sol=mumps.spsolve( scipy.sparse.csc_matrix(KFF[np.ix_(SolvedDofF,SolvedDofF)]-(omega*omega)*MFF[np.ix_(SolvedDofF,SolvedDofF)])
        #                   , np.array(FF[np.ix_(SolvedDofF)], dtype=float)
        #                   )

        press = np.zeros(fluid_ndof)
        press[np.ix_(SolvedDofF)]=sol[range(len(SolvedDofF))]
        press_save.append(press)
        frf.append(objFEM.getQuadraticPressure(
            fluid_nodes,
            fluid_elements,
            press
        ))
            

        #frf.append(xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press+0j,enrichment+0j,LevelSet,LevelSetTangent))


        if freq==freq_comparaison:
            PressTip=press[IdnodeTip-1]
            AllPressTip=[PressTip,fluid_nodes[IdnodeTip-1],freq]
            with open(cwd/(results_file+'_pressTip.frf'),'wb') as f:
                pickle.dump(AllPressTip, f)
            print('oui')
        #thetaRef=[]
        #for i in range(len(IdnodeTip)):
        #    x=fluid_nodes[IdnodeTip[i]-1][0]-0.6
        #    y=fluid_nodes[IdnodeTip[i]-1][1]-0.65
        #    thetaRef.append(np.arctan(x/y)*180.0/np.pi)
        #pl.figure(1)
        #pl.plot(thetaRef,20*np.log10(np.real(abs(press[IdnodeTip-1]))/20e-6),'ko-',label='Reference', linewidth=2)
        #pl.show()

    logger.info("time at the end of the FRF: {}".format(time.ctime()))
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)
    if (flag_write_gmsh_results==1) and (rank==0):
         msh2.mshWriter(
            cwd / (results_file +str(rank)+"_results_fluid_frf.msh"),
            fluid_nodes,
            [{"type": "TRI3", "connectivity": fluid_elements}],
            fields={
                "data": press_save,
                "type": "nodal",
                "dim": 1,
                "name": "pressure",
                "timesteps": frequencies,
            },
            append=True,
            )
        

    # save the FRF problem
    Allfrequencies=np.zeros(nb_freq_step)
    Allfrf=np.zeros(nb_freq_step)
    k=0
    if rank==0:
        for i in range(nproc):
            data = comm.recv(source=i, tag=11)
            for j in range(len(data[0])):
                Allfrequencies[k]=data[0][j]
                Allfrf[k]=data[1][j]
                k=k+1

        Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        Allfrfsave=[np.array(list(Allfrequencies)),np.array(list(Allfrf))]
        with open(cwd / (results_file+'_results.frf'),'wb') as f:
            pickle.dump(Allfrfsave, f)
        

    
    



