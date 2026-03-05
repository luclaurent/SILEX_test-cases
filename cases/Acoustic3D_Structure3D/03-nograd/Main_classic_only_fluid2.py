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


import mumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
acousticsFEM = silex_lib_fem.LinearAcousticsTET4()
acousticsXFEM = silex_lib_xfem.LinearAcousticsTET4()
structureFEM = silex_lib_fem.DKT()


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

#mesh_file='geom/cavity2'
#results_file='results/cavity2_damping'

import MakeCavityMesh3
angle=np.pi/8
MakeCavityMesh3.WriteGmshGeoFileAndMakeTheMesh('cavity3',[angle,1.5,2.0,3.0])

mesh_file='cavity3'
results_file='results/cavity3_damping'

freq_ini     = 20.0
freq_end     = 100.0
nb_freq_step_per_proc=100

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=340.0
rho=1.2
fluid_damping=(1.0+0.002j)
#fluid_damping=1.0

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements1,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1)
fluid_elements2,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2)
fluid_elements=np.vstack([fluid_elements1,fluid_elements2])
#fluid_elements_S2,IdNodesS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,2)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

print ("Number of fluid nodes:",fluid_nnodes)
print ("Number of fluid elements:",fluid_nelem)

#silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
#silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

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

P=np.zeros((fluid_ndof))
P[13-1]=1.0
#print(silex_lib_xfem_acou_tet4.forceonsurface.__doc__)
#P = silex_lib_xfem_acou_tet4.forceonsurface(fluid_nodes,fluid_elements_S2,1.0)


##############################################################
# FRF computation
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. {} / time at the beginning of the FRF: {}".format(rank, time.ctime()))

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
    #for freq in np.linspace(freq_ini,freq_end,nb_freq_step):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*np.pi*freq

        print ("proc number {} - frequency={}".format(rank,freq))

        F=np.array(omega**2*P , dtype='c16')
        #F=np.array(P , dtype='c16')

        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , np.array(F.todense() , dtype=complex) )
        sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )
        

        #press = np.zeros((fluid_ndof),dtype=float)
        press = np.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements2,fluid_nodes,press))
        #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
        i=i+1

        if rank==0:
            press_save.append(press.real)

    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if rank==0:
        silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

    print ("Proc. {} / time at the end of the FRF: {}".format(rank, time.ctime()))

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
        f=open(results_file+'_results.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()

