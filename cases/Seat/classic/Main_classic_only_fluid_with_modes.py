import string
import time
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
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt
import silex_lib_gmsh

import mumps

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4 python3 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='geom/sieges_classic_fluid_cavity'
results_file='results/classic_classic_fluid_cavity_with_modes'


freq_ini     = 20.0
freq_end     = 400.0
nb_freq_step_per_proc=50

nb_mode_F=50

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
celerity=340.0
rho=1.2
fluid_damping=(1.0+0.002j)

##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1)
fluid_elements_S2,IdNodesS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,2)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

#silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
#silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=list(range(fluid_ndof))
##################################################################
# Compute eigen modes of the fluid
##################################################################
tic = time.clock()
eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(KFF[SolvedDofF,:][:,SolvedDofF],nb_mode_F,MFF[SolvedDofF,:][:,SolvedDofF],sigma=0,which='LM')

freq_eigv_F=list(scipy.sqrt(eigen_values_F)/(2*scipy.pi))
eigen_vector_F_list=[]
for i in range(nb_mode_F):
    tmp=eigen_vectors_F[:,i].real
    eigen_vector_F_list.append(tmp[SolvedDofF])

toc = time.clock()

if rank==0:
    print ("fluid eigen frequencies : ",freq_eigv_F[0]," to ",freq_eigv_F[nb_mode_F-1])
    #silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes,fluid_elements,4,[[eigen_vector_F_list,'nodal',1,'pressure']])
    print ("time for computing the fluid modes:",toc-tic)


VK_diag_nn = eigen_values_F
VM_diag_nn = eigen_values_F/eigen_values_F
IIDnn = list(range(nb_mode_F))
JJDnn = list(range(nb_mode_F))
K_diag_nn= scipy.sparse.csc_matrix( (VK_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_F,nb_mode_F) )
M_diag_nn= scipy.sparse.csc_matrix( (VM_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_F,nb_mode_F) )

##################################################################
# Construct the whole system
##################################################################

##K=KFF[SolvedDofF,:][:,SolvedDofF]
##M=MFF[SolvedDofF,:][:,SolvedDofF]
K=K_diag_nn
M=M_diag_nn

# To impose the load on the fluid:
# fluid node number 1 is at (0,-ly/2,0)
# node number 1 is at (0,-ly/2,0)
#F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )

#P=scipy.zeros((fluid_ndof))
#P[0]=1.0
#print(silex_lib_xfem_acou_tet4.forceonsurface.__doc__)

P = silex_lib_xfem_acou_tet4.forceonsurface(fluid_nodes,fluid_elements_S2,1.0)
P = scipy.dot(eigen_vectors_F.T,P)



##############################################################
# FRF computation
##############################################################



Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
    #for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq

        print ("proc number",rank,"frequency=",freq)

        F=scipy.array(omega**2*P , dtype='c16')
        #F=scipy.array(P , dtype='c16')

        #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , scipy.array(F.todense() , dtype=complex) )
        sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega*omega)*M,dtype='c16') , F , comm=mycomm )
        #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype='float') , F , comm=mycomm )
        

        #press = scipy.zeros((fluid_ndof),dtype=float)
        press = scipy.zeros((fluid_ndof),dtype=complex)
        press[SolvedDofF]=scipy.dot(eigen_vectors_F,sol[list(range(nb_mode_F))])
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,press))
        #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
        i=i+1

        if rank==0:
            press_save.append(press.real)

    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if rank==0:
        silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

    print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

    # save the FRF problem
    Allfrequencies=scipy.zeros(nb_freq_step)
    Allfrf=scipy.zeros(nb_freq_step)
    k=0
    if rank==0:
        for i in range(nproc):
            data = comm.recv(source=i, tag=11)
            for j in range(len(data[0])):
                Allfrequencies[k]=data[0][j]
                Allfrf[k]=data[1][j]
                k=k+1

        Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
        f=open(results_file+'_results_with_damping.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()

