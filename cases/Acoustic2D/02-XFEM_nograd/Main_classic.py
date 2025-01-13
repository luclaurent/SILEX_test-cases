import string
import time
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
from SILEXlib import MeshField

#import xvibacoufo
#import shell_lib
import mumps
from mumps import DMumpsContext


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

freq_ini     = 100.5
freq_end     = 300.0
nb_freq_step_per_proc=50

nproc=comm.Get_size()
rank = comm.Get_rank()

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=1

freq_comparaison = 210.0

##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',2)
fluid_elements = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,1)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = len(np.unique(fluid_elements.flatten()))

fluid_elements_boun = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',1,2)
IdnodeS2=scipy.unique(fluid_elements_boun)

fluid_elements_tip = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',1,3)
IdnodeTip=[]
for e in range(len(fluid_elements_tip)):
    IdnodeTip.append(fluid_elements_tip[e][0])
IdnodeTip.append(fluid_elements_tip[len(fluid_elements_tip)-1][1])

IdnodeTip=np.array(IdnodeTip)


FF=np.zeros((fluid_ndof))

print "number of fluid nodes : ",fluid_nnodes
print "number of fluid elements : ",fluid_nelem

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_fluid_surface',fluid_nodes,fluid_elements,2)
    silex_lib_gmsh.WriteResults(results_file+'_fluid_boundary',fluid_nodes,fluid_elements_boun,1)
    silex_lib_gmsh.WriteResults(results_file+'_fluid_tip',fluid_nodes,fluid_elements_tip,1)

##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()
#print silex_acou_lib_tri3.globalacousticmatrices.__doc__

IIf,JJf,Vffk,Vffm=silex_acou_lib_tri3.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=np.setdiff1d(range(fluid_ndof),IdnodeS2-1)
#SolvedDofF=range(fluid_ndof)

##############################################################
# FRF computation of the FSI problem
##############################################################
enrichment=np.zeros((fluid_ndof))
LevelSet=np.zeros((fluid_ndof))+1.0
LevelSetTangent=np.zeros((fluid_ndof))-1.0

Flag_frf_analysis=1

FF=np.zeros((fluid_ndof))
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print "time at the beginning of the FRF:",time.ctime()

    press_save=[]
    disp_save=[]
    AllPressTipProcI=[]
    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*np.pi*freq

        FF[np.ix_(SolvedDofF)]=-(KFF[np.ix_(SolvedDofF,IdnodeS2-1)]-(omega*omega)*MFF[np.ix_(SolvedDofF,IdnodeS2-1)])*(np.zeros((len(IdnodeS2)))+1.0)
        print "proc number",rank,"frequency=",freq

        sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(KFF[np.ix_(SolvedDofF,SolvedDofF)]-(omega*omega)*MFF[np.ix_(SolvedDofF,SolvedDofF)])
                                          , np.array(FF[np.ix_(SolvedDofF)], dtype=float))
        
        #sol=mumps.spsolve( scipy.sparse.csc_matrix(KFF[np.ix_(SolvedDofF,SolvedDofF)]-(omega*omega)*MFF[np.ix_(SolvedDofF,SolvedDofF)])
        #                   , np.array(FF[np.ix_(SolvedDofF)], dtype=float)
        #                   )

        press = np.zeros(fluid_ndof)
        press[np.ix_(SolvedDofF)]=sol[range(len(SolvedDofF))]
        press_save.append(press)
        frf.append(silex_acou_lib_tri3.computequadratiquepressure(fluid_elements,fluid_nodes,press))
        #frf.append(xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press+0j,enrichment+0j,LevelSet,LevelSetTangent))


        if freq==freq_comparaison:
            PressTip=press[IdnodeTip-1]
            AllPressTip=[PressTip,fluid_nodes[IdnodeTip-1],freq]
            f=open(results_file+'_pressTip.frf','w')
            pickle.dump(AllPressTip, f)
            f.close()
        #thetaRef=[]
        #for i in range(len(IdnodeTip)):
        #    x=fluid_nodes[IdnodeTip[i]-1][0]-0.6
        #    y=fluid_nodes[IdnodeTip[i]-1][1]-0.65
        #    thetaRef.append(np.arctan(x/y)*180.0/np.pi)
        #pl.figure(1)
        #pl.plot(thetaRef,20*np.log10(np.real(abs(press[IdnodeTip-1]))/20e-6),'ko-',label='Reference', linewidth=2)
        #pl.show()

    print "time at the end of the FRF:",time.ctime()
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)
    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements,2,[[press_save,'nodal',1,'pressure']])

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
        f=open(results_file+'_results.frf','w')
        pickle.dump(Allfrfsave, f)
        f.close()

    
    



