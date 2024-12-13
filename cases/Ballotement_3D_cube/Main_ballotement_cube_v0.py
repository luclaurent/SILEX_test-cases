import string
import time
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')
import mumps
import gmsh
import utils as u
import utils_acoustics as ua


from SILEXlib import silex_lib_acou_tet4 as libF
#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_acou_tri3 as libFreeSurf

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
dataPb = dict()

# parallepipedic cavity with plane structure
mesh_file=Path(__file__).parent / 'cube_ballotement_tet4'
results_file=Path(__file__).parent / 'cube_ballotement_tet4'

dataPb['freq_ini'] = 150.0
dataPb['freq_ref'] = 237.0
dataPb['freq_end'] = 500.0
dataPb['nb_freq_step'] = 5
#deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# fluid
dataFluid = dict()
#data['celerity'] = 343.0
dataFluid['rho'] = 1000.0
#data['fluid_damping'] = 1.0

# structure
#dataStructure = dict()
#dataStructure['rho'] = 2700.0
#dataStructure['thickness'] = 5.0e-3
#dataStructure['young'] = 72000.0e6
#dataStructure['nu'] = 0.3

# LOAD on structure
#dataPb['loaddof'] = 0
#dataPb['valload'] = 1


##############################################################
# Load fluid mesh
##############################################################

tic = time.process_time()

gmsh.initialize()
gmsh.open(mesh_file.as_posix()+'.geo')
# mesh generation
gmsh.model.mesh.generate(3)

# get fluid nodes
_,_nodes,_ = gmsh.model.mesh.getNodes()
fluid_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))

# check physical groups
physical_groups = gmsh.model.getPhysicalGroups()
structure_surface = (2,20) # structure surface physical group
free_fluid_surface = (2,30) # free fluid surface physical group
fluid_volume = (3,10) # fluid physical group
if not structure_surface in physical_groups:
    raise ValueError('Structure surface physical group not found')
if not fluid_volume in physical_groups:
    raise ValueError('Fluid volume physical group not found')
if not free_fluid_surface in physical_groups:
    raise ValueError('Free fluid surface physical group not found')


# get entities per physical group
structure_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*structure_surface)
free_fluid_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*free_fluid_surface)
fluid_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*fluid_volume)

# get elements
structure_surface_elements = u.getElementsFromEntities(gmsh,structure_surface[0],structure_surface_entities)
free_fluid_surface_elements = u.getElementsFromEntities(gmsh,free_fluid_surface[0],free_fluid_surface_entities)
fluid_volume_elements = u.getElementsFromEntities(gmsh,fluid_volume[0],fluid_volume_entities)

print(' ----  structure_surface_elements : ',structure_surface_elements)
print(' ----  free_fluid_surface_elements : ',free_fluid_surface_elements)
print(' ----  fluid_volume_elements : ',fluid_volume_elements)

datafluidmesh = dict()
datafluidmesh['nodes'] = fluid_nodes
datafluidmesh['fluid_volume_elts'] = fluid_volume_elements
datafluidmesh['structure_surface_elts'] = structure_surface_elements
datafluidmesh['free_fluid_surface_elts'] = free_fluid_surface_elements
#datastructuremesh = dict()
#datastructuremesh['structure_surface_elts'] = structure_surface_elements
#datastructuremesh['nodes'] = fluid_nodes # !!!! attention aux tailles des matrices par la suite !!!!!

#silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
#silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

gmsh.finalize()


##############################################################
# Compute Standard Fluid Matrices : VOLUME
##############################################################

tic = time.process_time()
# get operators

IIf,JJf,Vffk,Vffm = libF.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                datafluidmesh['nodes'],
                                                1,
                                                1) # we put 1 for celerity and 1 for density

fluid_ndof = datafluidmesh['nodes'].shape[0]

HFF,MFFtmp,Ftmp = ua.assembleOP(IIf,JJf,Vffk,Vffm,
                  fluid_ndof,
                  0,
                  0)

##############################################################
# Compute Standard Fluid Matrices : FREE SURFACE
##############################################################



IIf,JJf,Vffktmp,VSFF=libFreeSurf.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                             datafluidmesh['nodes'][:,0:2],
                                                             1.0,1.0)

SFF,MFFtmp,Ftmp = ua.assembleOP(IIf,JJf,Vffktmp,VSFF,
                  fluid_ndof,
                  0,
                  0)

##############################################################
# Compute Standard Fluid load : rigid body motion of tank
##############################################################
print('blabla')
#print(libF.SloshImposedAcc.__doc__)

CF = libF.sloshimposedacc(
    np.array(datafluidmesh['nodes']),
                                 np.array(datafluidmesh['structure_surface_elts']),
                                 np.array([1.0,0.0,0.0]))

Stop


##############################################################
# Compute Standard Strcuture Matrices
##############################################################

#IIs,JJs,Vks,Vms=libS.stiffnessmatrix(datastructuremesh['nodes'],
#                                       datastructuremesh['structure_surface_elts'],
#                                       [dataStructure['young'],dataStructure['nu'],dataStructure['thickness'],dataStructure['rho']])
#
#KSS,MSS,FS = ua.assembleOP(IIs,JJs,Vks,Vms,
#                  fluid_ndof*6,
#                  dataPb['loaddofS'],
#                  dataPb['valload'])
#
## Il faut ici enlever les lignes et colonnes de KSS MSS FS qui ne sont pas utilisees ...
#NodesSlist = list(range(1,len(np.unique(datastructuremesh['structure_surface_elts']))+1))
#
#dofS=np.hstack([np.array(NodesSlist),np.array(NodesSlist)+1,np.array(NodesSlist)+2,np.array(NodesSlist)+3,np.array(NodesSlist)+4,np.array(NodesSlist)+5]) 

STOP
##############################################################
# Compute FRF
##############################################################
press=[]
frequencies=[]
SPL=[]
damping=None
#SolvedDofF=list(range(fluid_ndof))
print ("Time at the beginning of the FRF:",time.ctime())
for f in np.linspace(data['freq_ini'],
                    data['freq_end'],
                    data['nb_freq_step']):
    
    print('Solve freq : ',f)
    frequencies.append(f)
#    
#    press.append(ua.solveAcousticsOneFreq(f,
#                                            KFF,
#                                            MFF,
#                                            F,damping,FlagMumps=1))
#    
    omega=2*np.pi*f
    FF = omega**2*F
    forceType ='float'
#    if damping is None:
#        forceType ='float'
#    else:
#        forceType ='complex'
#        KFF = (1+damping*1j) * KFF
    # solve
#    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
    sol = mumps.spsolve( KFF-omega**2*MFF , FF.todense() , comm=mycomm )
    press.append(sol.copy())
    
    # compute Mean Quadaratic Pressure
    if damping is None:
        SPL.append(lib.computequadratiquepressure(datafluidmesh['control_volume_elts'],
                                                datafluidmesh['nodes'],
                                                press[-1]))
    else:
        SPL.append(lib.computecomplexquadratiquepressure(datafluidmesh['control_volume_elts'],
                                                        datafluidmesh['nodes'],
                                                        press[-1]))

frfsave=[np.array(frequencies),np.array(SPL)]

print('SPL : ',SPL)
#silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

print ("Time at the end of the FRF:",time.ctime())

f=open(results_file.as_posix() +'_results.frf','wb')
pickle.dump(frfsave, f)
f.close()

