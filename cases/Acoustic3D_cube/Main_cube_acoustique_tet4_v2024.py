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


def loadLib():
    from SILEXlib import silex_lib_acou_tet4 as lib
    return lib
try: 
    lib = loadLib()
except:
    lib = None    


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
data = dict()

# parallepipedic cavity with plane structure
mesh_file=Path(__file__).parent / 'cube_acou_tet4'
results_file=Path(__file__).parent / 'cube_acou_tet4'

data['freq_ini'] = 150.0
data['freq_ref'] = 237.0
data['freq_end'] = 500.0
data['nb_freq_step'] = 100
#deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

# air
data['celerity'] = 343.0
data['rho'] = 1.2
data['fluid_damping'] = 1.0

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
external_surface = (2,6) # external surface physical group
full_volume = (3,1) # physical group of the full volume (air: cavity + control volume)   
control_volume = (3,5) # physical group of the control volume (air: only control volume)
if not external_surface in physical_groups:
    raise ValueError('External surface physical group not found')
if not full_volume in physical_groups:
    raise ValueError('Full volume physical group not found')
if not control_volume in physical_groups:
    raise ValueError('Control volume physical group not found')

# get entities per physical group
external_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*external_surface)
full_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*full_volume)
control_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*control_volume)
# get elements
external_surface_elements = u.getElementsFromEntities(gmsh,external_surface[0],external_surface_entities)
full_volume_elements = u.getElementsFromEntities(gmsh,full_volume[0],full_volume_entities)
control_volume_elements = u.getElementsFromEntities(gmsh,control_volume[0],control_volume_entities)

print(' ----  full_volume_elements : ',full_volume_elements)

datafluidmesh = dict()
datafluidmesh['nodes'] = fluid_nodes
datafluidmesh['full_volume_elts'] = full_volume_elements
datafluidmesh['control_volume_elts'] = control_volume_elements
datafluidmesh['external_surface_elts'] = external_surface_elements

#fluid_elements = full_volume_elements # air, cavity + controlled volume
#fluid_elements5 = control_volume_elements # air, ONLY controlled volume
#,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
#fluid_elements_S6 = external_surface_elements# external surface
#,IdNodesS6 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,6) # external surface

#fluid_nnodes   = fluid_nodes.shape[0]
#fluid_nelem    = fluid_elements.shape[0]
#fluid_ndof     = fluid_nnodes

#fluid_nnodes5 = IdNodes5.shape[0]
#fluid_nnodes6 = IdNodesS6.shape[0]

#print ("Number of fluid nodes:",fluid_nnodes)
#print ("Number of fluid elements:",fluid_nelem)

#silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
#silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

gmsh.finalize()

# LOAD
data['loaddof'] = 0
data['valload'] = 1


##############################################################
# Compute Standard Fluid Matrices
##############################################################

tic = time.process_time()
# get operators

IIf,JJf,Vffk,Vffm = lib.globalacousticmatrices(datafluidmesh['full_volume_elts'],
                                                datafluidmesh['nodes'],
                                                data['celerity'],
                                                data['rho'])

fluid_ndof = datafluidmesh['nodes'].shape[0]

KFF,MFF,F = ua.assembleOP(IIf,JJf,Vffk,Vffm,
                  fluid_ndof,
                  data['loaddof'],
                  data['valload'])

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

