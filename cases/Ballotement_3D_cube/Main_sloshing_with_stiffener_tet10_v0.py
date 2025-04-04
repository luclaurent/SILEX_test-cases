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


from SILEXlib import silex_lib_acou_tet10 as libF
#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_acou_tri6 as libFreeSurf
from SILEXlib import silex_lib_gmsh

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
mesh_file=Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10'
results_file=Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10'

dataPb['freq_ini'] = 0.1
dataPb['freq_ref'] = 0.1
dataPb['freq_end'] = 2.0
dataPb['nb_freq_step'] = 300

# Imposed acceleration on tank and stiffener surfaces 
dataPb['U_dot_dot_imposed'] = np.array([1.0,0.0,0.0])

# Flags
dataPb['flag_eigen_vectors'] = 1
dataPb['flag_FRF'] = 1
dataPb['flag_write_gmsh_results'] = 1

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

gmsh.finalize()

print('node 8 :', datafluidmesh['nodes'][8-1])

# gmsh output to check
if dataPb['flag_write_gmsh_results'] == 1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],11)
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],9)
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],9)



##############################################################
# Compute Standard Fluid Matrices : VOLUME
##############################################################

tic = time.process_time()
# get operators

IIf,JJf,Vffk,Vffm = libF.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                datafluidmesh['nodes'],
                                                1.0,
                                                1.0) # we put 1 for celerity and 1 for density

fluid_ndof = datafluidmesh['nodes'].shape[0]

HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

##############################################################
# Compute Standard Fluid Matrices : FREE SURFACE
##############################################################

print(libFreeSurf.globalacousticmatrices.__doc__)

IIf,JJf,Vffktmp,VSFF=libFreeSurf.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                             datafluidmesh['nodes'][:,[0,1]],
                                                             1.0,1.0)

SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



##############################################################
# Compute Standard Fluid load : rigid body motion of tank
##############################################################

print(libF.sloshimposedacc.__doc__)

CF,VecNormalElts = libF.sloshimposedacc(
                                np.array(datafluidmesh['nodes']),
                                np.array(datafluidmesh['structure_surface_elts']),
                                dataPb['U_dot_dot_imposed'])

CF=CF*dataFluid['rho']



# check if normal vectors are pointing out of the fluid volume
if dataPb['flag_write_gmsh_results'] == 1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces',
                                datafluidmesh['nodes'],
                                datafluidmesh['structure_surface_elts'],
                                9,
                                [[VecNormalElts,'elemental',3,'Normal to tank elements']])


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


##############################################################
# Compute eigen modes
##############################################################

if dataPb['flag_eigen_vectors']==1:

    eigen_values,eigen_vectors= scipy.sparse.linalg.eigsh(HFF,
                                                              20,
                                                              SFF,
                                                              sigma=0,which='LM')

    freq_eigv_S=list(np.sqrt(eigen_values)/(2*np.pi))
    print('Classic eigen frequencies : ',freq_eigv_S)

    eigen_vector_list=[]
    for i in range(eigen_values.shape[0]):
        Q=np.zeros(fluid_ndof)
        Q=eigen_vectors[:,i]
        eigen_vector_list.append(Q)


    if dataPb['flag_write_gmsh_results'] == 1:
        silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_Eigen_modes',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    11,
                                    [[eigen_vector_list,'nodal',1,'modes']])



##############################################################
# Compute FRF
##############################################################
if dataPb['flag_FRF'] == 1:

    press=[]
    frequencies=[]
    QuantityOfInterest=[]
    damping=None

    print ("Time at the beginning of the FRF:",time.ctime())
    for f in np.linspace(dataPb['freq_ini'],
                        dataPb['freq_end'],
                        dataPb['nb_freq_step']):

        print('Solve freq : ',f)
        frequencies.append(f)

        omega=2*np.pi*f
        forceType ='float'

    #    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
        sol = mumps.spsolve( HFF-omega**2*SFF , CF , comm=mycomm )
        press.append(sol.copy())

        QuantityOfInterest.append(sol[8-1]) # upper corner



    frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

    if dataPb['flag_write_gmsh_results'] == 1:
        print('QuantityOfInterest : ',QuantityOfInterest)
        silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    11,
                                    [[press,'nodal',1,'pressure']]
                                    )

    print ("Time at the end of the FRF:",time.ctime())

    f=open(results_file.as_posix() +'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()

