import string
import time
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.sparse.construct

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
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_xfem
#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_acou_tri3 as libFreeSurf
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
mesh_file_fluid=Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet4'
mesh_file_stiffener=Path(__file__).parent / 'cube_xfem_sloshing_Stiffener_tet4'

results_file=Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet4'

dataPb['freq_ini'] = 0.1
dataPb['freq_ref'] = 0.1
dataPb['freq_end'] = 2.0
dataPb['nb_freq_step'] = 100

# Imposed acceleration on tank and stiffener surfaces 
dataPb['U_dot_dot_imposed'] = np.array([1.0,0.0,0.0])

# Flags
dataPb['flag_eigen_vectors'] = 0
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
# Load fluid mesh and tank
##############################################################

tic = time.process_time()

gmsh.initialize()
gmsh.open(mesh_file_fluid.as_posix()+'.geo')
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

fluid_ndof = datafluidmesh['nodes'].shape[0]

print('node 8 :', datafluidmesh['nodes'][8-1])

# gmsh output to check
if dataPb['flag_write_gmsh_results']==1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],4)
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],2)
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],2)

##############################################################
# Load stiffener mesh / compute Level Set / Get enriched nodes and elements
##############################################################

tic = time.process_time()

gmsh.initialize()
gmsh.open(mesh_file_stiffener.as_posix()+'.geo')
# mesh generation
gmsh.model.mesh.generate(3)

# get fluid nodes
_,_nodes,_ = gmsh.model.mesh.getNodes()
stiffener_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))

# check physical groups
physical_groups = gmsh.model.getPhysicalGroups()

stiffener_surface = (2,50) # stiffener surface physical group
stiffener_edge = (1,60) # stiffener edge physical group

if not stiffener_surface in physical_groups:
    raise ValueError('stiffener surface physical group not found')
if not stiffener_edge in physical_groups:
    raise ValueError('stiffener edge physical group not found')

# get entities per physical group
stiffener_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_surface)
stiffener_edge_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_edge)

# get elements
stiffener_surface_elements = u.getElementsFromEntities(gmsh,stiffener_surface[0],stiffener_surface_entities)
stiffener_edge_elements = u.getElementsFromEntities(gmsh,stiffener_edge[0],stiffener_edge_entities)

print(' ----  stiffener_surface_elements : ',stiffener_surface_elements)
print(' ----  stiffener_edge_elements : ',stiffener_edge_elements)


dataXfemStiffener = dict()
dataXfemStiffener['nodes'] = stiffener_nodes
dataXfemStiffener['stiffener_surface_elements'] = stiffener_surface_elements
dataXfemStiffener['stiffener_edge_elements'] = stiffener_edge_elements

gmsh.finalize()

# compute Level Set to stiffener
Stiffener_LS,Stiffener_distance = libF_xfem.computelevelset(datafluidmesh['nodes'],
                                                            dataXfemStiffener['nodes'],
                                                            dataXfemStiffener['stiffener_surface_elements'])
# make perpendicular mesh to the stiffener along the edge
Stiffener_tangent_nodes,Stiffener_tangent_mesh=libF_xfem.buildtangentedgemesh(dataXfemStiffener['nodes'],
                                                                         dataXfemStiffener['stiffener_surface_elements'],
                                                                         dataXfemStiffener['stiffener_edge_elements'])

# compute tangent Level Set to stiffener edge
Stiffener_tangent_LS,tmp = libF_xfem.computelevelset(datafluidmesh['nodes'],
                                                    Stiffener_tangent_nodes,
                                                    Stiffener_tangent_mesh)

if dataPb['flag_write_gmsh_results']==1:
    # gmsh output to check
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Stiffener_surface',
                                dataXfemStiffener['nodes'],
                                dataXfemStiffener['stiffener_surface_elements'],2)

    silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_LevelSet',
                                datafluidmesh['nodes'],
                                datafluidmesh['fluid_volume_elts'],
                                4,
                                [[[Stiffener_LS],'nodal',1,'Level set'],
                                 [[Stiffener_tangent_LS],'nodal',1,'Tangent Level set']])

# Get enriched nodes and elements
#LSEnrichedElementstmp,NbLSEnrichedElements=libF_xfem.getenrichedelementsfromlevelset(datafluidmesh['fluid_volume_elts'],
#                                                                                Stiffener_LS)
#
#LSEnrichedElements=LSEnrichedElementstmp[list(range(NbLSEnrichedElements))]

#EnrichedElementstmp1,NbEnrichedElements=libF_xfem.getsurfenrichedelements(dataXfemStiffener['nodes'],
#                                                                        dataXfemStiffener['stiffener_surface_elements'],
#                                                                        datafluidmesh['nodes'],
#                                                                        datafluidmesh['fluid_volume_elts'][LSEnrichedElements])
#EnrichedElementstmp2=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))])
#EnrichedElements=LSEnrichedElements[EnrichedElementstmp2-1]

# Get enriched nodes and elements directly from the stiffener surface mesh
EnrichedElementstmp1,NbEnrichedElements=libF_xfem.getsurfenrichedelements(dataXfemStiffener['nodes'],
                                                                        dataXfemStiffener['stiffener_surface_elements'],
                                                                        datafluidmesh['nodes'],
                                                                        datafluidmesh['fluid_volume_elts'])
EnrichedElements=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))]) # here, start with 1 (fortran indexing)


Enrichednodes = np.unique(datafluidmesh['fluid_volume_elts'][EnrichedElements-1])

# gmsh output to check
if dataPb['flag_write_gmsh_results']==1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Enriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][EnrichedElements-1],4)
    #silex_lib_gmsh.WriteResults(results_file.as_posix()+'_LSEnriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][LSEnrichedElements-1],4)

##############################################################
# Compute XFEM Fluid Matrices
##############################################################
#print(libF_xfem.computeedgeenrichment2.__doc__)
IIxf,JJxf,Vkaa,Vmaa,Vkfa,Vmfa = libF_xfem.computeedgeenrichment2(datafluidmesh['nodes'],
                                 datafluidmesh['fluid_volume_elts'],
                                 Stiffener_LS,
                                 Stiffener_tangent_LS,
                                 1.0,1.0)

HAA=scipy.sparse.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
HFA=scipy.sparse.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )

##############################################################
# Compute Standard Fluid Matrices : VOLUME
##############################################################

tic = time.process_time()
# get operators

IIf,JJf,Vffk,Vffm = libF.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                datafluidmesh['nodes'],
                                                1.0,
                                                1.0) # we put 1 for celerity and 1 for density


HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

##############################################################
# Compute Standard Fluid Matrices : FREE SURFACE
##############################################################



IIf,JJf,Vffktmp,VSFF=libFreeSurf.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                             datafluidmesh['nodes'][:,[0,1]],
                                                             1.0,1.0)

SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



##############################################################
# Compute Standard Fluid load : rigid body motion of tank
##############################################################

#print(libF.sloshimposedacc.__doc__)



CF,VecNormalEltsF = libF.sloshimposedacc(
                        np.array(datafluidmesh['nodes']),
                        np.array(datafluidmesh['structure_surface_elts']),
                        dataPb['U_dot_dot_imposed'])

CF=CF*dataFluid['rho']

# check if normal vectors are pointing out of the fluid volume
if dataPb['flag_write_gmsh_results']==1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces',
                            datafluidmesh['nodes'],
                            datafluidmesh['structure_surface_elts'],
                            2,
                            [[VecNormalEltsF,'elemental',3,'Normal to tank elements']])

##############################################################
# Compute XFEM Fluid load : rigid body motion of tank
##############################################################
print(libF_xfem.sloshimposedacc_xfem1.__doc__)

flag_write_quadrature_points_in_a_file=0

CA,VecNormalEltsA = libF_xfem.sloshimposedacc_xfem1(
                        np.array(datafluidmesh['nodes']),
                        np.array(dataXfemStiffener['nodes']),
                        np.array(datafluidmesh['fluid_volume_elts']),
                        np.array(dataXfemStiffener['stiffener_surface_elements']),
                        EnrichedElements,
                        dataPb['U_dot_dot_imposed'],
                        flag_write_quadrature_points_in_a_file)

CA=CA*dataFluid['rho']

# check if normal vectors are pointing out of the fluid volume
if dataPb['flag_write_gmsh_results']==1:
    silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_stiffener',
                            dataXfemStiffener['nodes'],
                            dataXfemStiffener['stiffener_surface_elements'],
                            2,
                            [[VecNormalEltsA,'elemental',3,'Normal to stiffener elements']])

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
# Assemble the whole system
##############################################################

SolvedDofF=list(range(fluid_ndof))
SolvedDofA=Enrichednodes-1

H=scipy.sparse.bmat( [ [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA]],
                                 [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA]]
                                 ] )
        
S=scipy.sparse.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],None],
                                 [None,HAA[SolvedDofA,:][:,SolvedDofA]*0.0]
                                 ] )
C = np.array([*CF[SolvedDofF], *CA[SolvedDofA]])
##############################################################
# Compute eigen modes
##############################################################

if dataPb['flag_eigen_vectors']==1:

    eigen_values,eigen_vectors= scipy.sparse.linalg.eigsh(H,
                                                              20,
                                                              S,
                                                              sigma=0,which='LM')

    freq_eigv_S=list(np.sqrt(eigen_values)/(2*np.pi))
    print('XFEM eigen frequencies : ',freq_eigv_S)
    eigen_vector_list=[]
    for i in range(eigen_values.shape[0]):
        Q=np.zeros(fluid_ndof)
        Q=eigen_vectors[SolvedDofF,i]
        eigen_vector_list.append(Q)


    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_Eigen_modes',
                                 datafluidmesh['nodes'],
                                 datafluidmesh['fluid_volume_elts'],
                                 4,
                                 [[eigen_vector_list,'nodal',1,'modes']])

##############################################################
# Compute FRF
##############################################################
if dataPb['flag_FRF']==1:
    press=[]
    frequencies=[]
    QuantityOfInterest=[]
    damping=None

    print ("Time at the beginning of the FRF: {}".format(time.ctime()))
    for f in np.linspace(dataPb['freq_ini'],
                        dataPb['freq_end'],
                        dataPb['nb_freq_step']):
        
        print('Solve freq : ',f)
        frequencies.append(f)

        omega=2*np.pi*f
        forceType ='float'

    #    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
        sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
        press.append(sol.copy())

        QuantityOfInterest.append(sol[8-1]) # upper corner
        


    frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

    if dataPb['flag_write_gmsh_results']==1:
        print('QuantityOfInterest : ',QuantityOfInterest)
        silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf',
                                datafluidmesh['nodes'],
                                datafluidmesh['fluid_volume_elts'],
                                4,
                                [[press,'nodal',1,'pressure']]
                                )

    print ("Time at the end of the FRF: {}".format(time.ctime()))

    f=open(results_file.as_posix() +'_results.frf','wb')
    pickle.dump(frfsave, f)
    f.close()

