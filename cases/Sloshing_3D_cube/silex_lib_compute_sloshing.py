import string
import time
import numpy as np
import scipy
import scipy.sparse as sps
import scipy.sparse.linalg  as spla
from loguru import logger

import pylab as pl
import pickle
import csv

import sys
from pathlib import Path

import pymumps as mumps
import gmsh
# import utils as u
# import utils_acoustics as ua
from meshRW import msh2

# useful tools
from SILEXrun import utils as misc_utils

# for tet10
from SILEXlib import silex_lib_acou_tet10 as libF_tet10
from SILEXlib import silex_lib_xfem_acou_tet10 as libF_tet10_xfem
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_tet4_xfem
from SILEXlib import silex_lib_acou_tri6 as libFreeSurf_tet10
# from SILEXlib import silex_lib_xfem_acou_tri6 as libFreeSurf_tet10_xfem

# for tet4
from SILEXlib import silex_lib_acou_tet4 as libF_tet4
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_tet4_xfem
from SILEXlib import silex_lib_acou_tri3 as libFreeSurf_tet4

# for tet4 and tet10 : level set
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_levelset

# for post-processing XFEM results
from SILEXlib import MeshField as lib

# tools
from SILEXlib.utils import utils


# for DKT
from SILEXlib import silex_lib_dkt as libDKT

#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_gmsh

import silex_lib_cube_tank_gmsh_geometry as mesher_cube_tank

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

def solve_linear(method, A, b, comm=None):
    if method == 'mumps':
        x = mumps.spsolve(A, b, comm=comm)
    elif method == 'scipy':
        x = spla.spsolve(A, b)
    return x


def sloshing_rigid_baffle_tet10_xfem(dataPb,dataFluid,mesh_file_fluid,mesh_file_stiffener,results_file):
    
    ##############################################################
    # Load fluid mesh and tank
    ##############################################################

    tic = time.process_time()
    
    #silex_lib_gmsh
    
    #gmsh.initialize()
    #gmsh.open(mesh_file_fluid.as_posix()+'.geo_unrolled')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #fluid_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))

    fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_fluid.as_posix()+'.msh',3)
    fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',11,10)
    
    structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,20)
    free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,30)

    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
    #structure_surface = (2,1) # structure surface physical group : rigid tank surface phys. 1
    #free_fluid_surface = (2,2) # free fluid surface physical group : phys. 2
    #fluid_volume = (3,3) # fluid physical group : phys 3
    #if not structure_surface in physical_groups:
    #    raise ValueError('Structure surface physical group not found')
    #if not fluid_volume in physical_groups:
    #    raise ValueError('Fluid volume physical group not found')
    #if not free_fluid_surface in physical_groups:
    #    raise ValueError('Free fluid surface physical group not found')
#
#
    ## get entities per physical group
    #structure_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*structure_surface)
    #free_fluid_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*free_fluid_surface)
    #fluid_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*fluid_volume)
#
    ## get elements
    #structure_surface_elements = u.getElementsFromEntities(gmsh,structure_surface[0],structure_surface_entities)
    #free_fluid_surface_elements = u.getElementsFromEntities(gmsh,free_fluid_surface[0],free_fluid_surface_entities)
    #fluid_volume_elements = u.getElementsFromEntities(gmsh,fluid_volume[0],fluid_volume_entities)

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

    #gmsh.finalize()

    fluid_ndof = datafluidmesh['nodes'].shape[0]

    print('node 8 :', datafluidmesh['nodes'][8-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],11)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],9)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],9)

    ##############################################################
    # Load stiffener mesh / compute Level Set / Get enriched nodes and elements
    ##############################################################

    tic = time.process_time()

    stiffener_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_stiffener.as_posix()+'.msh',3)
    stiffener_surface_elements,Idstiffener_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',2,50)
    stiffener_edge_elements,Idnode_stiffener_edge_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',1,60)


    #gmsh.initialize()
    #gmsh.open(mesh_file_stiffener.as_posix()+'.geo')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #stiffener_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))
#
    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
#
    #stiffener_surface = (2,50) # stiffener surface physical group
    #stiffener_edge = (1,60) # stiffener edge physical group
#
    #if not stiffener_surface in physical_groups:
    #    raise ValueError('stiffener surface physical group not found')
    #if not stiffener_edge in physical_groups:
    #    raise ValueError('stiffener edge physical group not found')
#
    ## get entities per physical group
    #stiffener_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_surface)
    #stiffener_edge_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_edge)
#
    ## get elements
    #stiffener_surface_elements = u.getElementsFromEntities(gmsh,stiffener_surface[0],stiffener_surface_entities)
    #stiffener_edge_elements = u.getElementsFromEntities(gmsh,stiffener_edge[0],stiffener_edge_entities)
#
    #print(' ----  stiffener_surface_elements : ',stiffener_surface_elements)
    #print(' ----  stiffener_edge_elements : ',stiffener_edge_elements)


    dataXfemStiffener = dict()
    dataXfemStiffener['nodes'] = stiffener_nodes
    dataXfemStiffener['stiffener_surface_elements'] = stiffener_surface_elements
    dataXfemStiffener['stiffener_edge_elements'] = stiffener_edge_elements

    #gmsh.finalize()

    # compute Level Set to stiffener
    Stiffener_LS,Stiffener_distance = libF_levelset.computelevelset(datafluidmesh['nodes'],
                                                                dataXfemStiffener['nodes'],
                                                                dataXfemStiffener['stiffener_surface_elements'])
    # make perpendicular mesh to the stiffener along the edge
    Stiffener_tangent_nodes,Stiffener_tangent_mesh=libF_levelset.buildtangentedgemesh(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            dataXfemStiffener['stiffener_edge_elements'])

    # compute tangent Level Set to stiffener edge
    Stiffener_tangent_LS,tmp = libF_levelset.computelevelset(datafluidmesh['nodes'],
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
                                    11,
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
    EnrichedElementstmp1,NbEnrichedElements=libF_levelset.getsurfenrichedelements(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            datafluidmesh['nodes'],
                                                                            datafluidmesh['fluid_volume_elts'][:,0:4])
    EnrichedElements=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))]) # here, start with 1 (fortran indexing)


    Enrichednodes = np.unique(datafluidmesh['fluid_volume_elts'][EnrichedElements-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Enriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][EnrichedElements-1],11)
        #silex_lib_gmsh.WriteResults(results_file.as_posix()+'_LSEnriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][LSEnrichedElements-1],4)



    ##############################################################
    # Compute XFEM Fluid Matrices
    ##############################################################
    #print(libF_xfem.globalxfemacousticmatrices.__doc__)

    IIxf,JJxf,Vkaa,_,Vkfa,_ = libF_tet10_xfem.globalxfemacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                                        datafluidmesh['nodes'],
                                                                        Stiffener_LS,
                                                                        Stiffener_tangent_LS*0.0-1.0,
                                                                        1.0,1.0)

    HAA=sps.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=sps.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )


    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,_ = libF_tet10.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

    ##############################################################
    # Compute XFEM Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

    ##############################################################
    # Compute Standard Fluid load : rigid body motion of tank
    ##############################################################

    #print(libF.sloshimposedacc.__doc__)



    CF,VecNormalEltsF = libF_tet10.sloshimposedacc(
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
    #print(libF_xfem.sloshimposedacc_xfem1.__doc__)

    flag_write_quadrature_points_in_a_file=0

    CA,VecNormalEltsA = libF_tet10_xfem.sloshimposedacc_xfem1(
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

    H=sps.construct.bmat( [ [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA]],
                                    [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )
    maxS = np.abs(SFF.max())    
    S=sps.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],None],
                                     [None,np.zeros((len(SolvedDofA),len(SolvedDofA)))]] )
    # matdense= SFF[SolvedDofF,:][:,SolvedDofF].todense()
    # matdenseAA = HAA[SolvedDofA,:][:,SolvedDofA].todense()*0.0
    # matdenseFA = HFA[SolvedDofF,:][:,SolvedDofA].todense()*0.0
    # S=sps.construct.bmat( [ [matdense+maxS/100,matdenseFA+maxS/100],
    #                                  [matdenseFA.transpose()+maxS/100,matdenseAA+maxS/100]])
    

    ##############################################################
    # Build a tet4 mesh from tet10 fluid mesh : just for plotting
    ##############################################################

    elemtet4=libF_tet10_xfem.tet10totet4(np.array(datafluidmesh['nodes']),
                                         np.array(datafluidmesh['fluid_volume_elts']))
            
    # Get enriched nodes and elements directly from the stiffener surface mesh
    EnrichedElementstet4tmp1,NbEnrichedElementstet4=libF_levelset.getsurfenrichedelements(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            datafluidmesh['nodes'],
                                                                            elemtet4)
    EnrichedElementstet4=np.unique(EnrichedElementstet4tmp1[list(range(NbEnrichedElementstet4))])  # here, start with 1 (fortran indexing)
    Enrichednodestet4 = np.unique(elemtet4[EnrichedElementstet4-1])
    # nodes enriched for tet10 / but no enriched for tet4
    SpecialTet4nodes = np.setdiff1d(Enrichednodes,Enrichednodestet4) 
    


    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume_tet10TOtet4',
                                    datafluidmesh['nodes'],
                                    elemtet4,4)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Enriched_Fluid_Elements_tet10TOtet4',
                                    datafluidmesh['nodes'],
                                    elemtet4[EnrichedElementstet4-1],4)





    ##############################################################
    # Compute eigen modes
    ##############################################################

    if dataPb['flag_eigen_vectors']==1:

        eigen_values,eigen_vectors= spla.eigsh(H,
                                                                dataPb['flag_nb_eigen_modes'],
                                                                S,
                                                                sigma=0,which='LM')

        freq_eigv_S=list(np.sqrt(eigen_values)/(2*np.pi))
        print('XFEM eigen frequencies : ',freq_eigv_S)
        eigen_vector_list=[]
        presslist=[]
        pressenrichlist=[]
        for i in range(eigen_values.shape[0]):
            #Q=np.zeros(fluid_ndof)
            #Q=eigen_vectors[SolvedDofF,i]
            press = np.zeros(fluid_ndof)
            Q = eigen_vectors[:,i].copy()
            press[SolvedDofF] = Q[SolvedDofF]
            enrichment = np.zeros(fluid_ndof)
            enrichment[SolvedDofA]= Q[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))]
            
            presslist.append(press)
            pressenrichlist.append(enrichment)
            
            CorrectedPressure=np.array(press)
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA].T+np.array(enrichment[SolvedDofA]*np.sign(Stiffener_LS[SolvedDofA]).T)
            eigen_vector_list.append(CorrectedPressure)
            #eigen_vector_list.append(Q)

        f=open(results_file.as_posix() +'_eigen_frequencies.pck','wb')
        pickle.dump(freq_eigv_S, f)
        f.close()





        if dataPb['flag_write_gmsh_results']==1:
            silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_Eigen_modes',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    11,
                                    [[eigen_vector_list,'nodal',1,'modes'],
                                     [presslist,'nodal',1,'press classic'],
                                     [pressenrichlist,'nodal',1,'press enrich']])
            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_eigenmodes_on_tet4mesh',
                                    datafluidmesh['nodes'],
                                    elemtet4,
                                    4,
                                    [[eigen_vector_list,'nodal',1,'pressure']]
                                    )

            objMesh = lib.MeshField(nodes=datafluidmesh['nodes'], 
                                    elems=elemtet4, 
                                    leveset=Stiffener_LS, 
                                    levelsetTg=Stiffener_tangent_LS)
            objMesh.addField(uncorrectedField=np.vstack(presslist).transpose(),
                             enrichmentField=np.vstack(pressenrichlist).transpose())
            datameshfield = objMesh.getData()
            
            # prepare fields
            dataW = []
            dataW.append({'data':datameshfield['LS'],'type':'nodal', 'name':'levelset'})
            dataW.append({'data':datameshfield['LST'],'type':'nodal','name':'tangent levelset'})
            # for it in range(len(press)):
            #     dataW.append({'data':datameshfield['fields'][:,it],'type':'nodal','name':'field '+str(it)+' (levelset)'})
            dataW = {
                'name': 'press',
                'nbsteps': len(presslist),
                'type': 'nodal',
                'data': datameshfield['fields']
            }
#            # export mesh
            msh2.mshWriter(
                filename= results_file.as_posix() +'_results_fluid_eigenmodes_meshfield3D.msh',
                nodes=datameshfield['nodes'],
                elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
                        {'type':'PRI6','connectivity':datameshfield['PRI6']}],
                fields=dataW,
                append=True
                )

    ##############################################################
    # Compute FRF
    ##############################################################
    if dataPb['flag_FRF']==1:
        Correctedpress=[]
        press=[]
        press_no_baffle=[]
        frequencies=[]
        enrichpress=[]
        QuantityOfInterest=[]
        damping=None

        print ("Time at the beginning of the FRF:",time.ctime())
        for f in np.linspace(dataPb['freq_ini'],
                            dataPb['freq_end'],
                            dataPb['nb_freq_step']):
            
            print('Solve freq : ',f)
            frequencies.append(f)

            omega=2*np.pi*f
            #forceType ='float'

        #    sol = spla.spsolve(KFF-omega**2*MFF,FF)
            C = np.array([*CF[SolvedDofF]*(-omega**2), *CA[SolvedDofA]*(-omega**2)])
            sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
            # sol = spla.spsolve(H-omega**2*S, C)
            
            CorrectedPressure = np.zeros(fluid_ndof)
            UncorrectedPressure = np.zeros(fluid_ndof)
            CorrectedPressure[SolvedDofF] = sol[SolvedDofF].copy()
            UncorrectedPressure[SolvedDofF] = sol[SolvedDofF].copy()
            enrichment = np.zeros(fluid_ndof)
            enrichment[SolvedDofA]= sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))].copy()
            enrichpress.append(enrichment)
            #CorrectedPressure=np.array(pressure)
            CorrectedPressure=CorrectedPressure+np.array(enrichment*np.sign(Stiffener_LS).T)
            Correctedpress.append(CorrectedPressure)
            
            # sol = spla.spsolve(HFF-omega**2*SFF, CF[SolvedDofF]*(-omega**2))
            # sol = mumps.spsolve( HFF-omega**2*SFF, CF[SolvedDofF]*(-omega**2), comm=mycomm )
            # press_no_baffle.append(sol)
            # soltmp = np.zeros(fluid_ndof) #np.zeros_like(sol[SolvedDofF])
            # soltmp[SolvedDofA] = sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))].copy()
            # soltmp[SpecialTet4nodes-1] = 0.0
            # enrichpress.append(soltmp)
            # pressSpecial= np.zeros(fluid_ndof)
            # pressSpecial[SolvedDofF]=sol[SolvedDofF].copy()
            # pressSpecial[SpecialTet4nodes-1]=pressSpecial[SpecialTet4nodes-1]+np.array(enrichment[SpecialTet4nodes-1]*np.sign(Stiffener_LS[SpecialTet4nodes-1]).T)
            # press.append(pressSpecial)
            #enrichpress.append(CorrectedPressure*0.0)
            press.append(UncorrectedPressure)

            QuantityOfInterest.append(sol[8-1]) # upper corner
            
            


        frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

        if dataPb['flag_write_gmsh_results']==1:
            print('QuantityOfInterest : ',QuantityOfInterest)

            print(libF_tet10_xfem.tet10totet4.__doc__)

            objMesh = lib.MeshField(nodes=datafluidmesh['nodes'], 
                                    elems=datafluidmesh['fluid_volume_elts'], 
                                    leveset=Stiffener_LS, 
                                    levelsetTg=Stiffener_tangent_LS)
            objMesh.addField(uncorrectedField=np.vstack(press).transpose(),
                             enrichmentField=np.vstack(enrichpress).transpose())
            datameshfield = objMesh.getData()
            
            # prepare fields
            dataW = []
            dataW.append({'data':datameshfield['LS'],'type':'nodal', 'name':'levelset'})
            dataW.append({'data':datameshfield['LST'],'type':'nodal','name':'tangent levelset'})
            dataW.append({'data':datameshfield['signLS'],'type':'nodal','name':'sign levelset'})
            # for it in range(len(press)):
            #     dataW.append({'data':datameshfield['fields'][:,it],'type':'nodal','name':'field '+str(it)+' (levelset)'})
            dataW.append({
                'name': 'press',
                'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['fields']
            })
            dataW.append({
                'name': 'uncorrected press',
                 'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['uncorrected']
            })
            dataW.append({
                'name': 'correction press',
                 'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['correction']
            })
            
            
            msh2.mshWriter(
                filename= results_file.as_posix() +'_results_fluid_frf_init.msh',
                nodes=datameshfield['nodes'],
                elements=[{'type':'TET10','connectivity':datafluidmesh['fluid_volume_elts'][2:5,:]}],
                fields=dataW,
                append=True
                )
#            # export mesh
            msh2.mshWriter(
                filename= results_file.as_posix() +'_results_fluid_frf_meshfield3D.msh',
                nodes=datameshfield['nodes'],
                elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
                          {'type':'PRI6','connectivity':datameshfield['PRI6']}],
                fields=dataW,
                append=True
                )

            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    11,
                                    [[Correctedpress,'nodal',1,'pressure'],
                                     [press,'nodal',1,'uncorrectpressure'],
                                     [enrichpress,'nodal',1,'correction'],
                                     [press_no_baffle,'nodal',1,'pressure_no_baffle']]
                                    )
            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf_on_tet4mesh',
                                    datafluidmesh['nodes'],
                                    elemtet4,
                                    4,
                                    [[Correctedpress,'nodal',1,'pressure'],
                                     [press,'nodal',1,'uncorrectpressure'],
                                     [enrichpress,'nodal',1,'correction']]
                                    )

        print ("Time at the end of the FRF:",time.ctime())

        f=open(results_file.as_posix() +'_results.frf','wb')
        pickle.dump(frfsave, f)
        f.close()
    return

def sloshing_rigid_baffle_tet4_xfem(dataPb,dataFluid,mesh_file_fluid,mesh_file_stiffener,results_file, convert_from_tet10=False):
    ##############################################################
    # Load fluid mesh and tank
    ##############################################################

    tic = time.process_time()

    if convert_from_tet10:
        fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_fluid.as_posix()+'.msh',3)
        fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',11,10)
        
        structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,20)
        free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,30)
        # convert tet10 to tet4
        fluid_volume_elements = fluid_volume_elements[:,0:4] # just take the first 4 nodes of the tet10 elements
        structure_surface_elements = structure_surface_elements[:,0:3]
        free_fluid_surface_elements = free_fluid_surface_elements[:,0:3]
        # remove unsed nodes
        nb_nodes_initial = fluid_nodes.shape[0]
        fluid_nodes = fluid_nodes[np.unique(fluid_volume_elements)-1,:] # -1
        # renumber elements
        renumbering = np.zeros(nb_nodes_initial)
        renumbering[np.unique(fluid_volume_elements)-1] = np.arange(1,fluid_nodes.shape[0]+1)
        fluid_volume_elements = renumbering[fluid_volume_elements-1].astype(int) # -1 for fortran indexing
        structure_surface_elements = renumbering[structure_surface_elements-1].astype(int) # -1 for fortran indexing
        free_fluid_surface_elements = renumbering[free_fluid_surface_elements-1].astype(int) # -1 for fortran indexing
        # update
        
    else:
        fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_fluid.as_posix()+'.msh',3)
        fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',4,10)
        
        structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',2,20)
        free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',2,30)



    #gmsh.initialize()
    #gmsh.open(mesh_file_fluid.as_posix()+'.geo')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #fluid_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))
#
    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
    #structure_surface = (2,20) # structure surface physical group
    #free_fluid_surface = (2,30) # free fluid surface physical group
    #fluid_volume = (3,10) # fluid physical group
    #if not structure_surface in physical_groups:
    #    raise ValueError('Structure surface physical group not found')
    #if not fluid_volume in physical_groups:
    #    raise ValueError('Fluid volume physical group not found')
    #if not free_fluid_surface in physical_groups:
    #    raise ValueError('Free fluid surface physical group not found')
#
#
    ## get entities per physical group
    #structure_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*structure_surface)
    #free_fluid_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*free_fluid_surface)
    #fluid_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*fluid_volume)
#
    ## get elements
    #structure_surface_elements = u.getElementsFromEntities(gmsh,structure_surface[0],structure_surface_entities)
    #free_fluid_surface_elements = u.getElementsFromEntities(gmsh,free_fluid_surface[0],free_fluid_surface_entities)
    #fluid_volume_elements = u.getElementsFromEntities(gmsh,fluid_volume[0],fluid_volume_entities)

    # print(' ----  structure_surface_elements : ',structure_surface_elements)
    # print(' ----  free_fluid_surface_elements : ',free_fluid_surface_elements)
    # print(' ----  fluid_volume_elements : ',fluid_volume_elements)

    datafluidmesh = dict()
    datafluidmesh['nodes'] = fluid_nodes
    datafluidmesh['fluid_volume_elts'] = fluid_volume_elements
    datafluidmesh['structure_surface_elts'] = structure_surface_elements
    datafluidmesh['free_fluid_surface_elts'] = free_fluid_surface_elements
    #datastructuremesh = dict()
    #datastructuremesh['structure_surface_elts'] = structure_surface_elements
    #datastructuremesh['nodes'] = fluid_nodes # !!!! attention aux tailles des matrices par la suite !!!!!

    #gmsh.finalize()

    fluid_ndof = datafluidmesh['nodes'].shape[0]

    print('node 8 :', datafluidmesh['nodes'][8-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume_tet4',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],4)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces_tet4',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],2)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface_tet4',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],2)

    ##############################################################
    # Load stiffener mesh / compute Level Set / Get enriched nodes and elements
    ##############################################################

    tic = time.process_time()

    stiffener_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_stiffener.as_posix()+'.msh',3)
    stiffener_surface_elements,Idstiffener_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',2,50)
    stiffener_edge_elements,Idnode_stiffener_edge_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',1,60)



    #gmsh.initialize()
    #gmsh.open(mesh_file_stiffener.as_posix()+'.geo')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #stiffener_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))
#
    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
#
    #stiffener_surface = (2,50) # stiffener surface physical group
    #stiffener_edge = (1,60) # stiffener edge physical group
#
    #if not stiffener_surface in physical_groups:
    #    raise ValueError('stiffener surface physical group not found')
    #if not stiffener_edge in physical_groups:
    #    raise ValueError('stiffener edge physical group not found')
#
    ## get entities per physical group
    #stiffener_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_surface)
    #stiffener_edge_entities = gmsh.model.getEntitiesForPhysicalGroup(*stiffener_edge)
#
    ## get elements
    #stiffener_surface_elements = u.getElementsFromEntities(gmsh,stiffener_surface[0],stiffener_surface_entities)
    #stiffener_edge_elements = u.getElementsFromEntities(gmsh,stiffener_edge[0],stiffener_edge_entities)

    # print(' ----  stiffener_surface_elements : ',stiffener_surface_elements)
    # print(' ----  stiffener_edge_elements : ',stiffener_edge_elements)


    dataXfemStiffener = dict()
    dataXfemStiffener['nodes'] = stiffener_nodes
    dataXfemStiffener['stiffener_surface_elements'] = stiffener_surface_elements
    dataXfemStiffener['stiffener_edge_elements'] = stiffener_edge_elements

    #gmsh.finalize()

    # compute Level Set to stiffener
    Stiffener_LS,Stiffener_distance = libF_levelset.computelevelset(datafluidmesh['nodes'],
                                                                dataXfemStiffener['nodes'],
                                                                dataXfemStiffener['stiffener_surface_elements'])
    # make perpendicular mesh to the stiffener along the edge
    Stiffener_tangent_nodes,Stiffener_tangent_mesh=libF_levelset.buildtangentedgemesh(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            dataXfemStiffener['stiffener_edge_elements'])

    # compute tangent Level Set to stiffener edge
    Stiffener_tangent_LS,tmp = libF_levelset.computelevelset(datafluidmesh['nodes'],
                                                        Stiffener_tangent_nodes,
                                                        Stiffener_tangent_mesh)

    # Stiffener_tangent_LS = Stiffener_tangent_LS*0.0-1.0 # we put -1 to have the same sign as the Level Set

    if dataPb['flag_write_gmsh_results']==1:
        # gmsh output to check
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Stiffener_surface_tet4',
                                    dataXfemStiffener['nodes'],
                                    dataXfemStiffener['stiffener_surface_elements'],2)

        silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_LevelSet_tet4',
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
    print('ici 1')
    EnrichedElementstmp1,NbEnrichedElements=libF_levelset.getsurfenrichedelements(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            datafluidmesh['nodes'],
                                                                            datafluidmesh['fluid_volume_elts'])
    EnrichedElements=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))]) # here, start with 1 (fortran indexing)


    Enrichednodes = np.unique(datafluidmesh['fluid_volume_elts'][EnrichedElements-1])

    print('ici 2')

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Enriched_Fluid_Elements_tet4',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][EnrichedElements-1],4)
        #silex_lib_gmsh.WriteResults(results_file.as_posix()+'_LSEnriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][LSEnrichedElements-1],4)

    ##############################################################
    # Compute XFEM Fluid Matrices
    ##############################################################
    #print(libF_xfem.computeedgeenrichment2.__doc__)
    IIxf,JJxf,Vkaa,Vmaa,Vkfa,Vmfa = libF_tet4_xfem.computeedgeenrichment2(datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    Stiffener_LS,
                                    Stiffener_tangent_LS,
                                    1.0,1.0)

    HAA=sps.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=sps.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,Vffm = libF_tet4.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet4.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



    ##############################################################
    # Compute Standard Fluid load : rigid body motion of tank
    ##############################################################

    #print(libF.sloshimposedacc.__doc__)



    CF,VecNormalEltsF = libF_tet4.sloshimposedacc(
                            np.array(datafluidmesh['nodes']),
                            np.array(datafluidmesh['structure_surface_elts']),
                            dataPb['U_dot_dot_imposed'])

    CF=CF*dataFluid['rho']

    # check if normal vectors are pointing out of the fluid volume
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces_tet4',
                                datafluidmesh['nodes'],
                                datafluidmesh['structure_surface_elts'],
                                2,
                                [[VecNormalEltsF,'elemental',3,'Normal to tank elements']])

    ##############################################################
    # Compute XFEM Fluid load : rigid body motion of tank
    ##############################################################
    #print(libF_xfem.sloshimposedacc_xfem1.__doc__)

    flag_write_quadrature_points_in_a_file=0

    CA,VecNormalEltsA = libF_tet4_xfem.sloshimposedacc_xfem1(
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
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_stiffener_tet4',
                                dataXfemStiffener['nodes'],
                                dataXfemStiffener['stiffener_surface_elements'],
                                2,
                                [[VecNormalEltsA,'elemental',3,'Normal to stiffener elements']])

    ##############################################################
    # Compute coupling terms around the stiffener
    ##############################################################
    IIc1,JJc1,Vc1=libF_tet4_xfem.computexfemcoupling1(datafluidmesh['nodes'],
                                                      dataXfemStiffener['nodes'],
                                                      datafluidmesh['fluid_volume_elts'],
                                                      dataXfemStiffener['stiffener_surface_elements'],
                                                      EnrichedElements)
    #IIc2,JJc2,Vc2=libF_tet10_xfem.computexfemcoupling2(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements,LevelSet)
    stiffener_ndof = dataXfemStiffener['nodes'].shape[0]*6
    CSA=sps.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(stiffener_ndof,fluid_ndof) ) 
    
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

    H=sps.construct.bmat( [ [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA]],
                                    [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )
            
    S=sps.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],None],
                                    [None,HAA[SolvedDofA,:][:,SolvedDofA]*0.0]
                                    ] )
    #C = np.array([*CF[SolvedDofF], *CA[SolvedDofA]])
    ##############################################################
    # Compute eigen modes
    ##############################################################

    if dataPb['flag_eigen_vectors']==1:

        eigen_values,eigen_vectors= spla.eigsh(H,
                                                                dataPb['flag_nb_eigen_modes'],
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
            silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_Eigen_modes_tet4',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    4,
                                    [[eigen_vector_list,'nodal',1,'modes']])

    ##############################################################
    # Compute FRF
    ##############################################################
    meanQI=0
    maxQI=0
    meanforce=0
    maxforce=0
    if dataPb['flag_FRF']==1:
        Correctedpress=[]
        press=[]
        force = []
        enrichpress = []
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

        #    sol = spla.spsolve(KFF-omega**2*MFF,FF)
            C = np.array([*CF[SolvedDofF]*(-omega**2), *CA[SolvedDofA]*(-omega**2)])
            sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
            
            CorrectedPressure = np.zeros(fluid_ndof)
            press.append(sol[SolvedDofF].copy())
            CorrectedPressure[SolvedDofF] = sol[SolvedDofF].copy()
            soltmp = np.zeros_like(sol[SolvedDofF])
            soltmp[SolvedDofA] = sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))].copy()

            enrichpress.append(soltmp)
            torseur = CSA@soltmp
            
            fxnodes=torseur[0::6]
            fynodes=torseur[1::6]
            fznodes=torseur[2::6]
            for kk in range(100):
                fxnodes[np.argmax(fxnodes)]=0.0
                fxnodes[np.argmin(fxnodes)]=0.0
                fynodes[np.argmax(fynodes)]=0.0
                fynodes[np.argmin(fynodes)]=0.0            
                fznodes[np.argmax(fznodes)]=0.0
                fznodes[np.argmin(fznodes)]=0.0                
            #print(max(fxnodes))
            #titi=np.argmax(fxnodes)
            #print(titi)
            #print(fxnodes[titi])
            fx = np.sum(fxnodes)
            fy = np.sum(fynodes)
            fz = np.sum(fznodes)
            #
            force.append(np.linalg.norm(np.array([fx,fy,fz])))
            
            CorrectedPressure=CorrectedPressure+np.array(soltmp*np.sign(Stiffener_LS).T)
            Correctedpress.append(CorrectedPressure)
            
            QuantityOfInterest.append(sol[8-1]) # upper corner
        maxQI = np.max(np.abs(QuantityOfInterest))
        meanQI = np.mean(np.abs(QuantityOfInterest))
        maxforce = np.max(force)
        meanforce = np.mean(force)
            


        frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

        if dataPb['flag_write_gmsh_results']==1:
            print('QuantityOfInterest : ',QuantityOfInterest)
            objMesh = lib.MeshField(nodes=datafluidmesh['nodes'], 
                                    elems=datafluidmesh['fluid_volume_elts'], 
                                    leveset=Stiffener_LS, 
                                    levelsetTg=Stiffener_tangent_LS)
            objMesh.addField(uncorrectedField=np.vstack(press).transpose(),
                             enrichmentField=np.vstack(enrichpress).transpose())
            datameshfield = objMesh.getData()
            # prepare fields
            dataW = []
            dataW.append({'data':datameshfield['LS'],'type':'nodal', 'name':'levelset'})
            dataW.append({'data':datameshfield['LST'],'type':'nodal','name':'tangent levelset'})
            # for it in range(len(press)):
            #     dataW.append({'data':datameshfield['fields'][:,it],'type':'nodal','name':'field '+str(it)+' (levelset)'})
            dataW.append({
                'name': 'press',
                'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['fields']
            })
            dataW.append({
                'name': 'uncorrected press',
                 'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['uncorrected']
            })
            dataW.append({
                'name': 'correction press',
                 'nbsteps': len(press),
                'type': 'nodal',
                'data': datameshfield['correction']
            })
            # export mesh
            msh2.mshWriter(
                filename= results_file.as_posix() +'_results_fluid_frf_tet4_meshfield3D.msh',
                nodes=datameshfield['nodes'],
                elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
                        {'type':'PRI6','connectivity':datameshfield['PRI6']}],
                fields=dataW,
                append=True
                )

            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf_tet4',
                                    datafluidmesh['nodes'],
                                    fluid_volume_elements,
                                    4,
                                    [[Correctedpress,'nodal',1,'pressure'],
                                     [press,'nodal',1,'uncorrectpressure'],
                                     [enrichpress,'nodal',1,'correction']]
                                    )

        print ("Time at the end of the FRF:",time.ctime())

        f=open(results_file.as_posix() +'_results.frf','wb')
        pickle.dump(frfsave, f)
        f.close()
    return meanQI, maxQI, meanforce, maxforce

def sloshing_rigid_baffle_tet4(dataPb,dataFluid,mesh_file,results_file):
    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()
    
    results_file.parent.mkdir(parents=True, exist_ok=True)

    fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file.as_posix()+'.msh',3)
    fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',4,10)
    
    structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',2,20)
    free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',2,30)

    #gmsh.initialize()
    #gmsh.open(mesh_file.as_posix()+'.geo')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #fluid_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))
#
    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
    #structure_surface = (2,20) # structure surface physical group
    #free_fluid_surface = (2,30) # free fluid surface physical group
    #fluid_volume = (3,10) # fluid physical group
    #if not structure_surface in physical_groups:
    #    raise ValueError('Structure surface physical group not found')
    #if not fluid_volume in physical_groups:
    #    raise ValueError('Fluid volume physical group not found')
    #if not free_fluid_surface in physical_groups:
    #    raise ValueError('Free fluid surface physical group not found')
#
#
    ## get entities per physical group
    #structure_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*structure_surface)
    #free_fluid_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*free_fluid_surface)
    #fluid_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*fluid_volume)
#
    ## get elements
    #structure_surface_elements = u.getElementsFromEntities(gmsh,structure_surface[0],structure_surface_entities)
    #free_fluid_surface_elements = u.getElementsFromEntities(gmsh,free_fluid_surface[0],free_fluid_surface_entities)
    #fluid_volume_elements = u.getElementsFromEntities(gmsh,fluid_volume[0],fluid_volume_entities)

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

    #gmsh.finalize()

    print('node 8 :', datafluidmesh['nodes'][8-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results'] == 1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],4)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],2)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],2)



    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,Vffm = libF_tet4.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density

    fluid_ndof = datafluidmesh['nodes'].shape[0]

    HFF=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet4.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



    ##############################################################
    # Compute Standard Fluid load : rigid body motion of tank
    ##############################################################

    #print(libF.sloshimposedacc.__doc__)

    CF,VecNormalElts = libF_tet4.sloshimposedacc(
        np.array(datafluidmesh['nodes']),
                                    np.array(datafluidmesh['structure_surface_elts']),
                                    dataPb['U_dot_dot_imposed'])

    CF=CF*dataFluid['rho']

    # check if normal vectors are pointing out of the fluid volume
    if dataPb['flag_write_gmsh_results'] == 1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['structure_surface_elts'],
                                    2,
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

        eigen_values,eigen_vectors= spla.eigsh(HFF,
                                                                dataPb['flag_nb_eigen_modes'],
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
                                        4,
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

        #    sol = spla.spsolve(KFF-omega**2*MFF,FF)
            sol = mumps.spsolve( HFF-omega**2*SFF , CF*(-omega**2) , comm=mycomm )
            press.append(sol.copy())

            QuantityOfInterest.append(sol[8-1]) # upper corner



        frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

        if dataPb['flag_write_gmsh_results'] == 1:
            print('QuantityOfInterest : ',QuantityOfInterest)
            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf',
                                        datafluidmesh['nodes'],
                                        datafluidmesh['fluid_volume_elts'],
                                        4,
                                        [[press,'nodal',1,'pressure']]
                                        )

        print ("Time at the end of the FRF:",time.ctime())

        f=open(results_file.as_posix() +'_results.frf','wb')
        pickle.dump(frfsave, f)
        f.close()

def sloshing_rigid_baffle_tet10(dataPb,dataFluid,mesh_file,results_file):
    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()
    fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file.as_posix()+'.msh',3)
    fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',11,10)
    
    structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',9,20)
    free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file.as_posix()+'.msh',9,30)

    #gmsh.initialize()
    #gmsh.open(mesh_file.as_posix()+'.geo')
    ## mesh generation
    #gmsh.model.mesh.generate(3)
#
    ## get fluid nodes
    #_,_nodes,_ = gmsh.model.mesh.getNodes()
    #fluid_nodes = np.reshape(_nodes,(int(len(_nodes)/3),3))
#
    ## check physical groups
    #physical_groups = gmsh.model.getPhysicalGroups()
    #structure_surface = (2,20) # structure surface physical group
    #free_fluid_surface = (2,30) # free fluid surface physical group
    #fluid_volume = (3,10) # fluid physical group
    #if not structure_surface in physical_groups:
    #    raise ValueError('Structure surface physical group not found')
    #if not fluid_volume in physical_groups:
    #    raise ValueError('Fluid volume physical group not found')
    #if not free_fluid_surface in physical_groups:
    #    raise ValueError('Free fluid surface physical group not found')
#
#
    ## get entities per physical group
    #structure_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*structure_surface)
    #free_fluid_surface_entities = gmsh.model.getEntitiesForPhysicalGroup(*free_fluid_surface)
    #fluid_volume_entities = gmsh.model.getEntitiesForPhysicalGroup(*fluid_volume)
#
    ## get elements
    #structure_surface_elements = u.getElementsFromEntities(gmsh,structure_surface[0],structure_surface_entities)
    #free_fluid_surface_elements = u.getElementsFromEntities(gmsh,free_fluid_surface[0],free_fluid_surface_entities)
    #fluid_volume_elements = u.getElementsFromEntities(gmsh,fluid_volume[0],fluid_volume_entities)

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

    #gmsh.finalize()

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

    IIf,JJf,Vffk,Vffm = libF_tet10.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density

    fluid_ndof = datafluidmesh['nodes'].shape[0]

    HFF=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################

    #print(libFreeSurf.globalacousticmatrices.__doc__)

    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



    ##############################################################
    # Compute Standard Fluid load : rigid body motion of tank
    ##############################################################

    #print(libF.sloshimposedacc.__doc__)

    CF,VecNormalElts = libF_tet10.sloshimposedacc(
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

        eigen_values,eigen_vectors= spla.eigsh(HFF,
                                                                dataPb['flag_nb_eigen_modes'],
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

        f=open(results_file.as_posix() +'_eigen_frequencies.pck','wb')
        pickle.dump(freq_eigv_S, f)
        f.close()


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

        #    sol = spla.spsolve(KFF-omega**2*MFF,FF)
            sol = mumps.spsolve( HFF-omega**2*SFF , CF*(-omega**2) , comm=mycomm )
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

def sloshing_flexible_baffle_tet10_xfem(dataPb,dataFluid,dataStructure,mesh_file_fluid,mesh_file_stiffener,results_file):
    
    ##############################################################
    # Load fluid mesh and tank
    ##############################################################

    tic = time.process_time()

    fluid_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_fluid.as_posix()+'.msh',3)
    fluid_volume_elements,Idfluid_volume_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',11,10)
    
    structure_surface_elements,Idnode_structure_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,20)
    free_fluid_surface_elements,Idnode_free_fluid_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_fluid.as_posix()+'.msh',9,30)

    print(' ----  structure_surface_elements : ',structure_surface_elements)
    print(' ----  free_fluid_surface_elements : ',free_fluid_surface_elements)
    print(' ----  fluid_volume_elements : ',fluid_volume_elements)

    datafluidmesh = dict()
    datafluidmesh['nodes'] = fluid_nodes
    datafluidmesh['fluid_volume_elts'] = fluid_volume_elements
    datafluidmesh['structure_surface_elts'] = structure_surface_elements
    datafluidmesh['free_fluid_surface_elts'] = free_fluid_surface_elements

    fluid_ndof = datafluidmesh['nodes'].shape[0]

    print('node 8 :', datafluidmesh['nodes'][8-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_volume',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'],11)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Tank_surfaces',datafluidmesh['nodes'],datafluidmesh['structure_surface_elts'],9)
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Fluid_Free_surface',datafluidmesh['nodes'],datafluidmesh['free_fluid_surface_elts'],9)

    ##############################################################
    # Load stiffener mesh / compute Level Set / Get enriched nodes and elements
    ##############################################################

    tic = time.process_time()

    stiffener_nodes=silex_lib_gmsh.ReadGmshNodes(mesh_file_stiffener.as_posix()+'.msh',3)
    stiffener_surface_elements,Idstiffener_surface_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',2,50)
    stiffener_edge_elements,Idnode_stiffener_edge_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',1,60)
    stiffener_base_edge_elements,Idnode_stiffener_base_edge_elements=silex_lib_gmsh.ReadGmshElements(mesh_file_stiffener.as_posix()+'.msh',1,70)

    dataXfemStiffener = dict()
    dataXfemStiffener['nodes'] = stiffener_nodes
    dataXfemStiffener['stiffener_surface_elements'] = stiffener_surface_elements
    dataXfemStiffener['stiffener_edge_elements'] = stiffener_edge_elements
    dataXfemStiffener['stiffener_imposed_acceleration_nodes'] = Idnode_stiffener_base_edge_elements

    stiffener_ndof = dataXfemStiffener['nodes'].shape[0]*6

    # compute Level Set to stiffener
    Stiffener_LS,Stiffener_distance = libF_levelset.computelevelset(datafluidmesh['nodes'],
                                                                dataXfemStiffener['nodes'],
                                                                dataXfemStiffener['stiffener_surface_elements'])
    # make perpendicular mesh to the stiffener along the edge
    Stiffener_tangent_nodes,Stiffener_tangent_mesh=libF_levelset.buildtangentedgemesh(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            dataXfemStiffener['stiffener_edge_elements'])

    # compute tangent Level Set to stiffener edge
    Stiffener_tangent_LS,tmp = libF_levelset.computelevelset(datafluidmesh['nodes'],
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
                                    11,
                                    [[[Stiffener_LS],'nodal',1,'Level set'],
                                    [[Stiffener_tangent_LS],'nodal',1,'Tangent Level set']])

    # Get enriched nodes and elements directly from the stiffener surface mesh
    EnrichedElementstmp1,NbEnrichedElements=libF_levelset.getsurfenrichedelements(dataXfemStiffener['nodes'],
                                                                            dataXfemStiffener['stiffener_surface_elements'],
                                                                            datafluidmesh['nodes'],
                                                                            datafluidmesh['fluid_volume_elts'][:,0:4])
    EnrichedElements=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))]) # here, start with 1 (fortran indexing)


    Enrichednodes = np.unique(datafluidmesh['fluid_volume_elts'][EnrichedElements-1])

    # gmsh output to check
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Enriched_Fluid_Elements',datafluidmesh['nodes'],datafluidmesh['fluid_volume_elts'][EnrichedElements-1],11)

    ##############################################################
    # Compute XFEM Fluid Matrices
    ##############################################################
    #print(libF_xfem.globalxfemacousticmatrices.__doc__)

    IIxf,JJxf,Vkaa,Vmaa,Vkfa,Vmfa = libF_tet10_xfem.globalxfemacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                                        datafluidmesh['nodes'],
                                                                        Stiffener_LS,
                                                                        Stiffener_tangent_LS,
                                                                        1.0,1.0)

    HAA=sps.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=sps.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )


    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,Vffm = libF_tet10.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################

    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

    ZEROSAA=HAA*0.0

    ##############################################################
    # Compute Standard Fluid load : rigid body motion of tank
    ##############################################################
    #print(libF.sloshimposedacc.__doc__)

    CF,VecNormalEltsF = libF_tet10.sloshimposedacc(
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
    #print(libF_xfem.sloshimposedacc_xfem1.__doc__)

    flag_write_quadrature_points_in_a_file=0

    CA,VecNormalEltsA = libF_tet10_xfem.sloshimposedacc_xfem1(
                            np.array(datafluidmesh['nodes']),
                            np.array(dataXfemStiffener['nodes']),
                            np.array(datafluidmesh['fluid_volume_elts']),
                            np.array(dataXfemStiffener['stiffener_surface_elements']),
                            EnrichedElements,
                            dataPb['U_dot_dot_imposed'],
                            flag_write_quadrature_points_in_a_file)

    ZEROCA=0.0*CA*dataFluid['rho']



    # check if normal vectors are pointing out of the fluid volume
    if dataPb['flag_write_gmsh_results']==1:
        silex_lib_gmsh.WriteResults(results_file.as_posix()+'_Mesh_Normal_to_stiffener',
                                dataXfemStiffener['nodes'],
                                dataXfemStiffener['stiffener_surface_elements'],
                                2,
                                [[VecNormalEltsA,'elemental',3,'Normal to stiffener elements']])


    ##############################################################
    # Compute Standard Structure Matrices
    ##############################################################

    IIs,JJs,Vks,Vms=libDKT.stiffnessmatrix(dataXfemStiffener['nodes'],
                                           dataXfemStiffener['stiffener_surface_elements'],
                                           [dataStructure['young'],dataStructure['nu'],dataStructure['thickness'],dataStructure['rho']])
    
    KSS=sps.csc_matrix( (Vks,(IIs,JJs)), shape=(stiffener_ndof,stiffener_ndof) ,dtype=float)
    MSS=sps.csc_matrix( (Vms,(IIs,JJs)), shape=(stiffener_ndof,stiffener_ndof) ,dtype=float)


    
    # Il faut ici enlever les lignes et colonnes de KSS MSS FS qui ne sont pas utilisees ...
    NodesSlist = np.array(range(1,len(np.unique(dataXfemStiffener['stiffener_surface_elements']))+1))
    
    dofS=np.hstack([(NodesSlist-1)*6,
                    (NodesSlist-1)*6+1,
                    (NodesSlist-1)*6+2,
                    (NodesSlist-1)*6+3,
                    (NodesSlist-1)*6+4,
                    (NodesSlist-1)*6+5]) 
    
    IdNodesFixed_x   =dataXfemStiffener['stiffener_imposed_acceleration_nodes']
    IdNodesFixed_y   =dataXfemStiffener['stiffener_imposed_acceleration_nodes']
    IdNodesFixed_z   =dataXfemStiffener['stiffener_imposed_acceleration_nodes']
    IdNodesFixed_rotx=dataXfemStiffener['stiffener_imposed_acceleration_nodes']
    IdNodesFixed_roty=dataXfemStiffener['stiffener_imposed_acceleration_nodes']
    IdNodesFixed_rotz=dataXfemStiffener['stiffener_imposed_acceleration_nodes']


    Fixed_Dof_S = np.hstack([(IdNodesFixed_x-1)*6,
                            (IdNodesFixed_y-1)*6+1,
                            (IdNodesFixed_z-1)*6+2,
                            (IdNodesFixed_rotx-1)*6+3,
                            (IdNodesFixed_roty-1)*6+4,
                            (IdNodesFixed_rotz-1)*6+5,                               
                               ])


    SolvedDofS = np.setdiff1d(dofS,Fixed_Dof_S)

    US_Imposed_acc=np.zeros(stiffener_ndof)
    US_Imposed_acc[(IdNodesFixed_x-1)*6]   = dataPb['U_dot_dot_imposed'][0]
    US_Imposed_acc[(IdNodesFixed_z-1)*6+1] = dataPb['U_dot_dot_imposed'][1]
    US_Imposed_acc[(IdNodesFixed_x-1)*6+2] = dataPb['U_dot_dot_imposed'][2]

    #for i in range(len(dataXfemStiffener['stiffener_imposed_acceleration_nodes'])):
    #    US_Imposed_acc[i*6]   = dataPb['U_dot_dot_imposed'][0]
    #    US_Imposed_acc[i*6+1] = dataPb['U_dot_dot_imposed'][1]
    #    US_Imposed_acc[i*6+2] = dataPb['U_dot_dot_imposed'][2]


    
    ##################################################################
    # Compute coupling FLEXIBLE BAFFLE STRUCTURE STIFFENER / FLUID terms
    ##################################################################

    #print(libF_tet10_xfem.computexfemcoupling1.__doc__)

    IIc1,JJc1,Vc1=libF_tet10_xfem.computexfemcoupling1(datafluidmesh['nodes'],
                                                       dataXfemStiffener['nodes'],
                                                       datafluidmesh['fluid_volume_elts'],
                                                       dataXfemStiffener['stiffener_surface_elements'],
                                                       EnrichedElements)
    #IIc2,JJc2,Vc2=libF_tet10_xfem.computexfemcoupling2(fluid_nodes1,struc_nodes,fluid_elements1,struc_elements,EnrichedElements,LevelSet)

    CSA=sps.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(stiffener_ndof,fluid_ndof) ) 
    CAS=CSA.T*dataFluid['rho']
    
    ##############################################################
    # Assemble the whole system
    ##############################################################

    SolvedDofF=list(range(fluid_ndof))
    SolvedDofA=Enrichednodes-1

#    H=sps.construct.bmat( [  [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA],None],
#                                        [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA],None],
#                                        [None,-CSA[SolvedDofS,:][:,SolvedDofA],KSS[SolvedDofS,:][:,SolvedDofS]]
#                                    ] )
#            
#    S=sps.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],  None,   None],
#                                     [None,HAA[SolvedDofA,:][:,SolvedDofA]*0.0,CAS[SolvedDofA,:][:,SolvedDofS]],
#                                     [None,None,MSS[SolvedDofS,:][:,SolvedDofS]]
#                                    ] )

    #C = np.array([*CF[SolvedDofF], *CA[SolvedDofA], *CS[SolvedDofA]])
    ##############################################################
    # Compute eigen modes
    ##############################################################

    if dataPb['flag_eigen_vectors']==1:

        # Fluid
        eigen_values_F,eigen_vectors_F= spla.eigsh(HFF,dataPb['flag_nb_eigen_modes'],SFF,sigma=0,which='LM')

        freq_eigv_F=list(np.sqrt(eigen_values_F)/(2*np.pi))
        print('XFEM eigen frequencies : ',freq_eigv_F)
        eigen_vector_list_F=[]
        for i in range(eigen_values_F.shape[0]):
            Q=np.zeros(fluid_ndof)
            Q=eigen_vectors_F[SolvedDofF,i]
            eigen_vector_list_F.append(Q)


        if dataPb['flag_write_gmsh_results']==1:
            silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_fluid_eigenmodes',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    4,
                                    [[eigen_vector_list_F,'nodal',1,'modes']])

        # stiffener, flexible baffle
        eigen_values_S,eigen_vectors_S= spla.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],dataPb['flag_nb_eigen_modes'],MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

        freq_eigv_S=list(np.sqrt(eigen_values_S)/(2*np.pi))

        eigen_vector_S_list=[]
        for i in range(eigen_values_S.shape[0]):
            Q=np.zeros(stiffener_ndof)
            Q[SolvedDofS]=eigen_vectors_S[:,i]
            disp=np.zeros((len(NodesSlist),3))
            disp[range(len(NodesSlist)),0]=Q[list(range(0,stiffener_ndof,6))]
            disp[range(len(NodesSlist)),1]=Q[list(range(1,stiffener_ndof,6))]
            disp[range(len(NodesSlist)),2]=Q[list(range(2,stiffener_ndof,6))]
            eigen_vector_S_list.append(disp)

        print ("structure eigen frequencies : ",freq_eigv_S)

        if dataPb['flag_write_gmsh_results']==1:
            silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_structure_eigenmodes',
                                         dataXfemStiffener['nodes'],
                                         dataXfemStiffener['stiffener_surface_elements'],
                                         2,
                                         [[eigen_vector_S_list,'nodal',3,'modes']])



    ##############################################################
    # Compute FRF
    ##############################################################
    if dataPb['flag_FRF']==1:
        PressSave=[]
        frequencies=[]
        QuantityOfInterest=[]
        #damping=None
        US_list=[]

        print ("Time at the beginning of the FRF:",time.ctime())
        for f in np.linspace(dataPb['freq_ini'],
                            dataPb['freq_end'],
                            dataPb['nb_freq_step']):
            
            print('Solve freq : ',f)
            frequencies.append(f)

            omega=2*np.pi*f
            forceType ='float'

            H=sps.construct.bmat( [  [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA],None],
                                                [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA],None],
                                                [None,CSA[SolvedDofS,:][:,SolvedDofA],KSS[SolvedDofS,:][:,SolvedDofS]]
                                            ] )

            S=sps.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],  None,   None],
                                             [None,ZEROSAA[SolvedDofA,:][:,SolvedDofA],-CAS[SolvedDofA,:][:,SolvedDofS]],
                                             [None,None,MSS[SolvedDofS,:][:,SolvedDofS]]
                                            ] )


            CS=(-KSS[SolvedDofS,:][:,Fixed_Dof_S]+omega**2*MSS[SolvedDofS,:][:,Fixed_Dof_S])*US_Imposed_acc[Fixed_Dof_S]
            C = np.array([*CF[SolvedDofF]*(-omega**2), *ZEROCA[SolvedDofA]*(-omega**2), *CS])
            #print(US_Imposed_acc)
            

            #sol = spla.spsolve(KFF-omega**2*MFF,FF)
            sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
            
            press      = np.zeros(fluid_ndof)
            press[SolvedDofF] = sol[SolvedDofF].copy()
            enrichment = np.zeros(fluid_ndof)
            enrichment[SolvedDofA]= sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))].copy()
            CorrectedPressure=np.array(press)
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA].T+np.array(enrichment[SolvedDofA]*np.sign(Stiffener_LS[SolvedDofA]).T)
            PressSave.append(CorrectedPressure)

            QuantityOfInterest.append(sol[8-1]) # upper corner

            Q=np.zeros(stiffener_ndof)
            Q[SolvedDofS]=sol[list(range(len(SolvedDofF)+len(SolvedDofA),len(SolvedDofF)+len(SolvedDofA)+len(SolvedDofS),1))]
            Q[Fixed_Dof_S]=US_Imposed_acc[Fixed_Dof_S]
            disp=np.zeros((len(NodesSlist),3))
            disp[range(len(NodesSlist)),0]=Q[list(range(0,stiffener_ndof,6))]
            disp[range(len(NodesSlist)),1]=Q[list(range(1,stiffener_ndof,6))]
            disp[range(len(NodesSlist)),2]=Q[list(range(2,stiffener_ndof,6))]
            US_list.append(disp)
            
            


        frfsave=[np.array(frequencies),np.array(QuantityOfInterest)]

        if dataPb['flag_write_gmsh_results']==1:
            print('QuantityOfInterest : ',QuantityOfInterest)
            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_fluid_frf',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    4,
                                    [[PressSave,'nodal',1,'pressure']]
                                    )
            silex_lib_gmsh.WriteResults2(results_file.as_posix() +'_results_struct_frf',
                                    dataXfemStiffener['nodes'],
                                    dataXfemStiffener['stiffener_surface_elements'],
                                    2,
                                    [[US_list,'nodal',3,'displacement']]
                                    )
        print ("Time at the end of the FRF:",time.ctime())




        f=open(results_file.as_posix() +'_results.frf','wb')
        pickle.dump(frfsave, f)
        f.close()
    return



class compute_sloshing():
    def __init__(self, dataPb, dataFluid, files):        
        self.dataPb = dataPb
        self.dataFluid = dataFluid
        self.mesh_file_struct = files.get('struct', None)
        self.mesh_file_fluid = files.get('fluid', None)
        self.results_file = files.get('results', None)
        self.results_file_basis = self.results_file
        #
        self._SolvedDofF = []
        self._SolvedDofA = []
        self.fluid_nodes = []
        self.fluid_elements = []
        self.fluid_id_elements = []
        self.struct_nodes = []
        self.struct_elements = []
        self.struct_id_elements = []
        # 
        self.op = dict()
        
        self._init()
        
    def _init(self):
        self.create_dir()
        # generate geometry and mesh files 
        self.load_fluid()
        # pre-processing
        self.pre_process()
        # compute offline operators
        self.compute_operators()
        # compute loads
        self.compute_loads()
        # 
        self.show_data()
    
    @property
    def SolvedDofF(self):
        if len(self._SolvedDofF) == 0:
            self._SolvedDofF = list(range(self.fluid_ndof))
        return self._SolvedDofF
    @property
    def SolvedDofA(self):
        if len(self._SolvedDofA) == 0:
            self._SolvedDofA = list(range(self.fluid_ndof))
        return self._SolvedDofA
    @property
    def nbSolvedDofF(self):
        return len(self.SolvedDofF)
    @property
    def nbSolvedDofA(self):
        return len(self.SolvedDofA)
    @property
    def fluid_ndof(self):
        return len(self.fluid_nodes)
    
    @property
    def shell_order(self):
        db = {'TRI3': 1, 'TRI6': 2}
        return db.get(self.dataPb.get('shell_element'))
    
    @property
    def fluid_order(self):
        db = {'TET4': 1, 'TET10': 2}
        return db.get(self.dataPb.get('fluid_element'))
    
    @property
    def method_sl(self):
        loadsolver = self.dataPb.get('method_sloshing', None)
        if loadsolver is None:
            loadsolver = 'mumps'
        return loadsolver
    
    @property
    def libFEM(self):
        if self.fluid_order == 1:
            return libF_tet4
        elif self.fluid_order == 2:
            return libF_tet10
        else:
            raise ValueError('Unsupported fluid element type')
    @property
    def libXFEM(self):
        if self.fluid_order == 1:
            return libF_tet4_xfem
        elif self.fluid_order == 2:
            return libF_tet10_xfem
        else:
            raise ValueError('Unsupported fluid element type')
    @property
    def libLS(self):
        if self.fluid_order == 1 or self.fluid_order == 2:
            return libF_tet4_xfem
        else:
            raise ValueError('Unsupported fluid element type')    
    
            
    @property
    def libFreeSurf(self):
        if self.fluid_order == 1:
            return libFreeSurf_tet4
        elif self.fluid_order == 2:
            return libFreeSurf_tet10
        else:
            raise ValueError('Unsupported fluid element type')

    @property
    def enrich(self):
        return self.dataPb.get('enrichment', False)
    
    def create_dir(self):
        if self.mesh_file_struct:
            self.mesh_file_struct.parent.mkdir(parents=True, exist_ok=True)
        if self.mesh_file_fluid:
            self.mesh_file_fluid.parent.mkdir(parents=True, exist_ok=True)
        if self.results_file_basis:
            self.results_file_basis.parent.mkdir(parents=True, exist_ok=True)
        
    @utils.timeit('Load fluid mesh')
    def load_fluid(self):        
        self.fluid_nodes = silex_lib_gmsh.ReadGmshNodes(self.mesh_file_fluid.as_posix()+'.msh',3)
        data = silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',11,10)
        self.fluid_elements, self.fluid_id_elements = data
        data = silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',9,30)
        self.free_fluid_elements, self.free_fluid_id_elements = data
        data = silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',9,20)
        self.fluid_bounds_elements, self.fluid_bounds_id_elements = data
        pass
    
    @utils.timeit('Load structure mesh')
    def load_struct(self):
        self.struct_nodes=silex_lib_gmsh.ReadGmshNodes(self.mesh_file_struct.as_posix()+'.msh',3)
        data=silex_lib_gmsh.ReadGmshElements(self.mesh_file_struct.as_posix()+'.msh',2,50)
        self.struct_elements, self.struct_id_elements = data
        data=silex_lib_gmsh.ReadGmshElements(self.mesh_file_struct.as_posix()+'.msh',1,60)
        self.struct_edge_elements, self.struct_edge_id_elements = data
        pass
        
    def show_data(self):
        logger.info('Nb fluid nodes: {}'.format(len(self.fluid_nodes)))
        logger.info('Nb fluid elements: {}'.format(len(self.fluid_elements)))
        logger.info('Nb free fluid elements: {}'.format(len(self.free_fluid_elements)))
        logger.info('Nb fluid boundary elements: {}'.format(len(self.fluid_bounds_elements)))
        if self.enrich and self.struct_nodes:
            logger.info('Nb structure nodes: {}'.format(len(self.structure_nodes)))
            logger.info('Nb structure elements: {}'.format(len(self.structure_elements)))
        pass
    
    def export(self,kind='fluid'):
        if isinstance(kind,str):
            kind=[kind]
        for k in kind:
            if k=='fluid':
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_volume',self.fluid_nodes,self.fluid_elements,11)
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Tank_surfaces',self.fluid_bounds_nodes,self.fluid_bounds_elements,9)
            elif k=='free_surface':
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_Free_surface',self.fluid_nodes,self.free_fluid_elements,9)
            elif k=='structure':
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Stiffener_surface',
                                            self.struct_nodes,
                                            self.structure_elements,2)
            elif k=='levelset':
                silex_lib_gmsh.WriteResults2(self.results_file.as_posix()+'_LevelSet',
                                            self.fluid_nodes,
                                            self.fluid_elements,
                                            11,
                                            [[[self.struct_LS],'nodal',1,'Level set'],
                                            [[self.struct_LS_tangent],'nodal',1,'Tangent Level set']])
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Enriched_Fluid_Elements',
                                            self.fluid_nodes,
                                            self.fluid_elements,11)
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Enriched_Fluid_Elements',
                                            self.fluid_nodes,
                                            self.fluid_elements[self.enriched_elements-1],11)
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Normal_to_stiffener',
                                self.struct_nodes,
                                self.struct_elements,
                                2,
                                [[self.vecNormalEltsA,'elemental',3,'Normal to stiffener elements']])
            elif k=='TET10toTET4':
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_volume_tet10TOtet4',
                                            self.fluid_nodes,
                                            self.dataTET10toTET4['elements'],4)
                silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Enriched_Fluid_Elements_tet10TOtet4',
                                            self.fluid_nodes,
                                            self.dataTET10toTET4['elements'][self.dataTET10toTET4['enriched_elements']-1],4)
            elif k=='eigenfrequencies':
                # write csv file
                with open(self.results_file.as_posix()+'_eigen_frequencies.csv', mode='w', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    writer.writerow(['#Mode number', 'Frequency (Hz)'])
                    for i,freq in enumerate(self.eigen_frequencies):
                        writer.writerow([i+1, freq])
                # write pickle
                with open(self.results_file.as_posix()+'_eigen_frequencies.pkl', 'wb') as f:
                    pickle.dump(self.eigen_frequencies, f)
            elif k=='eigenmodes':
                silex_lib_gmsh.WriteResults2(self.results_file.as_posix()+'_Eigen_modes',
                                             self.fluid_nodes,
                                             self.fluid_elements,
                                             11,
                                             [[self.eigen_vectors,'nodal',1,'modes'],
                                              [self.eigen_vectors_uncorrected,'nodal',1,'press classic'],
                                              [self.eigen_vectors_enrichment,'nodal',1,'press enrich']])
                silex_lib_gmsh.WriteResults2(self.results_file.as_posix() +'_results_fluid_eigenmodes_on_tet4mesh',
                                             self.fluid_nodes,
                                             self.dataTET10toTET4['elements'],
                                             4,
                                             [[self.eigen_vectors,'nodal',1,'pressure']]
                                             )

                objMesh = lib.MeshField(nodes=self.fluid_nodes, 
                                        elems=self.dataTET10toTET4['elements'], 
                                        leveset=self.struct_LS, 
                                        levelsetTg=self.struct_LS_tangent)
                objMesh.addField(uncorrectedField=self.eigen_vectors_uncorrected,
                                enrichmentField=self.eigen_vectors_enrichment)
                datameshfield = objMesh.getData()
            
                # prepare fields
                dataW = []
                dataW.append({'data':self.struct_LS,'type':'nodal', 'name':'levelset'})
                dataW.append({'data':self.struct_LS_tangent,'type':'nodal','name':'tangent levelset'})
                # 
                dataW.append({
                    'name': 'press',
                    'nbsteps': self.eigen_vectors.shape[1],
                    'type': 'nodal',
                    'data': datameshfield['fields']
                })
                dataW.append({
                    'name': 'uncorrected press',
                    'nbsteps': self.eigen_vectors.shape[1],
                    'type': 'nodal',
                    'data': datameshfield['uncorrected']
                })
                dataW.append({
                    'name': 'correction press',
                    'nbsteps': self.eigen_vectors.shape[1],
                    'type': 'nodal',
                    'data': datameshfield['correction']
                })
                # export mesh
                msh2.mshWriter(
                    filename= results_file.as_posix() +'_results_fluid_eigenmodes_meshfield3D.msh',
                    nodes=datameshfield['nodes'],
                    elements=[{'type': 'TET4', 'connectivity': datameshfield['TET4']},
                              {'type': 'PRI6', 'connectivity': datameshfield['PRI6']}],
                    fields=dataW,
                    append=True
                    )
        else:
            pass
        pass
    
    @utils.timeit('Compute Level Set')
    def compute_LS(self):
        # compute LS
        data = self.libLS.computelevelset(self.fluid_nodes,
                                             self.struct_nodes,
                                             self.struct_elements)
        self.struct_LS, self.struct_distance = data
        # compute tangent LS
        data = libF_levelset.buildtangentedgemesh(self.struct_nodes,
                                                  self.struct_elements,
                                                  self.struct_edge_elements)
        self.struct_LS_tangent_nodes, self.struct_LS_tangent_elements = data
        # compute tangent LS to edge
        data = libF_levelset.computelevelset(self.fluid_nodes,
                                             self.struct_LS_tangent_nodes,
                                             self.struct_LS_tangent_elements)
        self.struct_LS_tangent, _ = data
        # Get enriched nodes and elements directly from the stiffener surface mesh
        data = libF_levelset.getsurfenrichedelements(self.struct_nodes,
                                                     self.struct_elements,
                                                     self.fluid_nodes,
                                                     self.fluid_elements[:,0:4])
        enriched_elements_tmp, nb_enriched_elements = data
        self.enriched_elements = np.unique(enriched_elements_tmp[list(range(nb_enriched_elements))]) # here, start with 1 (fortran indexing)
        self.enriched_nodes = np.unique(self.fluid_elements[self.enriched_elements-1])
    
    @utils.timeit('Enforce computation of TET4 from TET10')
    def compute_LS_TET1OtoTET4(self):
        self.dataTET10toTET4['elements'] = libF_tet10_xfem.tet10totet4(self.fluid_nodes,
                                                                        self.fluid_elements)
            
        # Get enriched nodes and elements directly from the stiffener surface mesh
        data =libF_levelset.getsurfenrichedelements(self.struct_nodes,
                                                    self.struct_elements,
                                                    self.fluid_nodes,
                                                    self.dataTET10toTET4['elements'])
        EnrichedElementstet4tmp1,NbEnrichedElementstet4= data
        self.dataTET10toTET4['enriched_elements'] = np.unique(EnrichedElementstet4tmp1[list(range(NbEnrichedElementstet4))])  # here, start with 1 (fortran indexing)
        self.dataTET10toTET4['enriched_nodes'] = np.unique(self.dataTET10toTET4['elements'][self.dataTET10toTET4['enriched_elements']-1])
        # nodes enriched for tet10 / but no enriched for tet4
        self.dataTET10toTET4['special_nodes'] = np.setdiff1d(self.enriched_nodes,
                                                             self.dataTET10toTET4['enriched_nodes']) 
    
    @utils.timeit('Compute XFEM operators')
    def compute_operators_online(self):
        
        data = self.libXFEM.globalxfemacousticmatrices(self.fluid_elements,
                                                   self.fluid_nodes,
                                                   self.struct_LS,
                                                   self.struct_LS_tangent*0.0-1.0, # enforce all elements are considered
                                                   1.0,1.0)
        IIxf,JJxf,Vkaa,Vmaa,Vkfa,Vmfa = data
        # build matrices
        self.op['HAA'] = sps.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(self.fluid_ndof, self.fluid_ndof) )
        self.op['HFA'] = sps.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(self.fluid_ndof, self.fluid_ndof) )                                                 
                            
    @utils.timeit('Compute FEM operators')
    def compute_operators(self):
        ## compute volume matrix
        data = self.libFEM.globalacousticmatrices(self.fluid_elements,
                                                    self.fluid_nodes,
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density
        IIf,JJf,Vffk,Vffm = data
        # build matrix
        self.op['HFF']=sps.csc_matrix( (Vffk,(IIf,JJf)), shape=(self.fluid_ndof, self.fluid_ndof) )
        
        ## compute free surface matrix
        data = self.libFreeSurf.globalacousticmatrices(self.free_fluid_elements,
                                                       self.fluid_nodes[:,[0,1]],
                                                       1.0,1.0)
        IIf,JJf,_,VSFF= data
        # build matrix
        self.op['SFF']=1/self.dataPb['g']*sps.csc_matrix( (VSFF,(IIf,JJf)), shape=(self.fluid_ndof, self.fluid_ndof) )
    
    @utils.timeit('Compute loads operators')
    def compute_loads(self):
        # Standard Fluid load : rigid body motion of tank
        data = self.libFEM.sloshimposedacc(self.fluid_nodes,
                                      self.fluid_bounds_elements,
                                      self.dataPb['U_dot_dot_imposed'])
        CF,self.vecNormalEltsF = data
        self.op['CF'] = CF*self.dataFluid['rho']
        
    @utils.timeit('Compute XFEM loads operators')
    def compute_xfem_loads(self):
        # XFEM Fluid load : rigid body motion of tank
        data = self.libXFEM.sloshimposedacc_xfem1(
                            np.array(self.fluid_nodes),
                            np.array(self.struct_nodes),
                            np.array(self.fluid_elements),
                            np.array(self.struct_elements),
                            self.enriched_elements,
                            self.dataPb['U_dot_dot_imposed'],
                            self.dataPb['flag_write_quadrature_points'])
        CA,self.vecNormalEltsA = data
        self.op['CA']=CA*self.dataFluid['rho']
        
    @utils.timeit('Assemble system')
    def assemble(self):
        
        self.op['H'] = sps.construct.bmat([[self.op['HFF'][self.SolvedDofF,:][:,self.SolvedDofF],
                                            self.op['HFA'][self.SolvedDofF,:][:,self.SolvedDofA]],
                                           [self.op['HFA'][self.SolvedDofA,:][:,self.SolvedDofF],
                                            self.op['HAA'][self.SolvedDofA,:][:,self.SolvedDofA]]])
        self.op['S'] = sps.construct.bmat([[self.op['SFF'][self.SolvedDofF,:][:,self.SolvedDofF],None],
                                           [None,np.zeros((self.nbSolvedDofA,self.nbSolvedDofA))]])
        self.op['C'] = np.hstack([self.op['CF'][self.SolvedDofF],
                                  self.op['CA'][self.SolvedDofA]])
    
    @utils.timeit('Compute eigen modes')
    def compute_eigenmodes(self):
        # compute eigen values and modes
        self.eigen_values, self.raw_eigen_vectors = compute_eigen(self.op['H'],
                                                              self.dataPb['flag_nb_eigen_modes'],
                                                              self.op['S'],
                                                              sigma=0,
                                                              which='LM')
        # compute frequencies
        self.eigen_frequencies = np.sqrt(self.eigen_values)/(2*np.pi)
        logger.info('Eigen frequencies: {}/{}/{}/{} (min/avg/max/nb)'.format(self.eigen_frequencies.min(),
                                                                      self.eigen_frequencies.mean(),
                                                                      self.eigen_frequencies.max(),
                                                                      len(self.eigen_frequencies)))
        #
        self.eigen_vectors = np.zeros((self.fluid_ndof, len(self.eigen_values)))
        self.eigen_vectors_uncorrected = self.eigen_vectors.copy()
        self.eigen_vectors_enrichment = self.eigen_vectors.copy()
        # compute eigenvectors 
        self.eigen_vectors_uncorrected[self.SolvedDofF,:] = self.raw_eigen_vectors[self.SolvedDofF,:]
        self.eigen_vectors_enrichment[self.SolvedDofA,:] = self.raw_eigen_vectors[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1)),:]
        self.eigen_vectors = self.eigen_vectors_uncorrected.copy()
        self.eigen_vectors[self.SolvedDofA,:] = self.eigen_vectors[self.SolvedDofA,:] \
            + np.sign(self.struct_LS[self.SolvedDofA])*self.eigen_vectors_enrichment[self.SolvedDofA,:]
        
    def generate_formatted_id(self, paraval, paranames):
        # generate an id string based on parameter values
        id_formatted = ''
        for i,pname in enumerate(paranames): 
            parastr = f'{int(1e3*paraval[i]):+04d}'
            id_formatted += '{}_{}'.format(pname, parastr.replace('+','p').replace('-','m'))
            if i<len(paranames)-1:
                id_formatted += '_'
        return id_formatted
    
    def run_parametric(self, param_list=None):
        #
        if isinstance(param_list,dict):
            para_val = param_list.get('values', None)
            para_names = param_list.get('names', None)
            freq_list = param_list.get('freq_list', None)
        #
        if not isinstance(para_val, np.ndarray):
            para_val = np.vstack(param_list)
        if para_names is None:
            para_names = ['p{}'.format(i) for i in range(para_val.shape[1])]
        logger.info('Run parametric study for {} parameters and {} sets'.format(para_val.shape[1], para_val.shape[0]))
        # run along each parameter set
        results = []
        for i,pset in enumerate(para_val):
            logger.info('Run parametric set {}/{}: {}'.format(i+1, para_val.shape[0], pset))
            # update results file names
            id_formatted = self.generate_formatted_id(pset, para_names)
            self.results_file = self.results_file.parent / (self.results_file.stem + '_' + id_formatted)
            # run pre-process
            self.pre_process_online(pset)
            # build operators
            self.compute_operators_online()
            # build loads
            self.compute_xfem_loads()
            # assemble system
            self.assemble()
            # run frequencies
            data = self.run_frequencies(freq_list=freq_list)
            results.append(data)
        return results

    
    @utils.timeit('Run computation of FRF')
    def run_frequencies(self, freq_list=None):
        # get the list of frequencies
        if freq_list is None:
            freq_list = np.linspace(self.dataPb['freq_ini'],
                                    self.dataPb['freq_end'],
                                    self.dataPb['nb_freq_step'])
        press = []
        QoI = []
        for it,f in enumerate(freq_list):
            logger.info('Solve freq {}/{}: {:5g} Hz'.format(it+1,len(freq_list),f))
            data = self.run_one_freq(f)
            press.append(data[0])
            QoI.append(data[1])
        self.press = np.vstack(press).T
        self.QoI = np.reshape(np.array(QoI),shape=(len(QoI[0]),len(QoI)))
        return {'press': self.press,
                'QoI': self.QoI}
        
    def pre_process_offline(self, param_val):
        mesher_cube_tank.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,self.mesh_fluid_order,self.mesh_file_fluid)
        
    def pre_process_online(self, para_val):
        # build mesh file
        lx = self.dataPb.get('lx')
        ly = self.dataPb.get('ly')
        lz = self.dataPb.get('lz')
        struct_lx = self.dataPb.get('struct_lx')
        struct_lz = self.dataPb.get('struct_lz')
        struct_mesh_size = self.dataPb.get('struct_mesh_size')
        #
        lx_up = para_val[0]
        lx_down = para_val[1]
        if self.enrich:
            mesher_cube_tank.Stiffener_DKT(lx,ly,lz,
                                           struct_lx,
                                           lx_up,
                                           lx_down,
                                           struct_lz,
                                           struct_mesh_size,
                                           self.shell_order,
                                           self.mesh_file_struct)
            # load struct
            self.load_struct()
            # build level set
            self.compute_LS()
        
    def pre_process(self):
        # prepare data for post-processing
        self.dataPb['post-processing'] = {}
        # find id of QoI node(s)
        if self.dataPb.get('QoI_node_ids', None):
            self.dataPb['post-processing']['QoI_node_ids'] = self.dataPb.get('QoI_node_ids', [8])
        if self.dataPb.get('QoI_nodes_bbx', None):
            idNodes = misc_utils.getNodesBBX(self.fluid_nodes, self.dataPb['QoI_nodes_bbx'])
            self.dataPb['post-processing']['QoI_node_ids'] = idNodes
        pass
        
    def post_process(self, freq, sol, ):
        # extract QoI
        idNodes = self.dataPb['post-processing']['QoI_node_ids']
        QoI = sol[np.array(idNodes)]
        return QoI
        
    def run_one_freq(self, freq):
        #
        omega=2*np.pi*freq
        forceType ='float'
        # solve linear system                
        sol = solve_linear(self.method_sl, 
                            self.op['H']-omega**2*self.op['S'], 
                            -omega**2*self.op['C'] , comm=mycomm )
        if self.enrich:
            press      = np.zeros(self.fluid_ndof)
            press[self.SolvedDofF] = sol[self.SolvedDofF].copy()
            enrichment = np.zeros(self.fluid_ndof)
            enrichment[self.SolvedDofA]= sol[list(range(len(self.SolvedDofF),len(self.SolvedDofF)+len(self.SolvedDofA),1))].copy()
            CorrectedPressure=np.array(press)
            CorrectedPressure[self.SolvedDofA]=CorrectedPressure[self.SolvedDofA].T+np.array(enrichment[self.SolvedDofA]*np.sign(self.struct_LS[self.SolvedDofA]).T)
            press = CorrectedPressure
        else:
            press = np.zeros(self.fluid_ndof)
            press[self.SolvedDofF] = sol[self.SolvedDofF].copy()
        self.press = sol.copy()
        # run post-processing
        QoI = self.post_process(freq,sol)
        
        return press,QoI
            
        
       