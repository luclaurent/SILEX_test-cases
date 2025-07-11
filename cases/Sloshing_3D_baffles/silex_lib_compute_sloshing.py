import string
import time
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.sparse.construct
from loguru import logger as log

import pylab as pl
import pickle

import sys
from pathlib import Path
sys.path.append('/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/')
import pymumps as mumps
import gmsh
# import utils as u
# import utils_acoustics as ua
from meshRW import msh2

# for tet4
from SILEXlib import silex_lib_acou_tet4 as libF_tet
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_tet_xfem
from SILEXlib import silex_lib_acou_tri3 as libFreeSurf_tet

# for tet4 and tet10 : level set
from SILEXlib import silex_lib_xfem_acou_tet4 as libF_levelset

# for post-processing XFEM results
from SILEXlib import MeshField as lib

# for DKT
from SILEXlib import silex_lib_dkt as libDKT

#from SILEXlib import silex_lib_dkt as libS
from SILEXlib import silex_lib_gmsh

# for baffle geometry
import baffles_geometries as bgeo

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc = comm.Get_size()
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

results_dir = Path(__file__).parent / 'results'


class sloshing_baffle:
    def __init__(self, dataPb, dataFluid, mesh_file_fluid, baffle_function, results_file, convert_from_tet10=False):
        self.dataPb = dataPb
        self.dataFluid = dataFluid
        self.mesh_file_fluid = mesh_file_fluid
        self.baffle_function = baffle_function
        self.results_file = results_dir / results_file.name
        self.convert_from_tet10 = convert_from_tet10
        #
        self.variables = {
            'g': dataFluid.get('rho', 9.81),  # gravity
            'rho': dataFluid.get('rho', 1000.0),  # density of the fluid
            'c': dataFluid.get('c', 1482.0),  # celerity of the fluid
            'nb_eigen_modes': dataPb.get('nb_eigen_modes', 10),  # number of eigen modes to compute
                         }
        #
        self.flag_write_quadrature_points_in_a_file = False
        # 
        self.fluid_mesh = dict()
        self.baffle_mesh = dict()
        self.op = dict()
        self.eigen = dict()
        self.solution = dict()
        self.solutions = list()
        #
        self.tic = time.process_time()
        #
        self.build_fluid_data()
        self.build_standard_fluid_op()
        self.build_load_op()
        
        
    def build_fluid_data(self):        
        
        fluid_nodes=silex_lib_gmsh.ReadGmshNodes(self.mesh_file_fluid.as_posix()+'.msh',3)
        fluid_volume_elements,_=silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',4,10)
        structure_surface_elements,_=silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',2,20)
        free_fluid_surface_elements,_=silex_lib_gmsh.ReadGmshElements(self.mesh_file_fluid.as_posix()+'.msh',2,30)
        
        self.fluid_mesh['nodes'] = fluid_nodes
        self.fluid_mesh['fluid_volume_elts'] = fluid_volume_elements
        self.fluid_mesh['structure_surface_elts'] = structure_surface_elements
        self.fluid_mesh['free_fluid_surface_elts'] = free_fluid_surface_elements
        self.fluid_mesh['ndofs'] = fluid_nodes.shape[0]
        
        log.info(' ----  fluid mesh built from file : {}'.format(self.mesh_file_fluid))
        log.info(' ----  fluid nodes : {} nodes'.format(fluid_nodes.shape[0]))
        log.info(' ----  fluid dofs : {} dofs'.format(fluid_nodes.shape[0]))
        log.info(' ----  structure_surface_elements : {} elements'.format(structure_surface_elements.shape[0]))
        log.info(' ----  free_fluid_surface_elements : {} elements'.format(free_fluid_surface_elements.shape[0]))
        log.info(' ----  fluid_volume_elements : {} elements'.format(fluid_volume_elements.shape[0]))
        
        self.export('mesh_fluid')
        

    def build_baffle(self, parameters):
        # build geometry of the baffle
        self.baffle = self.baffle_function(parameters=parameters, active_parameters=self.dataPb.get('active_parameters', []))
        self.baffle.create_geometry(parameters)
        self.baffle_mesh['nodes'] = self.baffle.getNodes()
        self.baffle_mesh['elements'] = self.baffle.getElements()
        self.baffle_mesh['edge_elements'] = self.baffle.getElements(kind='edge')
        #
        log.info(' ----  baffle nodes : {} nodes'.format(self.baffle_mesh.get('nodes').shape[0]))
        log.info(' ----  baffle elements : {} elements'.format(self.baffle_mesh.get('elements').shape[0]))
        log.info(' ----  baffle edge elements : {} elements'.format(self.baffle_mesh.get('edge_elements').shape[0]))

        # compute level set to the baffle
        self.LS, self.distance = libF_levelset.computelevelset(self.fluid_mesh.get('nodes'),
                                                               self.baffle_mesh.get('nodes'),
                                                               self.baffle_mesh.get('elements'))
        # compute perpendicular mesh to the baffle edge
        baffle_tangent_nodes,baffle_tangent_mesh=libF_levelset.buildtangentedgemesh(self.baffle_mesh.get('nodes'),
                                                                            self.baffle_mesh.get('elements'),
                                                                            self.baffle_mesh.get('edge_elements'))
        # compute tangent level set to the baffle edge
        self.LStg,_ = libF_levelset.computelevelset(self.fluid_mesh.get('nodes'),
                                                    baffle_tangent_nodes,
                                                    baffle_tangent_mesh)
        # Get enriched nodes and elements directly from the stiffener surface mesh
        EnrichedElementstmp1,NbEnrichedElements=libF_levelset.getsurfenrichedelements(self.baffle_mesh.get('nodes'),
                                                                                self.baffle_mesh.get('elements'),
                                                                                self.fluid_mesh.get('nodes'),
                                                                                self.fluid_mesh.get('fluid_volume_elts')[:,0:4])
        self.enrichElements=np.unique(EnrichedElementstmp1[list(range(NbEnrichedElements))]) # here, start with 1 (fortran indexing)
        self.enrichNodes = np.unique(self.fluid_mesh.get('fluid_volume_elts')[self.enrichElements-1])
        #
        log.info(' ----  Enriched elements : {} elements'.format(self.enrichElements.shape[0]))
        log.info(' ----  Enriched nodes : {} nodes'.format(self.enrichNodes.shape[0]))
        self.export('baffle_mesh')
        self.export('enrich_data')
        
    def build_standard_fluid_op(self):
            ##############################################################
            # Compute Standard Fluid Matrices : VOLUME
            ##############################################################
            tic = time.process_time()
            # get operators

            IIf,JJf,Vffk,_ = libF_tet.globalacousticmatrices(self.fluid_mesh.get('fluid_volume_elts'),
                                                             self.fluid_mesh.get('nodes'),
                                                             1.0,
                                                             1.0) # we put 1 for celerity and 1 for density

            fluid_ndof = self.fluid_mesh.get('ndofs')
            self.op['HFF']=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )
            log.info(' ----  Standard Fluid Matrices : VOLUME computed in {:.2f} seconds'.format(time.process_time()-tic))
            # #############################################################
            # Compute Standard Fluid Matrices : FREE SURFACE
            # #############################################################
            # tic = time.process_time()
            # IIf,JJf,_,VSFF=libFreeSurf_tet.globalacousticmatrices(self.fluid_mesh.get('free_fluid_surface_elts'),
            #                                                       self.fluid_mesh.get('nodes')[:,[0,1]],
            #                                                       1.0,
            #                                                       1.0)
            # fluid_ndof = self.fluid_mesh.get('ndofs')
            # self.op['SFF']=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/self.variables.get('g')
            # log.info(' ----  Standard Fluid Matrices : FREE surface computed in {:.2f} seconds'.format(time.process_time()-tic))
            # ##############################################################
            # build list of standard dofs
            self.dofsF = list(range(self.fluid_mesh.get('ndofs')))  # convert to zero-based indexing
        
        
    def build_xfem_op(self):
        ##############################################################
        # Compute XFEM Fluid Matrices
        ##############################################################
        tic = time.process_time()
        # for tet10 only
        # IIxf,JJxf,Vkaa,_,Vkfa,_ = libF_tet_xfem.globalxfemacousticmatrices(self.fluid_mesh.get('fluid_volume_elts'),
        #                                                                    self.fluid_mesh.get('nodes'),
        #                                                                    self.LS,
        #                                                                    self.LStg*0.0-1.0,
        #                                                                    1.0,
        #                                                                    1.0)
        IIxf,JJxf,Vkaa,_,Vkfa,_ = libF_tet_xfem.computeedgeenrichment2(self.fluid_mesh.get('nodes'),
                                                                       self.fluid_mesh.get('fluid_volume_elts'),
                                                                       self.LS,
                                                                       self.LStg*0.0-1.0,
                                                                       1.0,
                                                                       1.0)
        fluid_ndof = self.fluid_mesh.get('ndofs')
        self.op['HAA']=scipy.sparse.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
        self.op['HFA']=scipy.sparse.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
        log.info(' ----  XFEM Fluid Matrices computed in {:.2f} seconds'.format(time.process_time()-tic))
        ##############################################################    
        # Compute XFEM Fluid Matrices : FREE SURFACE
        ##############################################################
        tic = time.process_time()
        IIf,JJf,_,VSFF=libFreeSurf_tet.globalacousticmatrices(self.fluid_mesh.get('free_fluid_surface_elts'),   
                                                                      self.fluid_mesh.get('nodes')[:,[0,1]],
                                                                      1.0,1.0)
        fluid_ndof = self.fluid_mesh.get('ndofs')
        self.op['SFF']=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/self.variables.get('g')
        log.info(' ----  XFEM Fluid Matrices : FREE surface computed in {:.2f} seconds'.format(time.process_time()-tic))
        
        # build list of enriched dofs
        self.dofsA = self.enrichNodes - 1  # convert to zero-based indexing

    def build_load_op(self):
        ##############################################################
        # Compute Standard Fluid load : rigid body motion of tank
        ##############################################################
        tic = time.process_time()
        CF,self.VecNormalEltsF = libF_tet.sloshimposedacc(
                            self.fluid_mesh.get('nodes'),
                            self.fluid_mesh.get('structure_surface_elts'),
                            self.dataPb.get('U_dot_dot_imposed'))

        self.op['CF'] = CF * self.variables.get('rho')
        log.info(' ----  Standard Fluid load computed in {:.2f} seconds'.format(time.process_time()-tic))
        self.export('load')
        
    def build_load_xfem_op(self):
        ##############################################################
        # Compute XFEM Fluid load : rigid body motion of tank
        ##############################################################
        tic = time.process_time()

        CA,self.VecNormalEltsA = libF_tet_xfem.sloshimposedacc_xfem1(
                                self.fluid_mesh.get('nodes'),
                                self.baffle_mesh.get('nodes'),
                                self.fluid_mesh.get('fluid_volume_elts'),
                                self.baffle_mesh.get('elements'),
                                self.enrichElements,
                                self.dataPb.get('U_dot_dot_imposed'),
                                self.flag_write_quadrature_points_in_a_file)

        self.op['CA']=CA * self.variables.get('rho')
        log.info(' ----  XFEM Fluid load computed in {:.2f} seconds'.format(time.process_time()-tic))
        self.export('load_xfem')
        
    def build_structure_op(self):
        # build structure operator
        # here, we assume that the structure is a rigid body
        # so we do not need to build a structure operator
        # but we can build a structure operator if needed
        self.op['HSS'] = None
        pass
    
    

    def assemble(self):
        # Assemble the global system
        tic = time.process_time()
        self.op['H'] = scipy.sparse.bmat([[self.op.get('HFF', None)[self.dofsF,:][:,self.dofsF],
                                           self.op.get('HFA', None)[self.dofsF,:][:,self.dofsA]],
                                          [self.op.get('HFA', None)[self.dofsA,:][:,self.dofsF], 
                                           self.op.get('HAA', None)[self.dofsA,:][:,self.dofsA]]],
                                         format='csc')
        #
        self.op['S'] = scipy.sparse.bmat([[self.op.get('SFF', None)[self.dofsF,:][:,self.dofsF], None],
                                          [None, self.op.get('SFF', None)[self.dofsA,:][:,self.dofsA]]],
                                         format='csc')
        log.info(' ----  Global system assembled in {:.2f} seconds'.format(time.process_time()-tic))
        
    def build_parametric_op(self, parameter):
        # build operators for parameter
        log.info(' ----  Building parametric operators for parameter : {}'.format(parameter))
        # build baffle geometry
        self.build_baffle(parameter)
        # build XFEM operators
        self.build_xfem_op()
        # build load operator
        self.build_load_xfem_op()
        # assemble the global system
        self.assemble()
    
    def compute_eigen_modes(self):
        tic = time.process_time()
        self.eigen['values'], self.eigen['vectors'] = eigen_values,eigen_vectors= scipy.sparse.linalg.eigsh(H,
                                                                self.variables.get('nb_eigen_modes', 10),
                                                                self.op['S'],
                                                                sigma=0,
                                                                which='LM')
        self.eigen['freq'] = np.sqrt(self.eigen['values'])/(2*np.pi)  # convert to frequency
        log.info(' ----  Eigen modes computed in {:.2f} seconds'.format(time.process_time()-tic))
        log.info(' ----  Eigen values : {}'.format(self.eigen['values']))
        self.export('eigen_freq')

    def linear_solve(self, A, b):
        if self.dataPb.get('mumps'):
            # use MUMPS solver
            log.info(' ----  Using MUMPS solver')
            tic = time.process_time()
            sol = mumps.spsolve( A, b , comm=mycomm )
            log.info(' ----  MUMPS solver solved in {:.2f} seconds'.format(time.process_time()-tic))
        else:
            # use scipy sparse solver
            log.info(' ----  Using scipy sparse solver')
            tic = time.process_time()
            sol = scipy.sparse.linalg.spsolve(A, b)
            log.info(' ----  scipy sparse solver solved in {:.2f} seconds'.format(time.process_time()-tic))
        return sol

    def solve(self, parameters):
        # build operators for parameter
        self.build_parametric_op(parameters)
        self.solution['parameters'] = parameters
        # solve the system
        tic = time.process_time()
        #
        freq_list = np.linspace(self.dataPb.get('freq_ini'),
                              self.dataPb.get('freq_end'),
                              self.dataPb.get('nb_freq_step'))
        for ixf,f in enumerate(freq_list):
            log.info(' ----  Solving for frequency ({}/{}) : {:.2f} Hz'.format(ixf,len(freq_list),f))
            #
            omega = 2 * np.pi * f  # convert frequency to angular frequency
            self.solution['freq'] = f
            self.solution['omega'] = omega
            # compute the right-hand side
            rhs = np.array([*self.op['CF'][self.dofsF]*(-omega**2), 
                            *self.op['CA'][self.dofsA]*(-omega**2)])
            # solve the system
            self.solution['raw'] = self.linear_solve(self.op.get('H')-omega**2*self.op['S'],rhs)
            # build fields based on XFEM raw solution
            self.solution['corrected'], self.solution['uncorrected'], self.solution['enrichment'] = self.build_fields(self.solution['raw'])
            # increment solutions
            self.solutions.append(self.solution.copy())
            self.export('solution')
        self.export('solutions')

                
    def solve_no_xfem(self):
        # solve the system without XFEM
        tic = time.process_time()
        #
        freq_list = np.linspace(self.dataPb.get('freq_ini'),
                              self.dataPb.get('freq_end'),
                              self.dataPb.get('nb_freq_step'))
        for ixf,f in enumerate(freq_list):
            log.info(' ----  Solving for frequency ({}/{}) : {:.2f} Hz'.format(ixf,len(freq_list),f))
            #
            omega = 2 * np.pi * f  # convert frequency to angular frequency
            self.solution['freq'] = f
            self.solution['omega'] = omega
            # compute the right-hand side
            rhs = np.array([*self.op['CF'][self.dofsF]*(-omega**2)])
            # solve the system
            self.solution['raw'] = self.linear_solve(self.op.get('HFF')-omega**2*self.op['SFF'],rhs)

        


    def build_fields(self, raw_solution):
        corrected_solution = np.zeros_like(self.dofsF)
        uncorrected_solution = np.zeros_like(self.dofsF)
        enrichment_solution = np.zeros_like(self.dofsF)
        #
        corrected_solution[self.dofsF] = raw_solution[self.dofsF].copy()
        uncorrected_solution[self.dofsF] = raw_solution[self.dofsF].copy()
        enrichment_solution[self.dofsA] = raw_solution[len(self.dofsF):len(self.dofsF)+len(self.dofsA)].copy()
        corrected_solution = corrected_solution + np.sign(self.LS) * enrichment_solution
        
        return corrected_solution, uncorrected_solution, enrichment_solution
        
    def post_process(self):
        # Post-process results
        pass
    
    def export(self, kind=''):
        if not self.results_file.parent.exists():
            self.results_file.parent.mkdir(parents=True, exist_ok=True)
        if kind == 'mesh_fluid':
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_volume',
                                        self.fluid_mesh.get('nodes'),
                                        self.fluid_mesh.get('fluid_volume_elts'),4)
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Tank_surfaces',
                                        self.fluid_mesh.get('nodes'),
                                        self.fluid_mesh.get('structure_surface_elts'),2)
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_Free_surface',
                                        self.fluid_mesh.get('nodes'),
                                        self.fluid_mesh.get('free_fluid_surface_elts'),2)
        elif kind == 'baffle_mesh':
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_baffle', 
                                        self.baffle_mesh.get('nodes'),
                                        self.baffle_mesh.get('elements'),2)
            silex_lib_gmsh.WriteResults2(self.results_file.as_posix()+'_LevelSet',
                                        self.fluid_mesh.get('nodes'),
                                        self.fluid_mesh.get('fluid_volume_elts'),
                                        4,
                                        [[[self.LS],'nodal',1,'Level set'],
                                        [[self.LStg],'nodal',1,'Tangent Level set']])
        elif kind == 'enrich_data':
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Enriched_Fluid_Elements',
                                        self.fluid_mesh.get('nodes'),
                                        self.fluid_mesh.get('fluid_volume_elts')[self.enrichElements-1],4)
        elif kind == 'load':
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces',
                                self.fluid_mesh.get('nodes'),
                                self.fluid_mesh.get('structure_surface_elts'),
                                2,
                                [[self.VecNormalEltsF,'elemental',3,'Normal to tank elements']])
        elif kind == 'load_xfem':
            silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Normal_to_baffle',
                        self.baffle_mesh.get('nodes'),
                        self.baffle_mesh.get('elements'),
                        2,
                        [[self.VecNormalEltsA,'elemental',3,'Normal to baffle elements']])
        elif kind == 'eigen_freq':
            # silex_lib_gmsh.WriteResults2(self.results_file.as_posix()+'_Eigen_modes',
            #                         self.fluid_mesh['nodes'],
            #                         self.fluid_mesh['fluid_volume_elts'],
            #                         11,
            #                         [[self.eigen['vectors_corrected'],'nodal',1,'modes']])
            with open(self.results_file.as_posix() +'_eigen_frequencies.pck','wb') as f:
                pickle.dump(self.eigen['freq'], f)
                f.close()
        elif kind == 'solution':
            # store the solution
            results_file = Path.joinpath(self.results_file.parent,'solution').with_suffix('.pck')
            with open(results_file, 'wb') as f:
                pickle.dump(self.solution, f)
                f.close()
                
            
        elif kind == 'solutions':
                        # store the solution
            results_file = Path.joinpath(self.results_file.parent,'solutions').with_suffix('.pck')
            with open(results_file, 'wb') as f:
                pickle.dump(self.solutions, f)
                f.close()
            objMesh = lib.MeshField(nodes=self.fluid_mesh.get('nodes'), 
                                    elems=self.fluid_mesh.get('fluid_volume_elts'), 
                                    leveset=self.LS, 
                                    levelsetTg=self.LStg)
            data_uncorrected = np.vstack([sol.get('uncorrected') for sol in self.solutions]).transpose()
            data_enrichment = np.vstack([sol.get('enrichment') for sol in self.solutions]).transpose()
            objMesh.addField(uncorrectedField=data_uncorrected,
                             enrichmentField=data_enrichment)
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
                'nbsteps': self.dataPb.get('nb_freq_step'),
                'type': 'nodal',
                'data': datameshfield['fields']
            })
            dataW.append({
                'name': 'uncorrected press',
                 'nbsteps': self.dataPb.get('nb_freq_step'),
                'type': 'nodal',
                'data': datameshfield['uncorrected']
            })
            dataW.append({
                'name': 'correction press',
                 'nbsteps': self.dataPb.get('nb_freq_step'),
                'type': 'nodal',
                'data': datameshfield['correction']
            })
            
            
            # msh2.mshWriter(
            #     filename= results_file.as_posix() +'_results_fluid_frf_init.msh',
            #     nodes=datameshfield['nodes'],
            #     elements=[{'type':'TET4','connectivity':self.fluid_mesh.get('fluid_volume_elts')[2:5,:]}],
            #     fields=dataW,
            #     append=True
            #     )
#            # export mesh
            msh2.mshWriter(
                filename= Path.joinpath(self.results_file.parent,'_results_fluid_frf_meshfield3D.msh'),
                nodes= datameshfield['nodes'],
                elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
                          {'type':'PRI6','connectivity':datameshfield['PRI6']}],
                fields=dataW,
                append=True
                )




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

    HAA=scipy.sparse.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=scipy.sparse.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )


    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,_ = libF_tet10.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

    ##############################################################
    # Compute XFEM Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

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

    H=scipy.sparse.construct.bmat( [ [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA]],
                                    [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )
    maxS = np.abs(SFF.max())    
    S=scipy.sparse.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],None],
                                     [None,np.zeros((len(SolvedDofA),len(SolvedDofA)))]] )
    # matdense= SFF[SolvedDofF,:][:,SolvedDofF].todense()
    # matdenseAA = HAA[SolvedDofA,:][:,SolvedDofA].todense()*0.0
    # matdenseFA = HFA[SolvedDofF,:][:,SolvedDofA].todense()*0.0
    # S=scipy.sparse.construct.bmat( [ [matdense+maxS/100,matdenseFA+maxS/100],
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

        eigen_values,eigen_vectors= scipy.sparse.linalg.eigsh(H,
                                                                10,
                                                                S,
                                                                sigma=0,which='LM')

        freq_eigv_S=list(np.sqrt(eigen_values)/(2*np.pi))
        print('XFEM eigen frequencies : ',freq_eigv_S)
        eigen_vector_list=[]
        for i in range(eigen_values.shape[0]):
            #Q=np.zeros(fluid_ndof)
            #Q=eigen_vectors[SolvedDofF,i]
            press = np.zeros(fluid_ndof)
            Q = eigen_vectors[:,i].copy()
            press[SolvedDofF] = Q[SolvedDofF]
            enrichment = np.zeros(fluid_ndof)
            enrichment[SolvedDofA]= Q[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))]
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
                                    [[eigen_vector_list,'nodal',1,'modes']])

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

        #    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
            C = np.array([*CF[SolvedDofF]*(-omega**2), *CA[SolvedDofA]*(-omega**2)])
            # sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
            sol = scipy.sparse.linalg.spsolve(H-omega**2*S, C)
            
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
            
            sol = scipy.sparse.linalg.spsolve(HFF-omega**2*SFF, CF[SolvedDofF]*(-omega**2))
            press_no_baffle.append(sol)
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

    print(' ----  stiffener_surface_elements : ',stiffener_surface_elements)
    print(' ----  stiffener_edge_elements : ',stiffener_edge_elements)


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

    Stiffener_tangent_LS = Stiffener_tangent_LS*0.0-1.0 # we put -1 to have the same sign as the Level Set

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

    HAA=scipy.sparse.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=scipy.sparse.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,Vffm = libF_tet4.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet4.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



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

    H=scipy.sparse.construct.bmat( [ [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA]],
                                    [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )
            
    S=scipy.sparse.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],None],
                                    [None,HAA[SolvedDofA,:][:,SolvedDofA]*0.0]
                                    ] )
    #C = np.array([*CF[SolvedDofF], *CA[SolvedDofA]])
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
            silex_lib_gmsh.WriteResults2(results_file.as_posix()+'_Eigen_modes_tet4',
                                    datafluidmesh['nodes'],
                                    datafluidmesh['fluid_volume_elts'],
                                    4,
                                    [[eigen_vector_list,'nodal',1,'modes']])

    ##############################################################
    # Compute FRF
    ##############################################################
    if dataPb['flag_FRF']==1:
        Correctedpress=[]
        press=[]
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

        #    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
            C = np.array([*CF[SolvedDofF]*(-omega**2), *CA[SolvedDofA]*(-omega**2)])
            sol = mumps.spsolve( H-omega**2*S , C , comm=mycomm )
            
            CorrectedPressure = np.zeros(fluid_ndof)
            press.append(sol[SolvedDofF].copy())
            CorrectedPressure[SolvedDofF] = sol[SolvedDofF].copy()
            soltmp = np.zeros_like(sol[SolvedDofF])
            soltmp[SolvedDofA] = sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA),1))].copy()

            enrichpress.append(soltmp)
            CorrectedPressure=CorrectedPressure+np.array(soltmp*np.sign(Stiffener_LS).T)
            Correctedpress.append(CorrectedPressure)
            
            QuantityOfInterest.append(sol[8-1]) # upper corner
            


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
    return

def sloshing_rigid_baffle_tet4(dataPb,dataFluid,mesh_file,results_file):
    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

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

    HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################



    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet4.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



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

        #    sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
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

    HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################

    #print(libFreeSurf.globalacousticmatrices.__doc__)

    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81



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

    HAA=scipy.sparse.csc_matrix( (Vkaa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )
    HFA=scipy.sparse.csc_matrix( (Vkfa,(IIxf,JJxf)), shape=(fluid_ndof, fluid_ndof) )


    ##############################################################
    # Compute Standard Fluid Matrices : VOLUME
    ##############################################################

    tic = time.process_time()
    # get operators

    IIf,JJf,Vffk,Vffm = libF_tet10.globalacousticmatrices(datafluidmesh['fluid_volume_elts'],
                                                    datafluidmesh['nodes'],
                                                    1.0,
                                                    1.0) # we put 1 for celerity and 1 for density


    HFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )

    ##############################################################
    # Compute Standard Fluid Matrices : FREE SURFACE
    ##############################################################

    IIf,JJf,Vffktmp,VSFF=libFreeSurf_tet10.globalacousticmatrices(datafluidmesh['free_fluid_surface_elts'],
                                                                datafluidmesh['nodes'][:,[0,1]],
                                                                1.0,1.0)

    SFF=scipy.sparse.csc_matrix( (VSFF,(IIf,JJf)), shape=(fluid_ndof, fluid_ndof) )/9.81

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
    
    KSS=scipy.sparse.csc_matrix( (Vks,(IIs,JJs)), shape=(stiffener_ndof,stiffener_ndof) ,dtype=float)
    MSS=scipy.sparse.csc_matrix( (Vms,(IIs,JJs)), shape=(stiffener_ndof,stiffener_ndof) ,dtype=float)


    
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

    CSA=scipy.sparse.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(stiffener_ndof,fluid_ndof) ) 
    CAS=CSA.T*dataFluid['rho']
    
    ##############################################################
    # Assemble the whole system
    ##############################################################

    SolvedDofF=list(range(fluid_ndof))
    SolvedDofA=Enrichednodes-1

#    H=scipy.sparse.construct.bmat( [  [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA],None],
#                                        [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA],None],
#                                        [None,-CSA[SolvedDofS,:][:,SolvedDofA],KSS[SolvedDofS,:][:,SolvedDofS]]
#                                    ] )
#            
#    S=scipy.sparse.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],  None,   None],
#                                     [None,HAA[SolvedDofA,:][:,SolvedDofA]*0.0,CAS[SolvedDofA,:][:,SolvedDofS]],
#                                     [None,None,MSS[SolvedDofS,:][:,SolvedDofS]]
#                                    ] )

    #C = np.array([*CF[SolvedDofF], *CA[SolvedDofA], *CS[SolvedDofA]])
    ##############################################################
    # Compute eigen modes
    ##############################################################

    if dataPb['flag_eigen_vectors']==1:

        # Fluid
        eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(HFF,20,SFF,sigma=0,which='LM')

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
        eigen_values_S,eigen_vectors_S= scipy.sparse.linalg.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],20,MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

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

            H=scipy.sparse.construct.bmat( [  [HFF[SolvedDofF,:][:,SolvedDofF],HFA[SolvedDofF,:][:,SolvedDofA],None],
                                                [HFA[SolvedDofA,:][:,SolvedDofF],HAA[SolvedDofA,:][:,SolvedDofA],None],
                                                [None,CSA[SolvedDofS,:][:,SolvedDofA],KSS[SolvedDofS,:][:,SolvedDofS]]
                                            ] )

            S=scipy.sparse.construct.bmat( [ [SFF[SolvedDofF,:][:,SolvedDofF],  None,   None],
                                             [None,ZEROSAA[SolvedDofA,:][:,SolvedDofA],-CAS[SolvedDofA,:][:,SolvedDofS]],
                                             [None,None,MSS[SolvedDofS,:][:,SolvedDofS]]
                                            ] )


            CS=(-KSS[SolvedDofS,:][:,Fixed_Dof_S]+omega**2*MSS[SolvedDofS,:][:,Fixed_Dof_S])*US_Imposed_acc[Fixed_Dof_S]
            C = np.array([*CF[SolvedDofF]*(-omega**2), *ZEROCA[SolvedDofA]*(-omega**2), *CS])
            #print(US_Imposed_acc)
            

            #sol = scipy.sparse.linalg.spsolve(KFF-omega**2*MFF,FF)
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

