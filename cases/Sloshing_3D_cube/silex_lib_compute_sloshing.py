import numpy as np
from typing import AnyStr

import scipy.sparse as sps
from loguru import logger

import joblib

from pathlib import Path

# useful tools
from SILEXrun import utils as misc_utils
from SILEXrun import baseFSI
from SILEXrun import solvers

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

        


class compute_sloshing(baseFSI.baseFSI):
    def __init__(self, dataPb, dataFluid, dataIO):
        super().__init__(dataPb, dataFluid, dataIO)        
        self._init()
        
    def _init(self):
        self.create_dir()
        # generate geometry and mesh files 
        if self.enrich:
            self.pre_process_offline()
            #
            self.load_fluid()
            # pre-processing
            self.pre_process()
            # compute offline operators
            self.compute_operators()
            # compute loads
            self.compute_loads()
            #
            self.export(kind=['fluid','free_surface']) 
            #
            self.show_data()
        
    def tag_db(self):
        return {
            'fluid_element': ['TET4','TET10'],
            'shell_element': ['TRI3','TRI6'],
            'enrichment': [True, False],
            'method_linear': ['mumps','pardiso','umfpack'],
            'tag_fluid_volume': 10,
            'tag_fluid_free_surface': 30,
            'tag_fluid_structure_surface': 20,
            'tag_struct_surface': 50,
            'tag_struct_edge': 60
        }
    
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
 
    
   
    
    @utils.timeit('Compute XFEM online operators')
    def compute_xfem_op_online(self):
        ## TODO: must be fixed
        if self.dataPb.get('fluid_element') == 'TET4':
            fun = lambda e,n,ls,lst,clty,rho: self.libXFEM.computeedgeenrichment2(n,e,ls,lst,clty,rho)
        elif self.dataPb.get('fluid_element') == 'TET10':
            fun = lambda e,n,ls,lst,clty,rho: self.libXFEM.globalxfemacousticmatrices(e,n,ls,lst,clty,rho)
        
        data = fun(self.fluid_elements,
                   self.fluid_nodes,
                   self.struct_LS,
                   self.struct_LS_tangent,#*0.0-1.0, # enforce all elements are considered
                   1.0,1.0)
        
        IIxf,JJxf,Vkaa,_,Vkfa,_ = data
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
        IIf,JJf,Vffk,_ = data
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
        # loads
        loadValue = None
        if self.dataPb.get('U_dot_dot_imposed', None) is not None:
            loadValue = self.dataPb.get('U_dot_dot_imposed')
            self.load_type = 'acceleration'
        if self.dataPb.get('U_dot_imposed', None) is not None:
            loadValue = self.dataPb.get('U_dot_imposed')
            self.load_type = 'velocity'
        if self.dataPb.get('U_imposed', None) is not None:
            loadValue = self.dataPb.get('U_imposed')
            self.load_type = 'displacement'
            
        # Standard Fluid load : rigid body motion of tank
        data = self.libFEM.sloshimposedacc(self.fluid_nodes,
                                      self.fluid_bounds_elements,
                                      loadValue)
        CF,self.vecNormalEltsF = data
        self.op['CF'] = CF*self.dataFluid['rho']

        
    @utils.timeit('Compute XFEM loads operators')
    def compute_xfem_loads(self):
        # loads
        if self.dataPb.get('U_dot_dot_imposed', None) is not None:
            loadValue = self.dataPb.get('U_dot_dot_imposed')
            self.load_type = 'acceleration'
        if self.dataPb.get('U_dot_imposed', None) is not None:
            loadValue = self.dataPb.get('U_dot_imposed')
            self.load_type = 'velocity'
        if self.dataPb.get('U_imposed', None) is not None:
            loadValue = self.dataPb.get('U_imposed')
            self.load_type = 'displacement'
            
        # XFEM Fluid load : rigid body motion of tank
        data = self.libXFEM.sloshimposedacc_xfem1(
                            np.array(self.fluid_nodes),
                            np.array(self.struct_nodes),
                            np.array(self.fluid_elements),
                            np.array(self.struct_elements),
                            self.enriched_elements,
                            loadValue,
                            self.dataPb['flag_write_quadrature_points'])
        CA,self.vecNormalEltsA = data
        self.op['CA']=CA*self.dataFluid['rho']
        
    @utils.timeit('Assemble system')
    def assemble(self):
        if self.enrich:
            self.op['H'] = sps.construct.bmat([[self.op['HFF'][self.SolvedDofF,:][:,self.SolvedDofF],
                                                self.op['HFA'][self.SolvedDofF,:][:,self.SolvedDofA]],
                                            [self.op['HFA'][self.SolvedDofA,:][:,self.SolvedDofF],
                                                self.op['HAA'][self.SolvedDofA,:][:,self.SolvedDofA]]])
            self.op['S'] = sps.construct.bmat([[self.op['SFF'][self.SolvedDofF,:][:,self.SolvedDofF],None],
                                            [None,np.zeros((self.nbSolvedDofA,self.nbSolvedDofA))]])
            self.op['C'] = np.hstack([self.op['CF'][self.SolvedDofF],
                                    self.op['CA'][self.SolvedDofA]])
        else:
            self.op['H'] = self.op['HFF']
            self.op['S'] = self.op['SFF']
            self.op['C'] = self.op['CF']
    
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
        self.para_names = para_names
        logger.info('Run parametric study for {} parameters and {} sets'.format(para_val.shape[1], para_val.shape[0]))
        # run along each parameter set
        results = []
        for i,pset in enumerate(para_val):
            logger.info('Run parametric set {}/{}: {}'.format(i+1, para_val.shape[0], pset))
            # save current parameter set and value
            self.current_param_set = pset            
            # update results file names
            id_formatted = misc_utils.generate_formatted_id(pset, para_names)
            self.results_file = self.results_file_basis.parent / (self.results_file_basis.stem + '_' + id_formatted)
            # run pre-process
            self.pre_process_online(pset)
            # load fluid
            if not self.enrich:
                self.load_fluid()
                self.pre_process()
                self.export(kind=['fluid','free_surface'])
            # build operators
            self.compute_operators_online()
            # build loads
            self.compute_loads_online()
            # export field
            if self.enrich:
                self.export(kind='load_LS')
                self.export(kind='levelset')
            self.export(kind='load')            
            # show data
            self.show_data()
            # assemble system
            self.assemble()
            # depending on cases
            if self.dataPb.get('flag_eigen_vectors', False):
                # compute eigen modes
                self.compute_eigenmodes()
                # export results
                self.export(kind=['eigenfrequencies','eigenmodes'])
            if self.dataPb.get('flag_FRF', False):
                # run frequencies
                data = self.run_frequencies(freq_list=freq_list)                
                results.append(data)
                # export results
                self.export(kind=['frf'])
        # export all parametric results
        self.export(kind='qoi')
        return results

    
    @utils.timeit('Run computation of FRF')
    def run_frequencies(self, freq_list=None):
        # get the list of frequencies
        if freq_list is None:
            freq_list = np.linspace(self.dataPb['freq_ini'],
                                    self.dataPb['freq_end'],
                                    self.dataPb['nb_freq_step'])
        self.press.clear()
        self.uncorrected.clear()
        self.enrichment.clear()
        self.qoi.clear()
        qoi = []
        if self.dataPb.get('nb_cpu', 1)>1:
            ## TODO: must be fixed
            with joblib.Parallel(n_jobs=self.dataPb.get('nb_cpu', 1), require=None) as parallel:
                nitf = len(freq_list)
                def run_one_freq_wrapper(itf, nitf, f):
                    logger.info('Solve freq: {:5g} Hz ({}/{})'.format(f, itf, nitf))
                    return self.run_one_freq(f)
                data = parallel( joblib.delayed(run_one_freq_wrapper)(itf,nitf,f) for (itf,f) in enumerate(freq_list) )
                qoi.append(data[1])
        else:
            for it,f in enumerate(freq_list):
                logger.info('Solve freq {}/{}: {:5g} Hz'.format(it+1,len(freq_list),f))
                data = self.run_one_freq(f)
                qoi.append(data[1])
        
        # export fields along frequencies
        self.export(kind='fields', show=False)
        if self.enrich:
            self.export(kind='fields_rebuilt', show=False)
        #
        format_QoI = np.reshape(np.array(qoi),shape=(len(qoi),len(qoi[0])))
        #
        id_str = 'FRF_'+misc_utils.generate_formatted_id(paranames=self.para_names, paraval=self.current_param_set)
        self.qoi.add(data=np.array(freq_list).reshape((len(freq_list),1)), name='freq_'+id_str)
        self.qoi.add(data=format_QoI, name=id_str)
        return {'press': self.press,
                'QoI': self.qoi}
        
    def pre_process_offline(self, meshing= True):
        # build mesh file
        lx = self.dataPb.get('lx')
        ly = self.dataPb.get('ly')
        lz = self.dataPb.get('lz')
        struct_mesh_size = self.dataPb.get('struct_mesh_size')
        if meshing:
            mesher_cube_tank.xfem_fluid_and_tank(
                lx,
                ly,
                lz,
                struct_mesh_size,
                self.fluid_order,
                self.mesh_file_fluid)
        
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
                                           1, # TODO: keep first order at this point ##self.shell_order,
                                           self.mesh_file_struct)
            # load struct
            self.load_struct()
            # build level set
            self.compute_LS()
        else:
            struct_thickness = self.dataPb['struct_thickness']
            mesher_cube_tank.classic_fluid_and_tank(lx,ly,lz,
                                                  struct_lx,
                                                  lx_up,
                                                  lx_down,
                                                  struct_lz,
                                                  struct_thickness,
                                                  struct_mesh_size,
                                                  self.fluid_order,
                                                  self.mesh_file_fluid)
        

        
    def run_one_freq(self, freq):
        #
        omega=2*np.pi*freq
        forceType ='float'
        # solve linear system 
        if self.load_type=='displacement':
            rhs = -omega**2*self.op['C']
        elif self.load_type=='velocity':
            rhs = omega*self.op['C']
        elif self.load_type=='acceleration':
            rhs = self.op['C']
        else:
            rhs = None
        #
        sol = solvers.solve_linear(self.method_sl, 
                            self.op['H']-omega**2*self.op['S'], 
                            rhs, comm=mycomm )
        if self.enrich:
            uncorrected      = np.zeros(self.fluid_ndof)
            uncorrected[self.SolvedDofF] = sol[self.SolvedDofF].copy()
            enrichment = np.zeros(self.fluid_ndof)
            enrichment[self.SolvedDofA]= sol[list(range(len(self.SolvedDofF),len(self.SolvedDofF)+len(self.SolvedDofA),1))].copy()
            CorrectedPressure=np.array(uncorrected)
            CorrectedPressure[self.SolvedDofA]=CorrectedPressure[self.SolvedDofA].T+np.array(enrichment[self.SolvedDofA]*np.sign(self.struct_LS[self.SolvedDofA]).T)
            press = CorrectedPressure
            #
            id_str = misc_utils.generate_formatted_id(freq=freq)
            self.press.add(data=press.copy(), name=id_str)
            self.enrichment.add(data=enrichment.copy(), name=id_str)
            self.uncorrected.add(data=uncorrected.copy(), name=id_str)
        else:
            id_str = misc_utils.generate_formatted_id(freq=freq)
            press = np.zeros(self.fluid_ndof)
            press[self.SolvedDofF] = sol[self.SolvedDofF].copy()
            self.press.add(data=press.copy(), name=id_str)
        # run post-processing
        QoI = self.post_process(freq,sol)

        
        return press,QoI
            
        
       