
import time
import numpy as np
from typing import Union,List,Dict,AnyStr
from datetime import datetime

import scipy.sparse as sps
import scipy.sparse.linalg  as spla
from loguru import logger

import pylab as pl
import pickle
import h5py
import csv
import pyvista as pv

import joblib

try:
    from scikits import umfpack
except ImportError:
    umfpack = None

import sys
from pathlib import Path

import mumps
import gmsh
# import utils as u
# import utils_acoustics as ua
from meshRW.msh2 import mshWriter
from meshRW.msh import mshReader
from meshRW.vtk import vtkWriter

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

def create_dir_sym(path_file:Path):
    # get Path parts
    parts = list(path_file.parts)
    # update parent dirname
    parts[-2] = parts[-2] + datetime.now().strftime('_%Y-%m-%d_%H-%M-%S')
    new_path_file = Path(*parts)
    new_path_dir = Path(new_path_file.parent)
    if not new_path_dir.exists():
        new_path_dir.mkdir(parents=True, exist_ok=True)
    orig_path_dir = Path(path_file.parent)
    if orig_path_dir.is_symlink():
        orig_path_dir.unlink()
    orig_path_dir.symlink_to(new_path_dir, target_is_directory=True)
    return new_path_file 

class enginePV():
    def __init__(self, pl, gf, items):
        self.items = items
        self.gf = gf
        self.pl = pl
        self.output = None
            

    def __call__(self, it):
        self.current = int(it)
        self.update()
        
    def update(self):
        scalars_default = self.items[self.current]
        mesh = self.pl.add_mesh(self.gf,
                         scalars = scalars_default,
                         show_edges=True,
                         show_scalar_bar=True,
                         cmap='thermal')
        self.output = mesh

class exportPV:
    def __init__(self, 
                 nodes, 
                 elements, 
                 nodalFields = None, 
                 elementalFields = None,
                 vectorFields = None,
                 title='',
                 show=True):
        self.items = []
        self.grid = None
        self.grid_fields = []
        self.arrows = None
        self.plotter = pv.Plotter()
        self.subPlotsId = []
        # load mesh
        self.loadMesh(nodes,elements)
        # load fields
        if nodalFields is not None:
            self.loadNodalFields(nodalFields)
        if elementalFields is not None:
            self.loadElementFields(elementalFields)
        if vectorFields is not None:
            self.loadVectorFields(vectorFields)
        if show:
            self.show()
        
    
    def getVTKelements(self, txt):
        db = {'QUA4':pv.CellType.QUAD,
              'TRI3': pv.CellType.TRIANGLE,
              'TRI6': pv.CellType.QUADRATIC_TRIANGLE,
              'TET4': pv.CellType.TETRA,
              'TET10': pv.CellType.QUADRATIC_TETRA,
              'PRI6': pv.CellType.WEDGE}
        return db.get(txt)
    
    def getNbPoints(self, txt):
        db = {'QUA4':4,
              'TRI3': 3,
              'TRI6': 6,
              'TET4': 4,
              'TET10': 10,
              'PRI6' : 6}
        return db.get(txt)
        
    def loadMesh(self, nodes, elements):
        if not isinstance(elements, list):
            elements = [elements]
        cell_list = []
        cell_list_types = []
        for e in elements:
            nbPts = self.getNbPoints(e.get('type'))
            connectivity = e.get('connectivity')-1
            cell_list.append(np.hstack((nbPts*np.ones((connectivity.shape[0],1),dtype=int),connectivity)))
            cell_list_types.append(connectivity.shape[0]*[self.getVTKelements(e.get('type'))])
        cells = np.hstack([c.flatten() for c in cell_list])
        cell_types = np.hstack(cell_list_types)
        m = pv.UnstructuredGrid(cells, 
                                cell_types,
                                nodes)
        self.grid = m
        
    def loadNodalFields(self, fields):
        if isinstance(fields, dict):
            fields = [fields]
        for f in fields:
            self.grid_fields.append(self.grid)
            self.items.append([])
            nbSteps = f.get('nbsteps')
            if nbSteps:
                for i in range(nbSteps):
                    txt = f.get('name')+f'-{i:05d}'
                    self.grid_fields[-1].point_data[txt] = f.get('data')[:,i]
                    self.items[-1].append(txt)
            else:
                txt = f.get('name')
                self.grid_fields[-1].point_data[txt] = f.get('data')
                self.items[-1].append(txt)
                
    def loadVectorFields(self, field):
        txt = field.get('name')
        self.grid_fields.append(self.grid)
        self.grid_fields[-1].cell_data[txt] = field.get('data')
        self.items.append([txt])
        arrows = self.grid_fields[-1].glyph(
            orient = txt,
            scale = txt,
            factor = 0.5
        )
        self.arrows = arrows
        
    def loadElementFields(self, fields):
        if isinstance(fields, dict):
            fields = [fields]
        for f in fields:
            self.grid_fields.append(self.grid)
            self.items.append([])
            nbSteps = f.get('nbsteps')
            if nbSteps:
                for i in range(nbSteps):
                    txt = f.get('name')+f'-{i:05d}'
                    self.grid_fields[-1].cell_data[txt] = f.get('data')[:,i]
                    self.items[-1].append(txt)
            else:
                txt = f.get('name')
                self.grid_fields[-1].cell_data[txt] = f.get('data')
                self.items[-1].append(txt)
                
    def getShapeSubplots(self, nbFields):
        self.currentSubplot = (1,1)
        if nbFields==0:
            nbFields = 1
        if nbFields<=3:
            shapeSubPlots = (1,nbFields)
        elif nbFields<=6:
            shapeSubPlots = (2,int(np.ceil(nbFields/2)))
        elif nbFields<=9:
            shapeSubPlots = (3,int(np.ceil(nbFields/3)))
        else:
            shapeSubPlots = (4,int(np.ceil(nbFields/4)))
        self.subPlotsId = []
        for i in range(shapeSubPlots[0]):
            for j in range(shapeSubPlots[1]):
                self.subPlotsId.append( (i,j) )
        return shapeSubPlots
                
        
    def iterateSubplots(self, it):
        return self.subPlotsId[it]
                
    def show(self, slider=True):
        nbFields = len(self.grid_fields)
        if self.arrows is not None:
            nbFields += 1
        self.shapeSubplot = self.getShapeSubplots(nbFields)
        pl = pv.Plotter(shape=self.shapeSubplot)
        itSubplot = 0
        # 
        if self.arrows is not None:
            pl.subplot(*self.iterateSubplots(itSubplot))
            itSubplot += 1
            pl.add_mesh(self.grid,
                        show_edges=True,
                        show_scalar_bar=True,
                        cmap='thermal')
            pl.add_mesh(self.arrows, color='red')
        for ig,gf in enumerate(self.grid_fields):
            pl.subplot(*self.iterateSubplots(itSubplot))
            itSubplot += 1
            scalars_default=None
            title = ''
            if slider:
                engine = enginePV(pl, gf, self.items[ig])                   
                pl.add_slider_widget(engine,[0,len(self.items[ig])-1],title='It')
            else:
                if len(self.items[ig])>0:
                    scalars_default = self.items[ig][0]
                    title = self.items[ig][0]
                pl.add_mesh(gf,
                            scalars = scalars_default,
                            show_edges=True,
                            show_scalar_bar=True,
                            cmap='thermal')
            pl.add_title(title)
        if itSubplot==0:
            pl.add_mesh(self.grid,
                        show_edges=True,
                        show_scalar_bar=True,
                        cmap='thermal')
        pl.show()
        

def solve_linear(method, A, b, comm=None):
    if method == 'mumps':
        x = mumps.spsolve(A, b, comm=comm)
    elif method == 'scipy':
        x = spla.spsolve(A, b)
    elif method == 'umfpack':
        x = umfpack.spsolve(A, b)        
    return x

class DB_results:
    def __init__(self):
        self.db = []
        self.names = []
    def clear(self):
        self.db = []
        self.names = []
    @property
    def nb_fields(self):
        return len(self.db)
    @property
    def db_format(self):
        if len(self.db[0].shape)==1:
            return np.vstack(self.db).T
        else:
            return np.hstack(self.db)
    def add(self, name, data):
        self.names.append(name)
        self.db.append(data)
    def export_csv(self, filename):
        header = ','.join(self.names)
        np.savetxt(filename, self.db_format, delimiter=',', header=header)
    
    def export_hdf5(self, filename):
        # write HDF5 file
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset('db',
                                data=self.db_format, 
                                chunks=True, 
                                compression="gzip",
                                fletcher32=True)
    
    def export(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump({'names': self.names, 'db': self.db}, f)


class compute_sloshing():
    def __init__(self, dataPb, dataFluid, dataIO):        
        self.dataPb = dataPb
        self.dataFluid = dataFluid
        self.mesh_file_struct = dataIO.get('struct', None)
        self.mesh_file_fluid = dataIO.get('fluid', None)
        self.results_file = dataIO.get('results', None)
        self.io = dataIO
        self.results_file_basis = self.results_file
        #
        self._SolvedDofF = []
        self._SolvedDofA = []
        self.fluid_nodes = []
        self.fluid_elements = []
        self.struct_nodes = []
        self.struct_elements = []
        self.struct_id_elements = []
        # 
        self.op = dict()
        #
        self.para_names = []
        self.current_para_set = []
        self.press = DB_results()
        self.disp = DB_results()
        self.enrichment = DB_results()
        self.uncorrected = DB_results()
        self.qoi = DB_results()
        
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
    
    def get_tag(self, name):
        db = self.tag_db()
        return db.get(name)
    
    @property
    def SolvedDofF(self):
        # if len(self._SolvedDofF) == 0:
        #     self._SolvedDofF = list(range(self.fluid_ndof))
        return list(range(self.fluid_ndof))
    @property
    def SolvedDofA(self):
        return self.enriched_nodes
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
    def callExport(self):
        formatExport = self.io.get('format', 'msh')
        if formatExport == 'msh':
            return mshWriter
        elif formatExport == 'vtk':
            return vtkWriter
        else:
            raise ValueError('Unsupported export format')
    @property
    def exportExtension(self):
        formatExport = self.io.get('format', 'msh')
        if formatExport == 'msh':
            return '.msh'
        elif formatExport == 'vtk':
            return '.vtu'
        else:
            raise ValueError('Unsupported export format')    
    
    
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
        objMeshFluid = mshReader(self.mesh_file_fluid.as_posix()+'.msh')
        self.fluid_nodes = objMeshFluid.getNodes()
        self.fluid_elements = objMeshFluid.getElements(tag=self.get_tag('tag_fluid_volume'),
                                                       type=self.dataPb.get('fluid_element'),
                                                       dictFormat=False)
        #
        self.free_fluid_elements = objMeshFluid.getElements(tag=self.get_tag('tag_fluid_free_surface'), dictFormat=False)
        self.fluid_bounds_elements = objMeshFluid.getElements(tag=self.get_tag('tag_fluid_structure_surface'), dictFormat=False)
        #
    
    @utils.timeit('Load structure mesh')
    def load_struct(self):
        objMeshStruct = mshReader(self.mesh_file_struct.as_posix()+'.msh')
        self.struct_nodes = objMeshStruct.getNodes()
        self.struct_elements = objMeshStruct.getElements(tag=self.get_tag('tag_struct_surface'), dictFormat=False)
        self.struct_edge_elements = objMeshStruct.getElements(tag=self.get_tag('tag_struct_edge'), dictFormat=False)
        pass
        
    def show_data(self):
        logger.info('Nb fluid nodes: {}'.format(len(self.fluid_nodes)))
        logger.info('Nb fluid elements: {}'.format(len(self.fluid_elements)))
        logger.info('Nb free fluid elements: {}'.format(len(self.free_fluid_elements)))
        logger.info('Nb fluid boundary elements: {}'.format(len(self.fluid_bounds_elements)))
        if self.enrich and len(self.struct_nodes)>0:
            logger.info('Nb structure nodes: {}'.format(len(self.struct_nodes)))
            logger.info('Nb structure elements: {}'.format(len(self.struct_elements)))
        pass
    
    def export(self,kind: Union[AnyStr,List[AnyStr]]='fluid', show=False):
        if isinstance(kind,str):
            kind=[kind]
        for k in kind:
            if k=='fluid':
                elements = {'connectivity': self.fluid_elements,
                           'type': self.dataPb.get('fluid_element')}
                title = f'Mesh of fluid volume ({self.dataPb.get("fluid_element")})'
                self.callExport(self.results_file.as_posix()+'_Mesh_Fluid_volume'+self.exportExtension,
                           nodes=self.fluid_nodes,
                           title=title,
                           elements=elements)
                if show:
                    exportPV(self.fluid_nodes,elements,show=True)
                
                elements = {'connectivity': self.fluid_bounds_elements,
                            'type': self.dataPb.get('shell_element')}
                title = f'Mesh of tank surfaces ({self.dataPb.get("shell_element")})'
                self.callExport(self.results_file.as_posix()+'_Mesh_Tank_surfaces'+self.exportExtension,
                           nodes=self.fluid_nodes,
                           title=title,
                           elements=elements)
                if show:
                    exportPV(self.fluid_nodes,elements,show=False)
                
            elif k=='free_surface':
                self.callExport(self.results_file.as_posix()+'_Mesh_Fluid_Free_surface'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                elements={'connectivity': self.free_fluid_elements,
                                          'type': self.dataPb.get('shell_element')},
                                title=f'Mesh of free fluid surface ({self.dataPb.get("shell_element")})')
            elif k=='structure':
                self.callExport(self.results_file.as_posix()+'_Mesh_Stiffener_surface'+self.exportExtension,
                                nodes=self.struct_nodes,
                                elements={'connectivity': self.structure_elements,
                                          'type': self.dataPb.get('shell_element')},
                                title=f'Mesh of structure ({self.dataPb.get("shell_element")})')
                
                # silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Stiffener_surface',
                #                             self.struct_nodes,
                #                             self.structure_elements,2)
            elif k=='levelset':
                dataW = list()
                dataW.append({'data':self.struct_LS,'type':'nodal', 'name':'Tevelset'})
                dataW.append({'data':self.struct_LS_tangent,'type':'nodal','name':'Tangent levelset'})
                self.callExport(self.results_file.as_posix()+'_LevelSet'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                elements={'connectivity': self.fluid_elements,
                                          'type': self.dataPb.get('fluid_element')},
                                title=f'Levelsets ({self.dataPb.get("fluid_element")})',
                                fields=dataW, append=True)
                self.callExport(self.results_file.as_posix()+'_Enriched_Fluid_Elements'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                elements={'connectivity': self.fluid_elements[self.enriched_elements-1],
                                          'type': self.dataPb.get('fluid_element')},
                                title=f'Levelsets ({self.dataPb.get("fluid_element")})')
            elif k=='load':
                elements = {'connectivity': self.fluid_bounds_elements,
                                          'type': self.dataPb.get('shell_element')}
                title = f'Mesh of structure ({self.dataPb.get("shell_element")})'
                fields = {'data': self.vecNormalEltsF,'type':'elemental', 'name':'Normal to tank elements'}
                self.callExport(self.results_file.as_posix()+'_Mesh_Normal_to_tank_surfaces'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                elements=elements,
                                title=title,
                                fields=fields,
                                append=True)
                if show:
                    exportPV(self.fluid_nodes,elements,elementalFields=fields,show=True)
            elif k=='load_LS':
                elements = {'connectivity': self.struct_elements,
                                          'type': self.dataPb.get('shell_element')}
                title = f'Mesh of structure ({self.dataPb.get("shell_element")})'
                fields = {'data': self.vecNormalEltsA,'type':'elemental', 'name':'Normal to stiffener elements'}
                self.callExport(self.results_file.as_posix()+'_Mesh_Normal_to_stiffener'+self.exportExtension,
                                nodes=self.struct_nodes,
                                elements=elements,
                                title=title,
                                fields=fields,
                                append=True)
                if show:
                    exportPV(self.struct_nodes,elements,vectorFields=fields,show=True)

            elif k=='TET10toTET4':
                self.callExport(self.results_file.as_posix()+'_Mesh_Fluid_volume_tet10TOtet4'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                title='TET4 fluid mesh converted from TET10',
                                elements={'connectivity': self.dataTET10toTET4['elements'],
                                          'type': 'TET4'})
                self.callExport(self.results_file.as_posix()+'_Enriched_Fluid_Elements_tet10TOtet4'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                title='TET4 enriched fluid elements converted from TET10',
                                elements={'connectivity': self.dataTET10toTET4['elements'][self.dataTET10toTET4['enriched_elements']-1],
                                          'type': 'TET4'})
                # silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Mesh_Fluid_volume_tet10TOtet4',
                #                             self.fluid_nodes,
                #                             self.dataTET10toTET4['elements'],4)
                # silex_lib_gmsh.WriteResults(self.results_file.as_posix()+'_Enriched_Fluid_Elements_tet10TOtet4',
                #                             self.fluid_nodes,
                #                             self.dataTET10toTET4['elements'][self.dataTET10toTET4['enriched_elements']-1],4)
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
                # write HDF5 file
                with h5py.File(self.results_file.as_posix()+'_eigen_frequencies.h5', 'w') as hf:
                    hf.create_dataset('eigen_frequencies',
                                      data=self.eigen_frequencies, 
                                      chunks=True, 
                                      compression="gzip",
                                      fletcher32=True)
                    
            elif k=='frf':
                self.qoi.export_hdf5(self.results_file.as_posix()+'_frf.h5')
                self.qoi.export(self.results_file.as_posix()+'.frf')

            elif k=='qoi':                
                pass
            
            elif k=='fields':
                dataW = list()
                dataW.append({'data':self.press.db_format, 
                              'nbsteps':self.press.nb_fields, 
                              'type':'nodal', 
                              'name':'Pressure' })
                if self.enrich:
                    dataW.append({'data':self.uncorrected.db_format, 
                              'nbsteps':self.uncorrected.nb_fields,  
                              'type':'nodal', 
                              'name':'Uncorrected Pressure' })
                    dataW.append({'data':self.enrichment.db_format, 
                              'nbsteps':self.enrichment.nb_fields,  
                              'type':'nodal', 
                              'name':'Enrichement' })
                elements = {'connectivity': self.fluid_elements,
                                          'type':self.dataPb.get("fluid_element")}
                self.callExport(self.results_file.as_posix()+'_fields'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                title='Fields',
                                append=True,
                                elements=elements,
                                fields=dataW)
                if show:
                    exportPV(self.fluid_nodes,elements,nodalFields=dataW,show=True)
            elif k=='fields_rebuilt':        
                objMesh = lib.MeshField(nodes=self.fluid_nodes, 
                                            elems=self.fluid_elements, 
                                            leveset=self.struct_LS, 
                                            levelsetTg=self.struct_LS_tangent)
                objMesh.addField(uncorrectedField=self.uncorrected.db_format,
                                enrichmentField=self.enrichment.db_format)
                datameshfield = objMesh.getData()
            
                # prepare fields
                dataW = []
                dataW.append({
                    'name': 'press',
                    'nbsteps': self.press.nb_fields,
                    'type': 'nodal',
                    'data': datameshfield['fields']
                })
                dataW.append({
                    'name': 'uncorrected press',
                    'nbsteps': self.press.nb_fields,
                    'type': 'nodal',
                    'data': datameshfield['uncorrected']
                })
                dataW.append({
                    'name': 'correction press',
                    'nbsteps': self.press.nb_fields,
                    'type': 'nodal',
                    'data': datameshfield['correction']
                })
                # elements
                elements = [{'type': 'TET4', 'connectivity': datameshfield['TET4']},
                              {'type': 'PRI6', 'connectivity': datameshfield['PRI6']}]
                # export mesh
                self.callExport(
                    filename= self.results_file.as_posix()+'_results_fluid_frf_meshfield3D'+self.exportExtension,
                    nodes=datameshfield['nodes'],
                    elements=elements,
                    fields=dataW,
                    append=True
                    )
                
                if show:
                    exportPV(datameshfield['nodes'],elements,nodalFields=dataW, show=True)
                
            elif k=='eigenmodes':
                self.callExport(self.results_file.as_posix()+'_Eigen_modes'+self.exportExtension,
                                nodes=self.fluid_nodes,
                                title='Eigen modes',
                                elements={'connectivity': self.fluid_elements,
                                          'type':self.dataPb.get("fluid_element")},
                                fields=[{'data': self.eigen_vectors, 'nbsteps': self.eigen_vectors.shape[1], 'type':'nodal', 'name':'Modes'},
                                        {'data': self.eigen_vectors_uncorrected, 'nbsteps': self.eigen_vectors.shape[1], 'type':'nodal', 'name':'Modes w/o correction'},
                                        {'data': self.eigen_vectors_enrichment, 'nbsteps': self.eigen_vectors.shape[1], 'type':'nodal', 'name':'Modes enrichment'},],
                                append=True)
                if len(self.dataTET10toTET4['elements'])>0:
                    self.callExport(self.results_file.as_posix()+'_results_fluid_eigenmodes_on_tet4mesh'+self.exportExtension,
                                    nodes=self.fluid_nodes,
                                    title='Eigen modes on TET4 mesh frome TET10',
                                    elements={'connectivity': self.dataTET10toTET4['elements'],
                                            'type': 'TET4'},
                                    fields ={'data': self.eigen_vectors, 'nbsteps': self.eigen_vectors.shape[1], 'type':'nodal', 'name':'Modes (pressure)'},
                                    append=True)
                
                
                # silex_lib_gmsh.WriteResults2(self.results_file.as_posix()+'_Eigen_modes',
                #                              self.fluid_nodes,
                #                              self.fluid_elements,
                #                              11,
                #                              [[self.eigen_vectors,'nodal',1,'modes'],
                #                               [self.eigen_vectors_uncorrected,'nodal',1,'press classic'],
                #                               [self.eigen_vectors_enrichment,'nodal',1,'press enrich']])
                # silex_lib_gmsh.WriteResults2(self.results_file.as_posix() +'_results_fluid_eigenmodes_on_tet4mesh',
                #                              self.fluid_nodes,
                #                              self.dataTET10toTET4['elements'],
                #                              4,
                #                              [[self.eigen_vectors,'nodal',1,'pressure']]
                #                              )

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
                self.callExport(
                    filename= self.results_file.as_posix()+'_results_fluid_eigenmodes_meshfield3D'+self.exportExtension,
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
        self.enriched_nodes = np.unique(self.fluid_elements[self.enriched_elements-1])-1
    
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
    
    def compute_operators_online(self):        
        if self.enrich:
            self.compute_xfem_op_online()
        else:
            self.compute_operators()
    
    
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
        
    def compute_loads_online(self):
        if self.enrich:
            self.compute_xfem_loads()
        else:
            self.compute_loads()
        
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
        
    def generate_formatted_id(self, paraval=None, paranames=None, freq=None):
        # generate an id string based on parameter values
        id_formatted = ''
        if paraval is not None and paranames is not None:
            for i,pname in enumerate(paranames): 
                parastr = f'{int(1e3*paraval[i]):+04d}'
                id_formatted += '{}_{}'.format(pname, parastr.replace('+','p').replace('-','m'))
                if i<len(paranames)-1:
                    id_formatted += '_'
        if freq is not None:
            freq_str = f'{freq:5g}'
            if id_formatted != '' and not id_formatted.endswith('_'):
                id_formatted += '_'
            id_formatted += 'freq_{}Hz'.format(freq_str)
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
        self.para_names = para_names
        logger.info('Run parametric study for {} parameters and {} sets'.format(para_val.shape[1], para_val.shape[0]))
        # run along each parameter set
        results = []
        for i,pset in enumerate(para_val):
            logger.info('Run parametric set {}/{}: {}'.format(i+1, para_val.shape[0], pset))
            # save current parameter set and value
            self.current_param_set = pset            
            # update results file names
            id_formatted = self.generate_formatted_id(pset, para_names)
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
        id_str = 'FRF_'+self.generate_formatted_id(paranames=self.para_names, paraval=self.current_param_set)
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
        if self.load_type=='displacement':
            rhs = -omega**2*self.op['C']
        elif self.load_type=='velocity':
            rhs = omega*self.op['C']
        elif self.load_type=='acceleration':
            rhs = self.op['C']
        else:
            rhs = None
        #
        sol = solve_linear(self.method_sl, 
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
            id_str = self.generate_formatted_id(freq=freq)
            self.press.add(data=press.copy(), name=id_str)
            self.enrichment.add(data=enrichment.copy(), name=id_str)
            self.uncorrected.add(data=uncorrected.copy(), name=id_str)
        else:
            id_str = self.generate_formatted_id(freq=freq)
            press = np.zeros(self.fluid_ndof)
            press[self.SolvedDofF] = sol[self.SolvedDofF].copy()
            self.press.add(data=press.copy(), name=id_str)
        # run post-processing
        QoI = self.post_process(freq,sol)

        
        return press,QoI
            
        
       