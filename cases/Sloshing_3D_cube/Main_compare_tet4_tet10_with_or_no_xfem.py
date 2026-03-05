import numpy as np


from pathlib import Path
# import utils as u
from  SILEXlib.tests import utils_acoustics as ua

import silex_lib_compute_sloshing as slib_sloshing
import silex_lib_cube_tank_gmsh_geometry as slib_geometry

#
##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################
dataPb = dict()

# gravity
dataPb['g'] = 9.81 # m/s2

# cube geom
dataPb['lx'] = 1.0
dataPb['ly'] = 0.8
dataPb['lz'] = 0.6

# Baffle position and geom
dataPb['struct_lx'] = 0.53333333333
dataPb['struct_lz'] = 0.32222333
dataPb['struct_thickness'] = 0.002 # only for classic conforming mesh
#size of elements
crit_mesh = 50 
dataPb['struct_mesh_size'] =  dataPb['lx']/crit_mesh #20
# sets of parameters for structure position
para_val = np.array([[+0.0,+0.0],
                     [+0.2,+0.2],
                     [-0.2,+0.2],
                     [+0.2,-0.2],
                     [-0.2,-0.2]])
para_dict ={'names': ['lx_baffle_shift_up','lx_baffle_shift_down'],
            'values': para_val}


dataPb['freq_ini'] = 0.4
dataPb['freq_ref'] = 0.5
dataPb['freq_end'] = 1.5 
dataPb['nb_freq_step'] = 10
dataPb['nb_cpu'] = 1 # set to 1 to deactivate parallelism in the frequency loop and thus be able to compare results with or without parallelism

# Imposed acceleration on tank and stiffener surfaces 
# dataPb['U_dot_dot_imposed'] = np.array([1.0e-3, 1.0e-3, 0.0])/np.sqrt(2.0)
# dataPb['U_dot_imposed'] = np.array([1.0e-3, 1.0e-3, 0.0])/np.sqrt(2.0)
dataPb['U_imposed'] = np.array([1.0e-3, 1.0e-3, 0.0])/np.sqrt(2.0)

# compute QoI 
epsbbx = 1e-10
dataPb['QoI_nodes_bbx'] = [0.-epsbbx ,0.+epsbbx, 0.8-epsbbx, 0.8+epsbbx, 0.6-epsbbx, 0.6+epsbbx]

# Flags
dataPb['flag_eigen_vectors'] = False
dataPb['flag_nb_eigen_modes'] = 10
dataPb['flag_FRF'] = True
dataPb['flag_write_gmsh_results'] = True
dataPb['flag_write_quadrature_points'] = True
dataPb['fluid_element'] = 'TET4'  # 'TET4' or 'TET10'
dataPb['shell_element'] = 'TRI3'  # 'TRI3' or 'TRI6'
dataPb['enrichment'] = True  # XFEM enrichment for structure or not

# fluid
dataFluid = dict()
#data['celerity'] = 343.0
dataFluid['rho'] = 1000.0
#data['fluid_damping'] = 1.0

################################################################
################################################################
import_export = dict()
dataPb['enrichment'] = True
import_export['fluid'] = Path(__file__).parent / 'xfem_tet4' / ('cube_xfem_sloshing_Fluid_and_Tank_tet4_h'+str(crit_mesh))
import_export['struct'] = Path(__file__).parent / 'xfem_tet4' / ('cube_sloshing_Stiffener_DKT_h'+str(crit_mesh))
import_export['results'] = Path(__file__).parent / 'xfem_tet4' / ('cube_sloshing_results_h'+str(crit_mesh))
import_export['format'] = 'msh'

slib_sloshing.create_dir_sym(import_export['fluid'])
objCompute = slib_sloshing.compute_sloshing(dataPb,dataFluid,import_export)
objCompute.run_parametric(para_dict)


################################################################
################################################################
import_export = dict()
dataPb['enrichment'] = False
import_export['fluid'] = Path(__file__).parent / 'classic_tet4' / ('cube_sloshing_Fluid_and_Tank_tet4_h'+str(crit_mesh))
import_export['struct'] = None
import_export['results'] = Path(__file__).parent / 'classic_tet4' / ('cube_sloshing_results_h'+str(crit_mesh))
import_export['format'] = 'msh'

slib_sloshing.create_dir_sym(import_export['fluid'])
objCompute = slib_sloshing.compute_sloshing(dataPb,dataFluid,import_export)
objCompute.run_parametric(para_dict)


################################################################
################################################################
dataPb['fluid_element'] = 'TET10'  # 'TET4' or 'TET10'
dataPb['shell_element'] = 'TRI6'  # 'TRI3' or 'TRI6'
dataPb['enrichment'] = True
import_export = dict()
import_export['fluid'] = Path(__file__).parent / 'xfem_tet10' / ('cube_xfem_sloshing_Fluid_and_Tank_tet10_h'+str(crit_mesh))
import_export['struct'] = Path(__file__).parent / 'xfem_tet10' / ('cube_sloshing_Stiffener_DKT_h'+str(crit_mesh))
import_export['results'] = Path(__file__).parent / 'xfem_tet10' / ('cube_sloshing_results_h'+str(crit_mesh))
import_export['format'] = 'msh'

slib_sloshing.create_dir_sym(import_export['fluid'])
objCompute = slib_sloshing.compute_sloshing(dataPb,dataFluid,import_export)
objCompute.run_parametric(para_dict)


################################################################
################################################################
import_export = dict()
dataPb['enrichment'] = False
import_export['fluid'] = Path(__file__).parent / 'classic_tet10' / ('cube_sloshing_Fluid_and_Tank_tet4_h'+str(crit_mesh))
import_export['struct'] = None
import_export['results'] = Path(__file__).parent / 'classic_tet10' / ('cube_sloshing_results_h'+str(crit_mesh))
import_export['format'] = 'msh'

slib_sloshing.create_dir_sym(import_export['fluid'])
objCompute = slib_sloshing.compute_sloshing(dataPb,dataFluid,import_export)
objCompute.run_parametric(para_dict)


