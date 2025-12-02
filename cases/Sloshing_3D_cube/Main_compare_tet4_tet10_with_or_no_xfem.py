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
dataPb['struct_lz'] = 0.52222333
dataPb['struct_thickness'] = 0.002 # only for classic conforming mesh
#size of elements
dataPb['struct_mesh_size'] =  dataPb['lx']/1 #20
# sets of parameters for structure position
para_val = np.array([[0.0,0.0],
                     [0.2,0.2],
                     [-0.2,0.2],
                     [0.2,-0.2],
                     [-0.2,-0.2]])
para_dict ={'names': ['lx_baffle_shift_up','lx_baffle_shift_down'],
            'values': para_val}


dataPb['freq_ini'] = 0.4
dataPb['freq_ref'] = 0.5
dataPb['freq_end'] = 1.5 
dataPb['nb_freq_step'] = 3#250

# Imposed acceleration on tank and stiffener surfaces 
dataPb['U_dot_dot_imposed'] = np.array([1.0e-3, 1.0e-3, 0.0])/np.sqrt(2.0)

# compute QoI 
epsbbx = 1e-10
dataPb['QoI_nodes_bbx'] = [0.-epsbbx ,0.+epsbbx, 0.8-epsbbx, 0.8+epsbbx, 0.6-epsbbx, 0.6+epsbbx]

# Flags
dataPb['flag_eigen_vectors'] = False
dataPb['flag_nb_eigen_modes'] = 10
dataPb['flag_FRF'] = True
dataPb['flag_write_gmsh_results'] = True
dataPb['flag_write_quadrature_points'] = True
dataPb['fluid_element'] = 'TET10'  # 'TET4' or 'TET10'
dataPb['shell_element'] = 'TRI6'  # 'TRI3' or 'TRI6'
dataPb['enrichment'] = True  # XFEM enrichment for structure or not

# fluid
dataFluid = dict()
#data['celerity'] = 343.0
dataFluid['rho'] = 1000.0
#data['fluid_damping'] = 1.0

files = dict()
files['fluid'] = Path(__file__).parent / 'xfem_tet10' / 'cube_xfem_sloshing_Fluid_and_Tank_tet10'
files['struct'] = Path(__file__).parent / 'xfem_tet10' / 'cube_xfem_sloshing_Stiffener_DKT'
files['results'] = Path(__file__).parent / 'xfem_tet10' / 'Main_compare_tet4_tet10_with_or_no_xfem_results'

objCompute = slib_sloshing.compute_sloshing(dataPb,dataFluid,files)
objCompute.run_parametric(para_dict)


# # parallepipedic cavity with plane structure
# mesh_file_fluid_tet10       =Path(__file__).parent / 'xfem_tet10' / 'cube_xfem_sloshing_Fluid_and_Tank_tet10'
# results_file_tet10          =Path(__file__).parent / 'xfem_tet10' / ('cube_xfem_sloshing_with_stiffener_tet10_'+id_formatted+'_h'+str(int(lx/h_fluid_elts)))
# mesh_file_fluid_tet10.parent.mkdir(parents=True, exist_ok=True)

# mesh_file_fluid_tet4        =Path(__file__).parent / 'xfem_tet4' / 'cube_xfem_sloshing_Fluid_and_Tank_tet4'
# results_file_tet4           =Path(__file__).parent / 'xfem_tet4' / ('cube_xfem_sloshing_with_stiffener_tet4_'+id_formatted+'_h'+str(int(lx/h_fluid_elts)))
# mesh_file_fluid_tet4.parent.mkdir(parents=True, exist_ok=True)

# mesh_file_stiffener         =Path(__file__).parent / 'cube_xfem_sloshing_Stiffener_DKT'

# mesh_file_tet10_classic     =Path(__file__).parent / 'classic_tet10' / 'cube_sloshing_with_stiffener_tet10'
# results_file_tet10_classic  =Path(__file__).parent / 'classic_tet10' / ('cube_sloshing_with_stiffener_tet10_'+id_formatted+'_h'+str(int(lx/h_fluid_elts)))
# mesh_file_tet10_classic.parent.mkdir(parents=True, exist_ok=True)
# #
# mesh_file_tet4_classic      = Path(__file__).parent / 'classic_tet4' / 'cube_sloshing_with_stiffener_tet4'
# results_file_tet4_classic   = Path(__file__).parent / 'classic_tet4' / ('cube_sloshing_with_stiffener_tet4_'+id_formatted+'_h'+str(int(lx/h_fluid_elts)))
# mesh_file_tet4_classic.parent.mkdir(parents=True, exist_ok=True)

# silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,2,mesh_file_fluid_tet10)
# silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,1,mesh_file_fluid_tet4)
# silex_lib_cube_tank_gmsh_geometry.Stiffener_DKT(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,h_fluid_elts*0.5,1,mesh_file_stiffener)
# silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,2,mesh_file_tet10_classic)
# silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,1,mesh_file_tet4_classic)





# # silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4_xfem(dataPb,dataFluid,mesh_file_fluid_tet4,mesh_file_stiffener,results_file_tet4)
# silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10_xfem(dataPb,dataFluid,mesh_file_fluid_tet10,mesh_file_stiffener,results_file_tet10)
# # silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4(dataPb,dataFluid,mesh_file_tet4_classic,results_file_tet4_classic)
# silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10(dataPb,dataFluid,mesh_file_tet10_classic,results_file_tet10_classic)
