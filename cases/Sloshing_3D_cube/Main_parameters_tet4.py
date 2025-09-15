import string
import time
import numpy as np


from loguru import logger as log

import pylab as pl
import pickle

import sys
from pathlib import Path

#sys.path.append("/home/legay/Codes/SILEXGIT/SILEXlib/SILEXlib/tests/")


import silex_lib_compute_sloshing
import silex_lib_cube_tank_gmsh_geometry


from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank = comm.Get_rank()


class comm_mumps_one_proc:
    rank = 0

    def py2f(self):
        return 0


mycomm = comm_mumps_one_proc()

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
class solve:
    def __init__(self):
        self.dataPb = dict()

        # cube geom
        self.lx = 1.0
        self.ly = 0.8
        self.lz = 0.6

        # Baffle position and geom
        self.lx_baffle = 0.53333333333
        self.lz_baffle = 0.422223333
        self.thickness_baffle = 0.002  # only for classic conforming mesh
        self.lx_baffle_shift_up = 0.18
        self.lx_baffle_shift_down = 0.08
        
        

        self.dataPb["freq_ini"] = 1.1498 #0.9
        self.dataPb["freq_ref"] = 0.5
        self.dataPb["freq_end"] = 1.15 #1.1
        self.dataPb["nb_freq_step"] = 5

        # Imposed acceleration on tank and stiffener surfaces
        self.dataPb["U_dot_dot_imposed"] = np.array([1.0e-3, 1.0e-3, 0.0])/np.sqrt(2.0)  

        # Flags
        self.dataPb["flag_eigen_vectors"] = 0
        self.dataPb["flag_FRF"] = 1
        self.dataPb["flag_write_gmsh_results"] = 0

        # fluid
        self.dataFluid = dict()
        # data['celerity'] = 343.0
        self.dataFluid["rho"] = 1000.0
        # data['fluid_damping'] = 1.0


        # size of elements
        self.h_fluid_elts = self.lx / 23

        # # parallepipedic cavity with plane structure
        # mesh_file_fluid_tet10       =Path(__file__).parent / 'cube_xfem_sloshing_Fluid_and_Tank_tet10'
        # results_file_tet10          =Path(__file__).parent / 'cube_xfem_sloshing_with_stiffener_tet10_h20'

        self.mesh_file_fluid_tet4 = Path(__file__).parent / "cube_xfem_sloshing_Fluid_and_Tank_tet4"
        self.results_file_tet4 = Path(__file__).parent / "cube_xfem_sloshing_with_stiffener_tet4"
        self.mesh_file_stiffener = Path(__file__).parent / "cube_xfem_sloshing_Stiffener_DKT"

        # mesh_file_tet10_classic     =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10'
        # results_file_tet10_classic  =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet10_h20'

        # mesh_file_tet4_classic      =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4'
        # results_file_tet4_classic   =Path(__file__).parent / 'cube_sloshing_with_stiffener_tet4_h20'

        # silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(lx,ly,lz,h_fluid_elts,2,mesh_file_fluid_tet10)
        silex_lib_cube_tank_gmsh_geometry.xfem_fluid_and_tank(self.lx, 
                                                              self.ly, 
                                                              self.lz, 
                                                              self.h_fluid_elts, 
                                                              1, 
                                                              self.mesh_file_fluid_tet4)
        
    def run(self, parameters):
        
        # default values
        lx_baffle = self.lx_baffle
        lx_baffle_shift_up = self.lx_baffle_shift_up
        lx_baffle_shift_down = self.lx_baffle_shift_down
        lz_baffle = self.lz_baffle
        h_baffle = self.h_fluid_elts * 0.5
        
        lx_baffle_shift_up = parameters[0]
        if len(parameters)>1:
            lx_baffle_shift_down = parameters[1]
        if len(parameters)>2:
            lx_baffle = parameters[2]
        if len(parameters)>3:
            lz_baffle = parameters[3]
        if len(parameters)>4:
            h_baffle = parameters[4]
        
        silex_lib_cube_tank_gmsh_geometry.Stiffener_DKT(
            self.lx,
            self.ly,
            self.lz,
            lx_baffle,
            lx_baffle_shift_up,
            lx_baffle_shift_down,
            lz_baffle,
            h_baffle,
            1,
            self.mesh_file_stiffener,
        )
        # silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,2,mesh_file_tet10_classic)
        # silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,1,mesh_file_tet4_classic)



        return silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4_xfem(self.dataPb, 
                                                                          self.dataFluid,
                                                                          self.mesh_file_fluid_tet4,
                                                                          self.mesh_file_stiffener, 
                                                                          self.results_file_tet4)
        # silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10_xfem(dataPb,dataFluid,mesh_file_fluid_tet10,mesh_file_stiffener,results_file_tet10)
        # silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4(dataPb,dataFluid,mesh_file_tet4_classic,results_file_tet4_classic)
        # silex_lib_compute_sloshing.sloshing_rigid_baffle_tet10(dataPb,dataFluid,mesh_file_tet10_classic,results_file_tet10_classic)

    def run_bis(self, parameters):
        
        # default values
        lx_baffle = self.lx_baffle
        ly_baffle = self.ly
        lz_baffle = self.lz_baffle
        lx_baffle_shift_up = self.lx_baffle_shift_up
        lx_baffle_shift_down = self.lx_baffle_shift_down
        lz_baffle = self.lz_baffle
        h_baffle = self.h_fluid_elts * 0.5
        
        if len(parameters)<6:
            log.error("You must provide exactly 6 parameters")
            raise ValueError("You must provide exactly 6 parameters")
        
        lxashift_up_A = parameters[0]
        lxashift_down_A = parameters[1]
        lxashift_up_B = parameters[2]
        lxashift_down_B = parameters[3]
        lzashift_up_A = 0.0
        lzashift_up_B = 0.0
        if len(parameters)>4:
            lzashift_up_A = parameters[4]
        if len(parameters)>5:
            lzashift_up_B = parameters[5]
        
        silex_lib_cube_tank_gmsh_geometry.Stiffener_DKT_bis(
            lx_baffle,
            ly_baffle,
            lz_baffle,
            lxashift_up_A,
            lxashift_down_A,
            lxashift_up_B,
            lxashift_down_B,
            lzashift_up_A,
            lzashift_up_B,
            h_baffle,
            1,
            self.mesh_file_stiffener,
        )
        # silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,2,mesh_file_tet10_classic)
        # silex_lib_cube_tank_gmsh_geometry.classic_fluid_and_tank(lx,ly,lz,lx_baffle,lx_baffle_shift_up,lx_baffle_shift_down,lz_baffle,thickness_baffle,h_fluid_elts,1,mesh_file_tet4_classic)



        return silex_lib_compute_sloshing.sloshing_rigid_baffle_tet4_xfem(self.dataPb, 
                                                                          self.dataFluid,
                                                                          self.mesh_file_fluid_tet4,
                                                                          self.mesh_file_stiffener, 
                                                                          self.results_file_tet4)

if __name__ == "__main__":

    obj = solve()
    # valtest = obj.run([0.09666666666666668, -0.1933333333333333])  # Run once to initialize the mesh and results
    # log.info(f"Test run completed with value: {valtest}")

    # X = np.linspace()
    nb_val = 11
    lup = np.linspace(-0.21, 0.21, nb_val)
    ldown = np.linspace(-0.21, 0.21, nb_val)
    X, Y = np.meshgrid(lup, ldown)
    val_p = np.zeros(nb_val * nb_val)
    val_f = np.zeros(nb_val * nb_val)

    for i,(xs,ys) in enumerate(zip(X.flatten(), Y.flatten())):
        log.info(f"Running for {xs:.2f} and {ys:.2f}")
        # return de run() : meanQI, maxQI, meanforce, maxforce
        val_p[i],_, val_f[i],_ = obj.run([xs, ys])
        log.info(f"Results: {val_p[i]}Pa and {val_f[i]}N")
        if val_p[i] > 1e12 or val_f[i] > 1e12:
            log.error(f"Error in computation for parameters {xs}, {ys}")
            raise ValueError(f"Computation failed for parameters {xs}, {ys}")
        


        
    results_file = Path(__file__).parent / "results_parametric_tet4.pck"
    log.info(f"Saving results to {results_file}")
    with open(results_file, "wb") as f:
        pickle.dump((X, Y, val_p, val_f), f)
