def WriteGmshGeoFileAndMakeTheMesh(filename,parameters):
    angle=parameters[0]

    f=open(filename+'.geo','w')
    f.write('// Parameters: acoustic cavity\n')
    f.write('lx1 = 7.0;\n')
    f.write('ly1 = 4.0;\n')
    f.write('lz1 = 2.5;\n')
    f.write('ly2 = 3.0;\n')
    f.write('lx2 = 3.5;\n')
    f.write('lx3 = 1.5;\n')
    f.write('ly3 = 2.0;\n')
    f.write('l4 = 3.0;\n')
    f.write('lz4 = 2.4;\n')
    f.write('r4   = 0.5;\n')
    f.write('e4   = 0.1;\n')
    f.write('deg=      %9.4g0;\n' % angle)
    f.write('angle = deg*Pi/180;\n')
    f.write('// porous on the moving wall\n')
    f.write('ee = 0.1;\n')
    f.write('// porous on the fixed wall\n')
    f.write('ee6 = 0.2;\n')
    f.write('dy6 = 0.1;\n')
    f.write('dz6 = 0.1;\n')
    f.write('ly6 = 6.8;\n')
    f.write('lz6 = 2.3;\n')
    f.write('\n')
    f.write('// size of elements\n')
    f.write('h = lx1/10;\n')
    f.write('h2 = ee/1;\n')
    f.write('h6 = ee6/1;\n')
    f.write('\n')
    f.write('lx5 = 6.0;\n')
    f.write('ly5 = 1.0;\n')
    f.write('lz5 = 1.5;\n')
    f.write('h5  = 1.0;\n')
    f.write('Mesh.CharacteristicLengthMax=2*h;\n')
    f.write('Mesh.ElementOrder = 1;\n')
    f.write('// Cavity: Corners\n')
    f.write('Point(1) = {0,     0  , 0, h};\n')
    f.write('Point(2) = {lx1,    0  , 0, h};\n')
    f.write('Point(3) = {lx1,    ly1 , 0, h};\n')
    f.write('Point(5) = {0,     0  , lz1, h};\n')
    f.write('Point(6) = {lx1,    0  , lz1, h};\n')
    f.write('Point(7) = {lx1,    ly1 , lz1, h};\n')
    f.write('Point(9) = {0,     ly1+ly2  , 0, h};\n')
    f.write('Point(10) = {lx2,     ly1+ly2  , 0, h};\n')
    f.write('Point(11) = {lx2,     ly1  , 0, h};\n')
    f.write('Point(12) = {0,     ly1+ly2  , lz1, h};\n')
    f.write('Point(13) = {lx2,     ly1+ly2  , lz1, h};\n')
    f.write('Point(14) = {lx2,     ly1  , lz1, h};\n')
    f.write('// Cavity: lines\n')
    f.write('Line(2) = {1, 9};\n')
    f.write('Line(3) = {9, 10};\n')
    f.write('Line(4) = {10, 11};\n')
    f.write('Line(5) = {11, 3};\n')
    f.write('Line(6) = {3, 2};\n')
    f.write('Line(7) = {2, 1};\n')
    f.write('Line(8) = {1, 5};\n')
    f.write('Line(10) = {9, 12};\n')
    f.write('Line(11) = {10, 13};\n')
    f.write('Line(12) = {11, 14};\n')
    f.write('Line(13) = {3, 7};\n')
    f.write('Line(14) = {2, 6};\n')
    f.write('Line(16) = {5, 12};\n')
    f.write('Line(17) = {12, 13};\n')
    f.write('Line(18) = {13, 14};\n')
    f.write('Line(19) = {14, 7};\n')
    f.write('Line(20) = {7, 6};\n')
    f.write('Line(21) = {6, 5};\n')
    f.write('// wall: corners\n')
    f.write('Point(15) = {lx3,     ly3-e4/2  , 0, h};\n')
    f.write('Point(16) = {lx3+l4,     ly3-e4/2  , 0, h};\n')
    f.write('Point(17) = {lx3,     ly3-e4/2  , lz4-r4 , h};\n')
    f.write('Point(18) = {lx3+l4,     ly3-e4/2  , lz4-r4, h};\n')
    f.write('Point(19) = {lx3+r4,     ly3-e4/2  , lz4-r4 , h};\n')
    f.write('Point(20) = {lx3+l4-r4,     ly3-e4/2  , lz4-r4, h};\n')
    f.write('Point(21) = {lx3+r4,     ly3-e4/2  , lz4 , h};\n')
    f.write('Point(22) = {lx3+l4-r4,     ly3-e4/2  , lz4, h};\n')
    f.write('Point(23) = {lx3,     ly3+e4/2  , 0, h2};\n')
    f.write('Point(24) = {lx3+l4,     ly3+e4/2  , 0, h2};\n')
    f.write('Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h2};\n')
    f.write('Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h2};\n')
    f.write('Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h2};\n')
    f.write('Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h2};\n')
    f.write('Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , h2};\n')
    f.write('Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, h2};\n')
    f.write('Point(123) = {lx3,     ly3+e4/2+ee  , 0, h2};\n')
    f.write('Point(124) = {lx3+l4,     ly3+e4/2+ee  , 0, h2};\n')
    f.write('Point(125) = {lx3,     ly3+e4/2+ee  , lz4-r4 , h2};\n')
    f.write('Point(126) = {lx3+l4,     ly3+e4/2+ee  , lz4-r4, h2};\n')
    f.write('Point(127) = {lx3+r4,     ly3+e4/2+ee  , lz4-r4 , h2};\n')
    f.write('Point(128) = {lx3+l4-r4,     ly3+e4/2+ee  , lz4-r4, h2};\n')
    f.write('Point(129) = {lx3+r4,     ly3+e4/2+ee  , lz4 , h2};\n')
    f.write('Point(130) = {lx3+l4-r4,     ly3+e4/2+ee  , lz4, h2};\n')
    f.write('// rotate wall\n')
    f.write('Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {\n')
    f.write('Point{22, 30, 18, 26, 20, 21, 28, 29, 19, 27, 17, 25, 16, 24, 15, 23,123,124,125,126,127,128,129,130};\n')
    f.write('}\n')
    f.write('// wall : lines\n')
    f.write('Line(22) = {15, 16};\n')
    f.write('Line(23) = {16, 18};\n')
    f.write('Line(24) = {22, 21};\n')
    f.write('Line(25) = {17, 15};\n')
    f.write('Circle(26) = {18, 20, 22};\n')
    f.write('Circle(27) = {21, 19, 17};\n')
    f.write('Line(28) = {23, 24};\n')
    f.write('Line(29) = {24, 26};\n')
    f.write('Line(30) = {30, 29};\n')
    f.write('Line(31) = {25, 23};\n')
    f.write('Circle(32) = {26, 28, 30};\n')
    f.write('Circle(33) = {29, 27, 25};\n')
    f.write('Line(128) = {123, 124};\n')
    f.write('Line(129) = {124, 126};\n')
    f.write('Line(130) = {130, 129};\n')
    f.write('Line(131) = {125, 123};\n')
    f.write('Circle(132) = {126, 128, 130};\n')
    f.write('Circle(133) = {129, 127, 125};\n')
    f.write('Line(34) = {16, 24};\n')
    f.write('Line(35) = {18, 26};\n')
    f.write('Line(36) = {22, 30};\n')
    f.write('Line(37) = {21, 29};\n')
    f.write('Line(38) = {17, 25};\n')
    f.write('Line(39) = {15, 23};\n')
    f.write('Line(134) = {123, 23};\n')
    f.write('Line(135) = {125, 25};\n')
    f.write('Line(136) = {129, 29};\n')
    f.write('Line(137) = {130, 30};\n')
    f.write('Line(138) = {126, 26};\n')
    f.write('Line(139) = {124, 24};\n')
    f.write('\n')
    f.write('// Porous, surfaces and volume, on the moving wall\n')
    f.write('Line Loop(140) = {32, 30, 33, 31, 28, 29};\n')
    f.write('Plane Surface(141) = {140};\n')
    f.write('Line Loop(142) = {138, 32, -137, -132};\n')
    f.write('Ruled Surface(143) = {142};\n')
    f.write('Line Loop(144) = {130, 136, -30, -137};\n')
    f.write('Plane Surface(145) = {144};\n')
    f.write('Line Loop(146) = {133, 135, -33, -136};\n')
    f.write('Ruled Surface(147) = {146};\n')
    f.write('Line Loop(148) = {134, -31, -135, 131};\n')
    f.write('Plane Surface(149) = {148};\n')
    f.write('Line Loop(150) = {134, 28, -139, -128};\n')
    f.write('Plane Surface(151) = {150};\n')
    f.write('Line Loop(152) = {139, 29, -138, -129};\n')
    f.write('Plane Surface(153) = {152};\n')
    f.write('Line Loop(154) = {128, 129, 132, 130, 133, 131};\n')
    f.write('Plane Surface(155) = {154};\n')
    f.write('Surface Loop(156) = {155, 151, 149, 141, 143, 153, 145, 147};\n')
    f.write('Volume(157) = {156};\n')
    f.write('\n')
    f.write('// porous on the fixed wall\n')
    f.write('Point(300) = {-ee6,  dy6  , dz6, h6};\n')
    f.write('Point(301) = {-ee6,  dy6  , dz6+lz6, h6};\n')
    f.write('Point(302) = {-ee6,  dy6+ly6  , dz6+lz6, h6};\n')
    f.write('Point(303) = {-ee6,  dy6+ly6  , dz6, h6};\n')
    f.write('Point(304) = {0,  dy6  , dz6, h6};\n')
    f.write('Point(305) = {0,  dy6  , dz6+lz6, h6};\n')
    f.write('Point(306) = {0,  dy6+ly6  , dz6+lz6, h6};\n')
    f.write('Point(307) = {0,  dy6+ly6  , dz6, h6};\n')
    f.write('\n')
    f.write('Line(188) = {304, 305};\n')
    f.write('Line(189) = {305, 306};\n')
    f.write('Line(190) = {306, 307};\n')
    f.write('Line(191) = {307, 304};\n')
    f.write('Line(192) = {300, 301};\n')
    f.write('Line(193) = {301, 302};\n')
    f.write('Line(194) = {302, 303};\n')
    f.write('Line(195) = {303, 300};\n')
    f.write('Line(196) = {300, 304};\n')
    f.write('Line(197) = {301, 305};\n')
    f.write('Line(198) = {302, 306};\n')
    f.write('Line(199) = {303, 307};\n')
    f.write('Line Loop(200) = {188, 189, 190, 191};\n')
    f.write('Plane Surface(201) = {200};\n')
    f.write('Line Loop(202) = {197, 189, -198, -193};\n')
    f.write('Plane Surface(203) = {202};\n')
    f.write('Line Loop(204) = {196, 188, -197, -192};\n')
    f.write('Plane Surface(205) = {204};\n')
    f.write('Line Loop(206) = {195, 196, -191, -199};\n')
    f.write('Plane Surface(207) = {206};\n')
    f.write('Line Loop(208) = {199, -190, -198, 194};\n')
    f.write('Plane Surface(209) = {208};\n')
    f.write('Line Loop(210) = {195, 192, 193, 194};\n')
    f.write('Plane Surface(211) = {210};\n')
    f.write('\n')
    f.write('Surface Loop(212) = {207, 211, 205, 203, 209, 201};\n')
    f.write('Volume(213) = {212};\n')
    f.write('\n')
    f.write('\n')
    f.write('// surfaces, cavity\n')
    f.write('Line Loop(170) = {6, 7, 2, 3, 4, 5};\n')
    f.write('Line Loop(171) = {22, 34, -139, -128, 134, -39};\n')
    f.write('Plane Surface(172) = {170, 171};\n')
    f.write('Line Loop(173) = {6, 14, -20, -13};\n')
    f.write('Plane Surface(174) = {173};\n')
    f.write('Line Loop(175) = {5, 13, -19, -12};\n')
    f.write('Plane Surface(176) = {175};\n')
    f.write('Line Loop(177) = {4, 12, -18, -11};\n')
    f.write('Plane Surface(178) = {177};\n')
    f.write('Line Loop(179) = {3, 11, -17, -10};\n')
    f.write('Plane Surface(180) = {179};\n')
    f.write('//Line Loop(181) = {2, 10, -16, -8};\n')
    f.write('//Plane Surface(182) = {181};\n')
    f.write('//Line Loop(181) = {8, 16, -10, -2};\n')
    f.write('//Plane Surface(182) = {200, 211};\n')
    f.write('\n')
    f.write('Line Loop(181) = {16, -10, -2, 8};\n')
    f.write('Plane Surface(182) = {200, 181};\n')
    f.write('\n')
    f.write('Line Loop(183) = {7, 8, -21, -14};\n')
    f.write('Plane Surface(184) = {183};\n')
    f.write('Line Loop(185) = {21, 16, 17, 18, 19, 20};\n')
    f.write('Plane Surface(186) = {185};\n')
    f.write('// surface, internal wall\n')
    f.write('Line Loop(158) = {39, -31, -38, 25};\n')
    f.write('Plane Surface(159) = {158};\n')
    f.write('Line Loop(160) = {38, -33, -37, 27};\n')
    f.write('Ruled Surface(161) = {160};\n')
    f.write('Line Loop(162) = {37, -30, -36, 24};\n')
    f.write('Plane Surface(163) = {162};\n')
    f.write('Line Loop(164) = {36, -32, -35, 26};\n')
    f.write('Ruled Surface(165) = {164};\n')
    f.write('Line Loop(166) = {35, -29, -34, 23};\n')
    f.write('Plane Surface(167) = {166};\n')
    f.write('Line Loop(168) = {22, 23, 26, 24, 27, 25};\n')
    f.write('Plane Surface(169) = {168};\n')
    f.write('\n')
    f.write('// control volume\n')
    f.write('Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, h};\n')
    f.write('Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, h};\n')
    f.write('Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, h};\n')
    f.write('Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, h};\n')
    f.write('Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, h};\n')
    f.write('Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, h};\n')
    f.write('Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, h};\n')
    f.write('Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, h};\n')
    f.write('Line(75) = {33, 37};\n')
    f.write('Line(76) = {37, 35};\n')
    f.write('Line(77) = {35, 31};\n')
    f.write('Line(78) = {31, 33};\n')
    f.write('Line(79) = {33, 34};\n')
    f.write('Line(80) = {37, 38};\n')
    f.write('Line(81) = {38, 36};\n')
    f.write('Line(82) = {36, 32};\n')
    f.write('Line(83) = {32, 34};\n')
    f.write('Line(84) = {34, 38};\n')
    f.write('Line(85) = {35, 36};\n')
    f.write('Line(86) = {31, 32};\n')
    f.write('Line Loop(87) = {75, 80, -84, -79};\n')
    f.write('Plane Surface(88) = {87};\n')
    f.write('Line Loop(89) = {76, 85, -81, -80};\n')
    f.write('Plane Surface(90) = {89};\n')
    f.write('Line Loop(91) = {77, 86, -82, -85};\n')
    f.write('Plane Surface(92) = {91};\n')
    f.write('Line Loop(93) = {78, 79, -83, -86};\n')
    f.write('Plane Surface(94) = {93};\n')
    f.write('Line Loop(95) = {75, 76, 77, 78};\n')
    f.write('Plane Surface(96) = {95};\n')
    f.write('Line Loop(97) = {84, 81, 82, 83};\n')
    f.write('Plane Surface(98) = {97};\n')
    f.write('\n')
    f.write('// controlled volume\n')
    f.write('Surface Loop(99) = {98, 88, 96, 90, 92, 94};\n')
    f.write('Volume(2) = {99};\n')
    f.write('\n')
    f.write('// volume: air cavity + controlled volume\n')
    f.write('Surface Loop(187) = {155, 153, 143, 145, 147, 149, 163, 161, 159, 172, 174, 184, 182, 180, 178, 176, 186, 167, 165, 169,201};\n')
    f.write('Volume(1) = {99, 187};\n')
    f.write('//Surface Loop(101) = {58, 54, 52, 50, 48, 46, 44, 42, 56, 64, 72, 60, 66, 62, 68, 70};\n')
    f.write('//Volume(1) = {99, 101};\n')
    f.write('\n')
    f.write('Physical Volume(1) = {1,2}; // air cavity + controlled air volume\n')
    f.write('Physical Volume(5) = {2}; // controlled air volume (small)\n')
    f.write('// Porous\n')
    f.write('Physical Volume(2) = {157,213}; // porous volume\n')
    f.write('Physical Surface(4) = {153, 143, 145, 147, 149, 151, 141,211,203,209,207,205}; // porous external surface\n')
    f.write('Physical Surface(3) = {155,201}; // [porous]-[air cavity] : interface surface\n')
    f.close()
    import os
    os.system('gmsh -3 %s.geo' % filename)


def ComputeFRF(parameters):
    import string
    import time
    import scipy
    import scipy.sparse
    import scipy.sparse.linalg

    #import os
    import pylab as pl
    import pickle

    import sys
    sys.path.append('../../../librairies')

    import silex_lib_xfem_acou_tet4
    import silex_lib_dkt_fortran
    import silex_lib_gmsh

    #import silex_lib_extra_fortran as silex_lib_extra

    import silex_lib_porous_tet4_fortran

    import mumps

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc= comm.Get_size()
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
    # To run it in sequentiel frequency per frequency with openblas in parallel:
    # export OPENBLAS_NUM_THREADS=10
    # python3.4 Main_toto.py
    #
    ##############################################################
    ##############################################################
    #              S T A R T   M A I N   P R O B L E M
    ##############################################################
    ##############################################################

    # parallepipedic cavity with plane structure
    mesh_file='geom/cavity5_with_porous'
    results_file='results/cavity5_with_porous'
    angle = parameters[0]
    import ComputePorousCavity5
    ComputePorousCavity5.WriteGmshGeoFileAndMakeTheMesh(mesh_file,[angle])


    freq_ini     = 1.0
    freq_end     = 1000.0
    nb_freq_step_per_proc=50

    nb_freq_step = nb_freq_step_per_proc*nproc
    deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

    # air
    celerity=343.0 # ok
    rho=1.21 # ok

    # porous properties
    E_sol = 286000.0*(3*1144.0+2*286.0)/(1144.0+286.0) #ok =800800.0
    nu_sol = 1144.0/(2*(1144.0+286.0)) #ok =0.4
    #ro_sol = 30.0/(1.0-0.9) #ok = 300.0
    ro_sol = 30.0

    to_por = 7.8 # ok
    po_por = 0.9 # ok
    #po_por = 0.1 # almost only solid
    sg_por = 25000.0 # ok
    lambda_por = 226e-6 # ok
    lambda_prime_por = 226e-6 # ok
    #lambda_por = 226e3 # ok
    #lambda_prime_por = 226e3 # ok

    ro_fl = 1.21 # ok
    visco_fl = 0.0000184 # ok
    pdtl_fl = 0.71 # ok
    gamma_fl = 1.4 # ok
    p0_fl = 101325.0
    ce_fl = celerity 


    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
    fluid_elements1,IdNodes1 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1) # air, cavity + controlled volume
    fluid_elements5,IdNodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,5) # air, ONLY controlled volume
    fluid_elements_S3,IdNodesS3 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,3)# air-porous interface surface
    #fluid_elements_S5,IdNodesS5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,5)
    #fluid_elements_S6,IdNodesS6 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,6)

    fluid_elements2,IdNodes2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2) # porous, volume
    fluid_elements_S4,IdNodesS4 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,4) # porous, external surface

    fluid_nnodes   = fluid_nodes.shape[0]
    fluid_nelem1   = fluid_elements1.shape[0]
    fluid_nelem2   = fluid_elements2.shape[0]
    fluid_nelem3   = fluid_elements_S3.shape[0]
    fluid_nelem4   = fluid_elements_S4.shape[0]
    fluid_nelem5   = fluid_elements5.shape[0]

    fluid_nnodes1 = IdNodes1.shape[0]
    fluid_nnodes2 = IdNodes2.shape[0]
    fluid_nnodes3 = IdNodesS3.shape[0]
    fluid_nnodes4 = IdNodesS4.shape[0]
    fluid_nnodes5 = IdNodes5.shape[0]

    fluid_ndof1     = fluid_nnodes1
    fluid_ndof2     = fluid_nnodes2*6

    print ("Number of nodes:",fluid_nnodes)
    print ("Number of elements in air:",fluid_nelem1)
    print ("Number of elements in porous:",fluid_nelem2)

    print ("Number of nodes in air:",fluid_nnodes1)
    print ("Number of nodes in porous:",fluid_nnodes2)
    print ("Number of nodes at interface:",fluid_nnodes3)


    ##############################################################
    # renumbering
    ##############################################################

    # renumbering air
    old = scipy.unique(fluid_elements1)
    new = list(range(1,len(scipy.unique(fluid_elements1))+1))
    new_nodes=fluid_nodes[scipy.unique(fluid_elements1)-1,:]

    dico1 = dict(zip(old,new))
    new_elements=scipy.zeros((fluid_nelem1,4),dtype=int)
    for e in range(fluid_nelem1):
        for i in range(4):
            new_elements[e,i]=dico1[fluid_elements1[e][i]]

    fluid_elements1 = new_elements
    fluid_nodes1    = new_nodes

    new_elements=scipy.zeros((fluid_nelem5,4),dtype=int)
    for e in range(fluid_nelem5):
        for i in range(4):
            new_elements[e,i]=dico1[fluid_elements5[e][i]]

    fluid_elements5 = new_elements

    # renumbering porous
    old = scipy.unique(fluid_elements2)
    new = list(range(1,len(scipy.unique(fluid_elements2))+1))
    new_nodes=fluid_nodes[scipy.unique(fluid_elements2)-1,:]

    dico2 = dict(zip(old,new))
    new_elements=scipy.zeros((fluid_nelem2,4),dtype=int)
    for e in range(fluid_nelem2):
        for i in range(4):
            new_elements[e,i]=dico2[fluid_elements2[e][i]]

    fluid_elements2 = new_elements
    fluid_nodes2    = new_nodes

    # renumbering porous external surfaces
    #print(IdNodesS4)
    for i in range(len(IdNodesS4)):
        IdNodesS4[i]=dico2[IdNodesS4[i]]
    ##for i in range(len(IdNodesS5)):
    ##    IdNodesS5[i]=dico2[IdNodesS5[i]]
    ##for i in range(len(IdNodesS6)):
    ##    IdNodesS6[i]=dico2[IdNodesS6[i]]

    #print(IdNodesS4)
    #stop

    IdNodesFixed_porous_us_x=IdNodesS4
    ##IdNodesFixed_porous_us_y=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS6]))
    ##IdNodesFixed_porous_us_z=scipy.unique(scipy.hstack([IdNodesS4,IdNodesS5]))
    IdNodesFixed_porous_us_y=IdNodesS4
    IdNodesFixed_porous_us_z=IdNodesS4
    IdNodesFixed_porous_uf_x=IdNodesS4
    IdNodesFixed_porous_uf_y=IdNodesS4
    IdNodesFixed_porous_uf_z=IdNodesS4

    Fixed_Dofs_porous = scipy.hstack([(IdNodesFixed_porous_us_x-1)*6,
                                      (IdNodesFixed_porous_us_y-1)*6+1,
                                      (IdNodesFixed_porous_us_z-1)*6+2,
                                      (IdNodesFixed_porous_uf_x-1)*6+3,
                                      (IdNodesFixed_porous_uf_y-1)*6+4,
                                      (IdNodesFixed_porous_uf_z-1)*6+5
                                      ])

    # get connectivity at interface
    IdNodesS3_for_1=scipy.zeros(fluid_nnodes3,dtype=int)
    IdNodesS3_for_2=scipy.zeros(fluid_nnodes3,dtype=int)
    for i in range(fluid_nnodes3):
        IdNodesS3_for_1[i]=dico1[IdNodesS3[i]]
        IdNodesS3_for_2[i]=dico2[IdNodesS3[i]]

    InterfaceConnectivity=scipy.zeros((fluid_nelem3,6),dtype=int)
    for e in range(fluid_nelem3):
        for i in range(3):
            InterfaceConnectivity[e,i]   = dico1[fluid_elements_S3[e,i]]
            InterfaceConnectivity[e,i+3] = dico2[fluid_elements_S3[e,i]]


    silex_lib_gmsh.WriteResults(results_file+'Mesh1',fluid_nodes1,fluid_elements1,4)
    silex_lib_gmsh.WriteResults(results_file+'Mesh2',fluid_nodes2,fluid_elements2,4)
    silex_lib_gmsh.WriteResults(results_file+'Mesh_surface3',fluid_nodes,fluid_elements_S3,2)
    silex_lib_gmsh.WriteResults(results_file+'Mesh_surface4',fluid_nodes,fluid_elements_S4,2)
    silex_lib_gmsh.WriteResults(results_file+'Mesh5',fluid_nodes1,fluid_elements5,4)

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.clock()

    IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements1,fluid_nodes1,celerity,rho)

    KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )
    MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof1,fluid_ndof1) )

    SolvedDofF=list(range(fluid_ndof1))

    ##############################################################
    # Compute Porous Matrices
    ##############################################################
    porous_material_prop=[E_sol,nu_sol,ro_sol,to_por,po_por,sg_por,lambda_por,lambda_prime_por,ro_fl,visco_fl,pdtl_fl,gamma_fl,p0_fl,ce_fl]
    #omega=2.0*scipy.pi*100.0
    #IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)

    #KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
    #MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )

    SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),Fixed_Dofs_porous)

    #SolvedDofP=scipy.setdiff1d(range(fluid_ndof2),range(fluid_ndof2))

    ##############################################################
    # Compute Coupling Porous-air Matrices
    ##############################################################

    #print(silex_lib_porous_tet4_fortran.computecouplingporousair.__doc__)

    IIpf,JJpf,Vpf=silex_lib_porous_tet4_fortran.computecouplingporousair(fluid_nodes1,InterfaceConnectivity,po_por)
    CPF=scipy.sparse.csc_matrix( (Vpf,(IIpf,JJpf)), shape=(fluid_ndof2,fluid_ndof1) )

    SolvedDof = scipy.hstack([SolvedDofF,SolvedDofP+fluid_ndof1])

    ##################################################################
    # Construct the whole system
    ##################################################################

    #K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],-CPF[SolvedDofP,:][:,SolvedDofF].T],
    #                                 [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
    #M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
    #                                 [CPF[SolvedDofP,:][:,SolvedDofF],MPP[SolvedDofP,:][:,SolvedDofP]] ] )

    # To impose the load on the fluid:
    # fluid node number 1
    #F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )
    #ff = 25e-8/(8e-3)! for surface load:/(4*nz*ny)
    UF = scipy.zeros(fluid_ndof1+fluid_ndof2,dtype=float)
    UF[11-1]=3.1250E-05

    #P=scipy.zeros((fluid_ndof))
    #P[13-1]=1.0
    #print(silex_lib_xfem_acou_tet4.forceonsurface.__doc__)
    #P = silex_lib_xfem_acou_tet4.forceonsurface(fluid_nodes,fluid_elements_S2,1.0)


    ##############################################################
    # FRF computation
    ##############################################################



    Flag_frf_analysis=1
    frequencies=[]
    frf=[]

    if (Flag_frf_analysis==1):
        print ("Proc. ",rank," / time at the beginning of the FRF:",time.ctime())

        press_save=[]
        disp_save=[]

        for i in range(nb_freq_step_per_proc):
        #for freq in scipy.linspace(freq_ini,freq_end,nb_freq_step):

            freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega=2*scipy.pi*freq

            print ("proc number",rank,"frequency=",freq)

            IIp,JJp,Vppk,Vppm=silex_lib_porous_tet4_fortran.stiffnessmassmatrix(fluid_nodes2,fluid_elements2,porous_material_prop,omega)
            KPP=scipy.sparse.csc_matrix( (Vppk,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )
            MPP=scipy.sparse.csc_matrix( (Vppm,(IIp,JJp)), shape=(fluid_ndof2,fluid_ndof2) )

            K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],-CPF[SolvedDofP,:][:,SolvedDofF].T],
                                             [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
            M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
                                             [CPF[SolvedDofP,:][:,SolvedDofF],MPP[SolvedDofP,:][:,SolvedDofP]] ] )
    ##        K=scipy.sparse.construct.bmat( [ [KFF[SolvedDofF,:][:,SolvedDofF],None],
    ##                                         [None,KPP[SolvedDofP,:][:,SolvedDofP]] ] )
    ##        M=scipy.sparse.construct.bmat( [ [MFF[SolvedDofF,:][:,SolvedDofF],None],
    ##                                         [None,MPP[SolvedDofP,:][:,SolvedDofP]] ] )

            F=scipy.array(omega**2*UF[SolvedDof] , dtype='c16')
            #F=scipy.array(P , dtype='c16')

            #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M,dtype=complex) , scipy.array(F , dtype=complex) )
            sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
            #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )
            

            #press = scipy.zeros((fluid_ndof),dtype=float)
            press1 = scipy.zeros((fluid_ndof1),dtype=complex)
            press1[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes1,press1))
            #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
            #i=i+1
            

            if rank==0:
                press_save.append(press1.real)

        frfsave=[frequencies,frf]
        comm.send(frfsave, dest=0, tag=11)

        if rank==0:
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes1,fluid_elements1,4,[[press_save,'nodal',1,'pressure']])

        print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())

        # save the FRF problem
        Allfrequencies=scipy.zeros(nb_freq_step)
        Allfrf=scipy.zeros(nb_freq_step)
        k=0
        if rank==0:
            for i in range(nproc):
                data = comm.recv(source=i, tag=11)
                for j in range(len(data[0])):
                    Allfrequencies[k]=data[0][j]
                    Allfrf[k]=data[1][j]
                    k=k+1

            Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
            Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            f.close()


        
