def WriteGmshGeoFileAndMakeTheMesh(filename,parameters):
    angle=parameters[0]
    lx3 = parameters[1]
    ly3 = parameters[2]
    l4 = parameters[3]

    f=open(filename+'.geo','w')
    f.write('// Parameters: acoustic cavity\n')
    f.write('lx1 = 7.0;\n')
    f.write('ly1 = 4.0;\n')
    f.write('lz1 = 2.8;\n')

    f.write('ly2 = 3.0;\n')
    f.write('lx2 = 3.5;\n')

    f.write('lx3 = %9.4g;\n' % lx3)
    f.write('ly3 = %9.4g;\n' % ly3)
    f.write('l4 = %9.4g;\n' %l4)
    f.write('lz4 = 2.0;\n')
    f.write('r4   = 0.5;\n')
    f.write('e4   = 0.2;\n')
    f.write('angle = %9.4g;\n' % angle)
    f.write('h = lx1/10;\n')

    f.write('lx5 = 6.0;\n')
    f.write('ly5 = 1.0;\n')
    f.write('lz5 = 1.5;\n')
    f.write('h5  = 0.5;\n')

    f.write('Mesh.CharacteristicLengthMax=h;\n')
    f.write('Mesh.ElementOrder = 1;\n')

    f.write('// Cavity: Corners\n')
    f.write('Point(1) = {0,     0  , 0, h};\n')
    f.write('Point(2) = {lx1,    0  , 0, h};\n')
    f.write('Point(3) = {lx1,    ly1 , 0, h};\n')
    f.write('Point(4) = {0 ,    ly1 , 0, h};\n')

    f.write('Point(5) = {0,     0  , lz1, h};\n')
    f.write('Point(6) = {lx1,    0  , lz1, h};\n')
    f.write('Point(7) = {lx1,    ly1 , lz1, h};\n')
    f.write('Point(8) = {0 ,    ly1 , lz1, h};\n')

    f.write('Point(9) = {0,     ly1+ly2  , 0, h};\n')
    f.write('Point(10) = {lx2,     ly1+ly2  , 0, h};\n')
    f.write('Point(11) = {lx2,     ly1  , 0, h};\n')
    f.write('Point(12) = {0,     ly1+ly2  , lz1, h};\n')
    f.write('Point(13) = {lx2,     ly1+ly2  , lz1, h};\n')
    f.write('Point(14) = {lx2,     ly1  , lz1, h};\n')

    f.write('// Cavity: lines\n')
    f.write('Line(1) = {1, 4};\n')
    f.write('Line(2) = {4, 9};\n')
    f.write('Line(3) = {9, 10};\n')
    f.write('Line(4) = {10, 11};\n')
    f.write('Line(5) = {11, 3};\n')
    f.write('Line(6) = {3, 2};\n')
    f.write('Line(7) = {2, 1};\n')
    f.write('Line(8) = {1, 5};\n')
    f.write('Line(9) = {4, 8};\n')
    f.write('Line(10) = {9, 12};\n')
    f.write('Line(11) = {10, 13};\n')
    f.write('Line(12) = {11, 14};\n')
    f.write('Line(13) = {3, 7};\n')
    f.write('Line(14) = {2, 6};\n')
    f.write('Line(15) = {5, 8};\n')
    f.write('Line(16) = {8, 12};\n')
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

    f.write('Point(23) = {lx3,     ly3+e4/2  , 0, h};\n')
    f.write('Point(24) = {lx3+l4,     ly3+e4/2  , 0, h};\n')
    f.write('Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h};\n')
    f.write('Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h};\n')
    f.write('Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h};\n')
    f.write('Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h};\n')
    f.write('Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , h};\n')
    f.write('Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, h};\n')

    f.write('// rotate wall\n')
    f.write('Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {\n')
    f.write('Point{22, 30, 18, 26, 20, 21, 28, 29, 19, 27, 17, 25, 16, 24, 15, 23};\n')
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

    f.write('Line(34) = {16, 24};\n')
    f.write('Line(35) = {18, 26};\n')
    f.write('Line(36) = {22, 30};\n')
    f.write('Line(37) = {21, 29};\n')
    f.write('Line(38) = {17, 25};\n')
    f.write('Line(39) = {15, 23};\n')

    f.write('// surfaces, cavity\n')
    f.write('Line Loop(40) = {1, 2, 3, 4, 5, 6, 7};\n')
    f.write('Line Loop(41) = {22, 34, -28, -39};\n')
    f.write('Plane Surface(42) = {40, 41};\n')
    f.write('Line Loop(43) = {4, 12, -18, -11};\n')
    f.write('Plane Surface(44) = {43};\n')
    f.write('Line Loop(45) = {11, -17, -10, 3};\n')
    f.write('Plane Surface(46) = {45};\n')
    f.write('Line Loop(47) = {10, -16, -9, 2};\n')
    f.write('Plane Surface(48) = {47};\n')
    f.write('Line Loop(49) = {9, -15, -8, 1};\n')
    f.write('Plane Surface(50) = {49};\n')
    f.write('Line Loop(51) = {8, -21, -14, 7};\n')
    f.write('Plane Surface(52) = {51};\n')
    f.write('Line Loop(53) = {14, -20, -13, 6};\n')
    f.write('Plane Surface(54) = {53};\n')
    f.write('Line Loop(55) = {13, -19, -12, 5};\n')
    f.write('Plane Surface(56) = {55};\n')
    f.write('Line Loop(57) = {20, 21, 15, 16, 17, 18, 19};\n')
    f.write('Plane Surface(58) = {57};\n')

    f.write('// surface, internal wall\n')
    f.write('Line Loop(59) = {34, 29, -35, -23};\n')
    f.write('Plane Surface(60) = {59};\n')
    f.write('Line Loop(61) = {30, -37, -24, 36};\n')
    f.write('Plane Surface(62) = {61};\n')
    f.write('Line Loop(63) = {31, -39, -25, 38};\n')
    f.write('Plane Surface(64) = {63};\n')
    f.write('Line Loop(65) = {35, 32, -36, -26};\n')
    f.write('Ruled Surface(66) = {65};\n')
    f.write('Line Loop(67) = {37, 33, -38, -27};\n')
    f.write('Ruled Surface(68) = {67};\n')
    f.write('Line Loop(69) = {22, 23, 26, 24, 27, 25};\n')
    f.write('Plane Surface(70) = {69};\n')
    f.write('Line Loop(71) = {29, 32, 30, 33, 31, 28};\n')
    f.write('Plane Surface(72) = {71};\n')

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

    f.write('// volume\n')
    f.write('Surface Loop(99) = {98, 88, 96, 90, 92, 94};\n')
    f.write('Volume(2) = {99};\n')
    f.write('Surface Loop(101) = {58, 54, 52, 50, 48, 46, 44, 42, 56, 64, 72, 60, 66, 62, 68, 70};\n')
    f.write('Volume(1) = {99, 101};\n')

    f.write('Physical Volume(1) = {1};\n')
    f.write('Physical Volume(2) = {2};\n')
    f.close()
    import os
    os.system('gmsh -3 %s.geo' % filename)

    return

def ComputeFRF(parameters):
    import string
    import time
    import scipy
    import scipy.sparse
    import scipy.sparse.linalg
    import pylab as pl
    import pickle
    import sys
    sys.path.append('../../../librairies')
    import silex_lib_xfem_acou_tet4
    import silex_lib_dkt_fortran
    import silex_lib_gmsh
    import mumps
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    nproc=comm.Get_size()
    rank = comm.Get_rank()
    class comm_mumps_one_proc:
        rank = 0
        def py2f(self):
            return 0
    mycomm=comm_mumps_one_proc()
    angle = parameters[0]
    lx3   = parameters[1]
    ly3   = parameters[2]
    l4    = parameters[3]

    mesh_file='cavity3'
    results_file='results/cavity3_damping'
    import ComputeCavity3
    ComputeCavity3.WriteGmshGeoFileAndMakeTheMesh('cavity3',[angle,1.5,2.0,3.0])

    freq_ini     = 20.0
    freq_end     = 100.0
    nb_freq_step_per_proc=200

    nb_freq_step = nb_freq_step_per_proc*nproc
    deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

    # air
    celerity=340.0
    rho=1.2
    fluid_damping=(1.0+0.02j)
    #fluid_damping=1.0

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'.msh',3)
    fluid_elements1,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,1)
    fluid_elements2,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',4,2)
    fluid_elements=scipy.vstack([fluid_elements1,fluid_elements2])
    #fluid_elements_S2,IdNodesS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'.msh',2,2)

    fluid_nnodes   = fluid_nodes.shape[0]
    fluid_nelem    = fluid_elements.shape[0]
    fluid_ndof     = fluid_nnodes

    print ("Number of fluid nodes:",fluid_nnodes)
    print ("Number of fluid elements:",fluid_nelem)

    #silex_lib_gmsh.WriteResults(results_file+'Mesh',fluid_nodes,fluid_elements,4)
    #silex_lib_gmsh.WriteResults(results_file+'Mesh_surface',fluid_nodes,fluid_elements_S2,2)

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.clock()

    IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

    KFF=scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    MFF=scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    SolvedDofF=list(range(fluid_ndof))

    ##################################################################
    # Construct the whole system
    ##################################################################

    K=KFF[SolvedDofF,:][:,SolvedDofF]
    M=MFF[SolvedDofF,:][:,SolvedDofF]

    # To impose the load on the fluid:
    # fluid node number 1 is at (0,-ly/2,0)
    # node number 1 is at (0,-ly/2,0)
    #F = csc_matrix( ([1],([0],[0])), shape=(len(SolvedDofS)+len(SolvedDofF),1) )

    P=scipy.zeros((fluid_ndof))
    P[13-1]=1.0
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

            F=scipy.array(omega**2*P , dtype='c16')
            #F=scipy.array(P , dtype='c16')

            #sol=scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype=complex) , scipy.array(F.todense() , dtype=complex) )
            sol = mumps.spsolve( scipy.sparse.csc_matrix(fluid_damping*K-(omega**2)*M,dtype='c16') , F , comm=mycomm )
            #sol = mumps.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float') , F , comm=mycomm )
            

            #press = scipy.zeros((fluid_ndof),dtype=float)
            press = scipy.zeros((fluid_ndof),dtype=complex)
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements2,fluid_nodes,press))
            #frf[i]=silex_lib_xfem_acou_tet4.computequadratiquepressure(fluid_elements,fluid_nodes,press)
            i=i+1

            if rank==0:
                press_save.append(press.real)

        frfsave=[frequencies,frf]
        comm.send(frfsave, dest=0, tag=11)

        if rank==0:
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])

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

    return Allfrfsave
