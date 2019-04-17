###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# Libraries
###########################################################
import getopt
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io

import pickle

import sys
import os
from shutil import copyfile
sys.path.append('../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_gmsh
#import silex_lib_dkt_fortran as silex_lib_dkt

#import silex_lib_tet4_fortran_test as silex_lib_tet4

#import silex_lib_porous_tet4_fortran

# manage mumps librairie
import mumps
from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank = comm.Get_rank()

from mpi4py import MPI


def mpiInfo():
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    return nproc, rank, comm


class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm = comm_mumps_one_proc()

# distribution of frequencies per processor
def computeFreqPerProc(nbStep, nbProc, freqInit, freqEnd):
    # compute integer number of freq per proc and remaining steps
    nbFreqProc = nbStep // nbProc
    nbFreqProcRemain = nbStep % nbProc
    # compute frequencies steps
    varCase = 1
    if nbFreqProcRemain == 0:
        varCase = 0
    listFreq = scipy.zeros((nbFreqProc+varCase, nbProc))
    listAllFreq = scipy.linspace(freqInit, freqEnd, nbStep)
    # print(scipy.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1

    # print(listFreq)
    return listFreq

#load structure mesh fo values of parameters
def buildStructMesh(fileOrig,destFile,paraVal):
    #copy original file to the used one
    copyfile(fileOrig+'.geo',destFile+'.geo')
    #change value of parameters in the new file
    for key,value in enumerate(paraVal):
        oldText="<val##"+str(key)+">"
        newText='%g'%value
        #print(oldText)
        #print(newText) 
        cmdSed="sed -i 's/"+oldText+"/"+newText+"/g' "+destFile+'.geo'
        #print(cmdSed)
        os.system(cmdSed)
        
    #run gmsh to build the mesh
    #os.system('gmsh -3 -format msh2 '+destFile+'.geo')


###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7-v4.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################


def RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm, paraVal,gradValRequire=[],saveResults=1):#, caseDefine):

    print("##################################################")
    print("##################################################")
    print("##################################################")
    print("##    Start SILEX vibro-acoustics computation   ##")
    if len(gradValRequire)>0:
        print("##          (with gradients computation)        ##")
    print("##################################################")
    print("##################################################")

    # load 3D geometry
    orig_mesh_file = 'geom/cavity_acou3D_struc_3D_v3_para'
    mesh_file = 'geom/cavity_acou3D_struc_3D_v3'
    results_file_ini = 'results/cavity_acou3D_struc_3D_v3'

    listFreqPerProc = computeFreqPerProc(nbStep, nbProc, freqMin, freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity = 340.0
    rho = 1.2
    fluid_damping = (1+0.01j)

    nproc = comm.Get_size()
    rank = comm.Get_rank()

    flag_write_gmsh_results = saveResults

    flag_edge_enrichment = 0

    # number of parameters
    nbPara = len(paraVal)
    # prepare save file
    file_extension = "{:.{}E}".format(paraVal[0], 2)
    if (nbPara > 1) and (rank==0):
        for i in range(1, nbPara):
            file_extension = file_extension+'_'+"{:.{}E}".format(paraVal[i], 2)

    results_file = results_file_ini+'_'+file_extension
    print(results_file)

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_air.msh', 3)
    fluid_elements1, IdNodes1 = silex_lib_gmsh.ReadGmshElements(
        mesh_file+'_air.msh', 4, 1)  # air, cavity + controlled volume
    fluid_elements5, IdNodes5 = silex_lib_gmsh.ReadGmshElements(
        mesh_file+'_air.msh', 4, 5)  # air, ONLY controlled volume

    fluid_nnodes = fluid_nodes.shape[0]
    fluid_nelem1 = fluid_elements1.shape[0]
    fluid_nelem5 = fluid_elements5.shape[0]

    fluid_nnodes1 = IdNodes1.shape[0]
    fluid_nnodes5 = IdNodes5.shape[0]

    fluid_ndof = fluid_nnodes

    if rank == 0:
        print("Number of nodes:", fluid_nnodes)
        print("Number of elements in air:", fluid_nelem1)
        print("Number of nodes in air:", fluid_nnodes1)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults(
            results_file+'_air_cavity_Mesh1', fluid_nodes, fluid_elements1, 4)
        silex_lib_gmsh.WriteResults(
            results_file+'_air_controlled_volume_Mesh5', fluid_nodes, fluid_elements5, 4)

    ##############################################################
    # Load structure mesh
    ##############################################################
    #change parameters values and build mesh
    buildStructMesh(orig_mesh_file+'_struc',mesh_file+'_struc',paraVal)

    struc_nodes = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh', 3)
    struc_elements, Idnodes_S_air_interface = silex_lib_gmsh.ReadGmshElements(
        mesh_file+'_struc.msh', 2, 2)

    struc_nnodes = struc_nodes.shape[0]
    struc_nelem = struc_elements.shape[0]

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults2(
            results_file+'_struc_surface', struc_nodes, struc_elements, 2)

    if rank == 0:
        print("nnodes for structure=", struc_nnodes)
        print("nelem for structure=", struc_nelem)

    ##################################################################
    # compute level set
    ##################################################################

    tic = time.clock()

    # LS from a structure mesh
    LevelSet_from_Mesh,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)

    # LS from a simple analytic shape (sphere)
    lx3 = paraVal[0] #2.0 # Xc
    ly3 = paraVal[1] #2.0 # Yc
    R = paraVal[2] #1.0 # sphere radius
    #
    lz3 = 0.0 # fixed


    # analytic LS
    LevelSet=scipy.sqrt((fluid_nodes[:,0]-lx3)**2+(fluid_nodes[:,1]-ly3)**2+(fluid_nodes[:,2]-lz3)**2)-R
    #temprorary levelset gradients
    LevelSet_gradient_tmp=[]
    NameParaTmp=['X','Y','R']
    #Compute LS gradient according to Xc
    LevelSet_gradient_tmp.append((lx3-fluid_nodes[:,0])/(LevelSet+R))
    #Compute LS gradient according to Yc
    LevelSet_gradient_tmp.append((ly3-fluid_nodes[:,1])/(LevelSet+R))
    #Compute LS gradient according to R
    LevelSet_gradient_tmp.append(fluid_nodes[:,0]*0.-1)

    #load require levelSet gradients
    LevelSetGradient=[]
    NamePara=[]
    if len(gradValRequire)>0:
        for it in gradValRequire:
            LevelSetGradient.append(LevelSet_gradient_tmp[it])
            NamePara.append(NameParaTmp[it])

    #number of parameters 
    nbPara=len(NamePara)
  

    toc = time.clock()
    if rank == 0:
        print("time to compute level set:", toc-tic)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults2(
            results_file+'_struc_air_interface', struc_nodes, struc_elements, 2)
        #export levelset and levelset gradient
        dataW=list()
        dataW.append([[LevelSet],'nodal',1,'Level set'])
        itP=0
        for iN in NamePara:
            dataW.append([[LevelSetGradient[itP]],'nodal',1,'Level Set Grad '+iN])
            itP=itP+1

        silex_lib_gmsh.WriteResults2(results_file+'_LS_data', fluid_nodes,
                                    fluid_elements1, 4, dataW)

    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.clock()

    LSEnrichedElements, NbLSEnrichedElements = silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(
        fluid_elements1, LevelSet)
    LSEnrichedElements = LSEnrichedElements[list(range(NbLSEnrichedElements))]
    silex_lib_gmsh.WriteResults2(results_file+'_LS_enriched_elements',
                                fluid_nodes, fluid_elements1[LSEnrichedElements], 4)
    # EnrichedElements=LSEnrichedElements#[EnrichedElements-1]
    LSEnrichednodes = scipy.unique(fluid_elements1[LSEnrichedElements])

    tmp = []
    for i in LSEnrichednodes:
        for j in range(4):
            tmpp = scipy.where(fluid_elements1[:, j] == i)[0]
            for k in range(len(tmpp)):
                tmp.append(tmpp[k])
    # tmp.append(scipy.where(fluid_elements1[:,1]==i))
    # tmp.append(scipy.where(fluid_elements1[:,2]==i))
    # tmp.append(scipy.where(fluid_elements1[:,3]==i))

    tmp = scipy.unique(scipy.array(tmp))
    # tmp1,elttest0,tmp2=scipy.intersect1d(fluid_elements1[:,0],LSEnrichednodes,return_indices=True)
    # silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements_test0',fluid_nodes,fluid_elements1[tmp],4)
    #[75804, 97252, 97253,34973, 93135, 93137, 93248,83787, 93136,93525]
    EnrichedElements0, NbEnrichedElements = silex_lib_xfem_acou_tet4.getsurfenrichedelements(
        struc_nodes, struc_elements, fluid_nodes, fluid_elements1[tmp])
    EnrichedElements0 = scipy.unique(
        EnrichedElements0[list(range(NbEnrichedElements))])
    EnrichedElements0 = EnrichedElements0-1
    EnrichedElements = tmp[EnrichedElements0]
    toc = time.clock()
    if rank == 0:
        print("time to find enriched elements:", toc-tic)

    tic = time.clock()

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults2(
            results_file+'_enriched_elements', fluid_nodes, fluid_elements1[EnrichedElements], 4)
        LS_moins_enriched = scipy.setdiff1d(LSEnrichedElements, EnrichedElements)
        enriched_moins_LS = scipy.setdiff1d(EnrichedElements, LSEnrichedElements)
        silex_lib_gmsh.WriteResults2(
            results_file+'_LS_moins_enriched', fluid_nodes, fluid_elements1[LS_moins_enriched], 4)
        silex_lib_gmsh.WriteResults2(
                results_file+'_enriched_moins_LS', fluid_nodes, fluid_elements1[enriched_moins_LS], 4)
    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.clock()

    IIf, JJf, Vffk, Vffm = silex_lib_xfem_acou_tet4.globalacousticmatrices(
        fluid_elements1, fluid_nodes, celerity, rho)

    KFF = scipy.sparse.csc_matrix(
        (Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix(
        (Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = list(range(fluid_ndof))


    ##################################################################
    # Compute Heaviside enrichment
    ##################################################################
    tic = time.clock()

    Enrichednodes = scipy.unique(fluid_elements1[EnrichedElements])

    IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(
        fluid_elements1, fluid_nodes, LevelSet, celerity, rho)

    KAA = scipy.sparse.csc_matrix(
        (Vaak, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
    MAA = scipy.sparse.csc_matrix(
        (Vaam, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
    KAF = scipy.sparse.csc_matrix(
        (Vafk, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))
    MAF = scipy.sparse.csc_matrix(
        (Vafm, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofA = Enrichednodes-1

    toc = time.clock()
    if rank == 0:
        print("time to compute Heaviside enrichment:", toc-tic)

    ##################################################################
    # Construct the whole system
    #################################################################

    K = scipy.sparse.construct.bmat([
        [fluid_damping*KFF[SolvedDofF, :][:, SolvedDofF], fluid_damping*KAF[SolvedDofF, :][:, SolvedDofA]],
        [fluid_damping*KAF[SolvedDofA, :][:, SolvedDofF], fluid_damping*KAA[SolvedDofA, :][:, SolvedDofA]]])

    M = scipy.sparse.construct.bmat([
        [MFF[SolvedDofF, :][:, SolvedDofF], MAF[SolvedDofF, :][:, SolvedDofA]],
        [MAF[SolvedDofA, :][:, SolvedDofF], MAA[SolvedDofA, :][:, SolvedDofA]]])

    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 1
    UF = scipy.zeros(2*fluid_ndof, dtype=float)
    UF[9-1] = 3.1250E-05

    SolvedDof = scipy.hstack([SolvedDofF, SolvedDofA+fluid_ndof])

 
    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################
    #print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
    dK=list()
    dM=list()
    for itP in range(0,nbPara):
        print(' Build gradient matrices for parameter '+NamePara[itP])
        #
        IIf,JJf,Vfak_gradient,Vfam_gradient=\
            silex_lib_xfem_acou_tet4.globalacousticgradientmatrices(fluid_nodes,\
                fluid_elements1,LevelSet,celerity,rho,LevelSetGradient[itP])
        dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        #build full stiffness and mass gradient matrices
        dK.append(scipy.sparse.construct.bmat( [
                    [None,fluid_damping*dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                    [fluid_damping*dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]] ))
        dM.append(scipy.sparse.construct.bmat( [
                    [None,dMFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                    [dMFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]] ))

    ##############################################################
    # FRF computation
    ##############################################################

    Flag_frf_analysis = 1
    frequencies = []
    frf = []
    frfgradient=list()
    for it in range(0,nbPara):
        frfgradient.append([])

    if (Flag_frf_analysis == 1):
        print("Proc. ", rank, " / time at the beginning of the FRF:", time.ctime())

        if rank == 0:
            print('nb of total dofs: ', len(SolvedDofF)+len(SolvedDofA))

        press_save = []
        enrichment_save = []
        uncorrectedpress_save = []
        disp_save = []
        dpress_save=list()
        for it in range(0,nbPara):
            dpress_save.append([])

        #extract frequencies for the associated processors
        freqCompute=listFreqPerProc[:,rank]
        freqCompute=freqCompute[freqCompute>0]
        it=0
        itmax=len(freqCompute)         
        for freq in freqCompute:
            it=it+1
            #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2*scipy.pi*freq

            print("Freq. step ",it," proc number", rank, "frequency=", freq)

            tic = time.clock()

            F = scipy.array(omega**2*UF[SolvedDof], dtype='c16')

            if rank>=0:
                print(K)
                print(M)
                print(omega)
                sol = mumps.spsolve(scipy.sparse.csc_matrix(
                    K-(omega**2)*M, dtype='c16'), F)
                sol
                #sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(
                #     K-(omega**2)*M, dtype='c16'), F)
            stop
            ## pressure field without enrichment
            press1 = scipy.zeros((fluid_ndof), dtype=complex)
            press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))]
            ## enrichment field
            enrichment = scipy.zeros((fluid_nnodes), dtype=complex)
            enrichment[SolvedDofA] = sol[list(
                range(len(SolvedDofF), len(SolvedDofF)+len(SolvedDofA)))]
            ## correction of the pressure field with enrichment
            CorrectedPressure = press1.copy()
            CorrectedPressure[SolvedDofA] = press1[SolvedDofA] + \
                enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
            ## compute and store FRF on the test volume
            # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
            frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(
                fluid_elements5, fluid_nodes, press1, enrichment, LevelSet, LevelSet*0-1.0))

            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(CorrectedPressure)
                enrichment_save.append(enrichment)
                uncorrectedpress_save.append(press1)

            #####################
            #####################
            ######################
            #####################
            Dpress_Dtheta = scipy.zeros([fluid_ndof,nbPara],dtype=complex)
            DCorrectedPressure_Dtheta=scipy.array(Dpress_Dtheta)
            Denrichment_Dtheta = scipy.zeros([fluid_ndof,nbPara],dtype=complex)
            #####################
            #####################
            ## compute gradients
            for itP in range(0,nbPara):
                ## solve gradient problem
                tmp=-(dK[itP]-(omega**2)*dM[itP])*sol
                Dsol_Dtheta_RAW = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
                # Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
                #####################
                #####################
                ## gradient of the pressure field without enrichment
                Dpress_Dtheta[SolvedDofF,itP] = Dsol_Dtheta_RAW[list(range(len(SolvedDofF)))].copy()
                #####################
                #####################
                ## gradient of the enrichment field
                Denrichment_Dtheta[SolvedDofA,itP]= Dsol_Dtheta_RAW[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
                #####################
                #####################
                #compute the corrected gradient pressure field (via enrichment)
                DCorrectedPressure_Dtheta[:,itP]=scipy.array(Dpress_Dtheta[:,itP].copy())
                DCorrectedPressure_Dtheta[SolvedDofA,itP]=DCorrectedPressure_Dtheta[SolvedDofA,itP].T+ \
                    scipy.array(Denrichment_Dtheta[SolvedDofA,itP]*scipy.sign(LevelSet[SolvedDofA]).T)
                #####################
                #####################
                #store gradients
                frfgradient[itP].append(\
                    silex_lib_xfem_acou_tet4.computegradientcomplexquadratiquepressure(\
                        fluid_elements5,\
                            fluid_nodes,\
                                CorrectedPressure+0j,\
                                    Dpress_Dtheta[:,itP]+0j,\
                                        LevelSet))                
                #####################
                #####################
                dpress_save[itP].append(DCorrectedPressure_Dtheta[:,itP].copy())

        frfsave=[frequencies,frf,frfgradient]
        if rank!=0:
            comm.send(frfsave, dest=0, tag=11)

        print("Proc. ", rank, " / time at the end of the FRF:", time.ctime())

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW=list()
            #prepare pressure field
            dataW.append([scipy.real(press_save),'nodal',1,'pressure (real)'])
            dataW.append([scipy.imag(press_save),'nodal',1,'pressure (imaginary)'])
            dataW.append([scipy.absolute(press_save),'nodal',1,'pressure (norm)'])
            #prepare gradient pressure field
            itG=0
            for itP in NamePara:
                dataW.append([scipy.real(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (real)'])
                dataW.append([scipy.imag(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (imaginary)'])
                dataW.append([scipy.absolute(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (norm)'])
                itG=itG+1
            print("Write pressure field and gradients in msh file")
            silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',
                                        fluid_nodes, fluid_elements1, 4,dataW)
            print(">>> Done!!")

            #export results with discontinuities on .pos files
            varExport=scipy.vstack(uncorrectedpress_save).transpose()
            varExportC=scipy.vstack(press_save).transpose()
            varExportB=scipy.vstack(enrichment_save).transpose()
            print("Write pressure field in pos file")
            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,varExport.real,varExportB.real,'press_plus.pos')
            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,varExport.real,-varExportB.real,'press_moins.pos')
            print(">>> Done!!")
            #
            itG=0
            for itP in NamePara:
                varExport=scipy.vstack(dpress_save[itG]).transpose()
                print("Write gradient of pressure field in pos file (",itP,")")
                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,varExport.real,varExportB.real,'Gpress_plus_'+itP+'.pos')
                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,varExport.real,-varExportB.real,'Gpress_moins_'+itP+'.pos')
                print(">>> Done!!")
                itG=itG+1

        #####################
        #####################
        # save the FRF problem
        Allfrequencies=scipy.zeros(nbStep)
        Allfrf=scipy.zeros(nbStep)
        Allfrfgradient=scipy.zeros([nbStep,nbPara])
        k=0
        if rank==0:
            for i in range(nproc):
                if i==0:
                    data=frfsave
                    print(data)
                else:
                    print(i)
                    data=comm.recv(source=i, tag=11)
                    #data=data_buffer

                for j in range(len(data[0])):
                    Allfrequencies[k]=data[0][j]
                    Allfrf[k]=data[1][j]
                    for itP in range(0,nbPara):
                        Allfrfgradient[k,itP]=data[2][itP][j]
                    k=k+1
            #####################
            #####################zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            IXsort=scipy.argsort(Allfrequencies)
            AllfreqSorted=scipy.zeros(nbStep)
            AllfrfSorted=scipy.zeros(nbStep)
            AllfrfgradientSorted=scipy.zeros([nbStep,nbPara])
            for itS in range(0,nbStep):
                AllfreqSorted[itS]=Allfrequencies[IXsort[itS]]
                AllfrfSorted[itS]=Allfrf[IXsort[itS]]
                for itP in range(0,nbPara):
                    AllfrfgradientSorted[itS,itP]=Allfrfgradient[IXsort[itS],itP]

            #Allfrequencies, Allfrf,Allfrfgradient = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            Allfrfsave=list()
            Allfrfsave.append(AllfreqSorted)
            Allfrfsave.append(AllfrfSorted)
            for itP in range(0,nbPara):
                Allfrfsave.append(AllfrfgradientSorted[:,itP])

            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            print(Allfrfsave)
            f.close()
            #####################
            #####################
            #save on mat file
            scipy.io.savemat(results_file+'_results.mat',mdict={'AllFRF': Allfrfsave})
            scipy.io.savemat(results_file_ini+'results.mat',mdict={'AllFRF': Allfrfsave})
            #####################
            #####################
            return Allfrfsave


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#function for dealing with options
def manageOpt(argv,dV):
    #load default values
    freqMin     = dV.freqMin
    freqMax     = dV.freqMax
    nbStep      = dV.nbStep
    paraVal  = scipy.array(dV.paraVal)
    gradCompute = scipy.array(dV.gradCompute)
    #caseDefine = dV.caseDef
    
    #load info from MPI
    nbProc,rank,comm=mpiInfo()
    #load options
    opts,args = getopt.getopt(argv,"p:s:F:f:hp:c:")
    for opt,arg in opts:
        if opt == "-s":
            nbStep  = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-p":
            tmp = scipy.array(arg.split(','),dtype=scipy.float32)
            paraVal=tmp
        elif opt == "-c":
            caseDefine=str(arg)
        elif opt == "-g":
            tmp = scipy.array(arg.split(','),dtype=scipy.int)
            gradCompute=tmp
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    print ("Number of processors: ",nbProc)
    print ("Number of frequency steps: ",nbStep)
    print ("Maximum frequency: ",freqMax)
    print ("Minimum frequency: ",freqMin)
    print ("Components of grad: ",gradCompute)
    #print ("Case: ",caseDefine)
    it=0
    for itP in paraVal:
        print ('Parameter num '+str(it)+': '+str(itP))
        it=it+1
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,gradCompute)#,caseDefine)

#usage definition
def usage():
    dV=defaultV
    print("Usage: ",sys.argv[0],"-psFfh [+arg]")
    print("\t -p : number of processors (default value ",dV.nbProc,")")
    print("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    print("\t -F : maximum frequency (default value ",dV.freqMax,")")
    print("\t -f : minimum frequency (default value ",dV.freqMin,")")
    print("\t -g : Components of grad (default value ",dV.gradCompute,")")

#default values
class defaultV:
    freqMin     = 10.0
    freqMax     = 200.0
    nbStep      = 5
    paraVal   = [1.5,1.,1.]
    gradCompute = [0,1,2]
    nbProc=1
    #caseDef= 'thick_u'

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)

