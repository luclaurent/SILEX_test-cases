###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# Libraries
###########################################################
import string
import time
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io
import getopt
import os
from shutil import copyfile
from pathlib import Path
from loguru import logger

import pylab as pl
import pickle


import pymumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
acousticsFEM = silex_lib_fem.LinearAcousticsTET4()
acousticsXFEM = silex_lib_xfem.LinearAcousticsTET4()
# structureFEM = silex_lib_fem.DKT()


from mpi4py import MPI

comm = MPI.COMM_WORLD

nproc = comm.Get_size()
rank = comm.Get_rank()
log_format = (
    "<cyan> R{extra[rank]}</cyan> |"
    "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> | "
    "<level>{message}</level>"
)

logger.remove()
logger.configure(extra={"rank": 0})  # Default values
logger.add(
    sys.stdout,
    level="DEBUG",
    format=log_format,
    colorize=True,
    backtrace=True,
    diagnose=True,
)
logger = logger.bind(rank=rank)

# mpirun -np 2 python Main_xfem.py
logger.info("START")


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
    listFreq = np.zeros((nbFreqProc+varCase, nbProc))
    listAllFreq = np.linspace(freqInit, freqEnd, nbStep)
    # logger.info(np.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1

    # logger.info(listFreq)
    return listFreq

#load structure mesh fo values of parameters
def buildStructMesh(fileOrig,destFile,paraVal):
    #copy original file to the used one
    copyfile(fileOrig+'.geo',destFile+'.geo')
    #change value of parameters in the new file
    for key,value in enumerate(paraVal):
        oldText="<val##"+str(key)+">"
        newText='%g'%value
        #logger.info(oldText)
        #logger.info(newText) 
        cmdSed="sed -i 's/"+oldText+"/"+newText+"/g' "+destFile+'.geo'
        #logger.info(cmdSed)
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

    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##################################################")
    logger.info("##    Start SILEX vibro-acoustics computation   ##")
    if len(gradValRequire)>0:
        logger.info("##          (with gradients computation)        ##")
    logger.info("##################################################")
    logger.info("##################################################")

    # load 3D geometry
    orig_mesh_file = 'geom/cavity_acou3D_struc_3D_v3_para'
    mesh_file = 'geom/cavity_acou3D_struc_3D_v3'
    results_file_ini = 'results/cavity_acou3D_struc_3D_v3'
    cwd = Path(__file__).resolve().parent

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
    logger.info(results_file)

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

    mesh = msh.mshReader(cwd / (mesh_file + "_air.msh"))

    fluid_nodes = mesh.getNodes()
    fluid_elements1 = mesh.getElements(tag=1)["TET4"]  # air, cavity
    fluid_elements5 = mesh.getElements(tag=5)["TET4"]  # air, control volum

    IdNodes1 = np.unique(fluid_elements1.flatten())
    IdNodes5 = np.unique(fluid_elements5.flatten())

    fluid_nnodes = fluid_nodes.shape[0]
    fluid_nelem1 = fluid_elements1.shape[0]
    fluid_nelem5 = fluid_elements5.shape[0]

    fluid_nnodes1 = IdNodes1.shape[0]
    fluid_nnodes5 = IdNodes5.shape[0]

    fluid_ndof = fluid_nnodes

    if rank == 0:
        logger.info("Number of nodes:", fluid_nnodes)
        logger.info("Number of elements in air:", fluid_nelem1)
        logger.info("Number of nodes in air:", fluid_nnodes1)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
        cwd / (results_file + "_air_cavity_Mesh1.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1},
        )
        msh2.mshWriter(
            cwd / (results_file + "_air_controlled_volume_Mesh5.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements5},
        )

    # ##############################################################
    # # Load structure mesh
    # ##############################################################
    # #change parameters values and build mesh
    # buildStructMesh(orig_mesh_file+'_struc',mesh_file+'_struc',paraVal)

    # struc_nodes = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh', 3)
    # struc_elements, Idnodes_S_air_interface = silex_lib_gmsh.ReadGmshElements(
    #     mesh_file+'_struc.msh', 2, 2)

    # struc_nnodes = struc_nodes.shape[0]
    # struc_nelem = struc_elements.shape[0]

    # if (flag_write_gmsh_results == 1) and (rank == 0):
    #     silex_lib_gmsh.WriteResults2(
    #         results_file+'_struc_surface', struc_nodes, struc_elements, 2)

    # if rank == 0:
    #     logger.info("nnodes for structure=", struc_nnodes)
    #     logger.info("nelem for structure=", struc_nelem)

    ##################################################################
    # compute level set
    ##################################################################

    tic = time.process_time()

    # LS from a structure mesh
    # LevelSet_from_Mesh,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)

    # LS from a simple analytic shape (sphere)
    lx3 = paraVal[0] #2.0 # Xc
    ly3 = paraVal[1] #2.0 # Yc
    lz3 = paraVal[2] #0.0 # YZ
    R = paraVal[3] #1.0 # sphere radius
    #
    logger.info("Parameters values")
    logger.info("Xc ",lx3," Yc ",ly3," Zc ",lz3," R ",R)
    # analytic LS
    LevelSet=np.sqrt((fluid_nodes[:,0]-lx3)**2+(fluid_nodes[:,1]-ly3)**2+(fluid_nodes[:,2]-lz3)**2)-R
    #temprorary levelset gradients
    LevelSet_gradient_tmp=[]
    NameParaTmp=['X','Y','Z','R']
    #Compute LS gradient according to Xc
    LevelSet_gradient_tmp.append((lx3-fluid_nodes[:,0])/(LevelSet+R))
    #Compute LS gradient according to Yc
    LevelSet_gradient_tmp.append((ly3-fluid_nodes[:,1])/(LevelSet+R))
    #Compute LS gradient according to Zc
    LevelSet_gradient_tmp.append((lz3-fluid_nodes[:,2])/(LevelSet+R))
    #Compute LS gradient according to R
    LevelSet_gradient_tmp.append(fluid_nodes[:,0]*0.-1.)

    #load require levelSet gradients
    LevelSetGradient=[]
    NamePara=[]
    if len(gradValRequire)>0:
        for it in gradValRequire:
            LevelSetGradient.append(LevelSet_gradient_tmp[it])
            NamePara.append(NameParaTmp[it])

    #number of parameters 
    nbPara=len(NamePara)
  

    toc = time.process_time()
    if rank == 0:
        logger.info("time to compute level set: {}".format(toc-tic))

    if (flag_write_gmsh_results == 1) and (rank == 0):
        # silex_lib_gmsh.WriteResults2(
        #     results_file+'_struc_air_interface', struc_nodes, struc_elements, 2)
        #export levelset and levelset gradient
        dataW=list()
        dataW.append({'data':LevelSet,
                      'type':'nodal',
                      'name':'Level set'})
        itP=0
        for iN in NamePara:
            dataW.append({'data':LevelSetGradient[itP],
                         'type':'nodal',
                         'name':'Level Set Grad '+iN})
            itP=itP+1

        msh2.mshWriter(
            cwd / (results_file + "_LS_data.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1},
            fields=dataW,
            append=True,
        )


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.process_time()

    LSEnrichedElements = acousticsXFEM.getEnrichedElements(
        fluidElements=fluid_elements1, levelset=LevelSet
    )
    msh2.mshWriter(
        cwd / (results_file + "_LSenriched_elements.msh"),
        fluid_nodes,
        {"type": "TET4", "connectivity": fluid_elements1[LSEnrichedElements]},
    )
    # EnrichedElements=LSEnrichedElements#[EnrichedElements-1]
    LSEnrichednodes = np.unique(fluid_elements1[LSEnrichedElements])

    tmp = []
    for i in LSEnrichednodes:
        for j in range(4):
            tmpp = np.where(fluid_elements1[:, j] == i)[0]
            for k in range(len(tmpp)):
                tmp.append(tmpp[k])
    # tmp.append(np.where(fluid_elements1[:,1]==i))
    # tmp.append(np.where(fluid_elements1[:,2]==i))
    # tmp.append(np.where(fluid_elements1[:,3]==i))

    tmp = np.unique(np.array(tmp))
    # tmp1,elttest0,tmp2=scipy.intersect1d(fluid_elements1[:,0],LSEnrichednodes,return_indices=True)
    # silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements_test0',fluid_nodes,fluid_elements1[tmp],4)
    #[75804, 97252, 97253,34973, 93135, 93137, 93248,83787, 93136,93525]
    # EnrichedElements0, NbEnrichedElements = silex_lib_xfem_acou_tet4.getsurfenrichedelements(
    #     struc_nodes, struc_elements, fluid_nodes, fluid_elements1[tmp])
    # EnrichedElements0 = np.unique(
    #     EnrichedElements0[list(range(NbEnrichedElements))])
    # EnrichedElements0 = EnrichedElements0-1
    # EnrichedElements = tmp[EnrichedElements0]

    EnrichedElements=LSEnrichedElements


    toc = time.process_time()
    if rank == 0:
        logger.info("time to find enriched elements: {}".format(toc-tic))

    tic = time.process_time()

    LS_moins_enriched = np.setdiff1d(LSEnrichedElements, EnrichedElements)
    enriched_moins_LS = np.setdiff1d(EnrichedElements, LSEnrichedElements)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_enriched_elements.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1[EnrichedElements]},
        )
        msh2.mshWriter(
            cwd / (results_file + "_LS_moins_enriched.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1[LS_moins_enriched]},
        )
        msh2.mshWriter(
            cwd / (results_file + "_enriched_moins_LS.msh"),
            fluid_nodes,
            {"type": "TET4", "connectivity": fluid_elements1[enriched_moins_LS]},
        )


    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.process_time()

    IIf, JJf, Vffk, Vffm = acousticsFEM.getMatrices(
        fluid_nodes, fluid_elements1, [celerity, rho]
    )

    KFF = scipy.sparse.csc_matrix(
        (Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix(
        (Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = list(range(fluid_ndof))


    ##################################################################
    # Compute Heaviside enrichment
    ##################################################################
    tic = time.process_time()

    Enrichednodes = np.unique(fluid_elements1[EnrichedElements])

   
    IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = acousticsXFEM.getMatrices(
        fluid_nodes, fluid_elements1, LevelSet, [celerity, rho]
    )

    KAA = scipy.sparse.csc_matrix(
        (Vaak, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
    MAA = scipy.sparse.csc_matrix(
        (Vaam, (IIaa, JJaa)), shape=(fluid_ndof, fluid_ndof))
    KAF = scipy.sparse.csc_matrix(
        (Vafk, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))
    MAF = scipy.sparse.csc_matrix(
        (Vafm, (IIaf, JJaf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofA = Enrichednodes-1

    toc = time.process_time()
    if rank == 0:
        logger.info("time to compute Heaviside enrichment: {}".format(toc-tic))

    ##################################################################
    # Construct the whole system
    #################################################################

    K = scipy.sparse.bmat([
        [fluid_damping*KFF[SolvedDofF, :][:, SolvedDofF], fluid_damping*KAF[SolvedDofF, :][:, SolvedDofA]],
        [fluid_damping*KAF[SolvedDofA, :][:, SolvedDofF], fluid_damping*KAA[SolvedDofA, :][:, SolvedDofA]]])

    M = scipy.sparse.bmat([
        [MFF[SolvedDofF, :][:, SolvedDofF], MAF[SolvedDofF, :][:, SolvedDofA]],
        [MAF[SolvedDofA, :][:, SolvedDofF], MAA[SolvedDofA, :][:, SolvedDofA]]])

    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 1
    UF = np.zeros(2*fluid_ndof, dtype=float)
    UF[9-1] = 3.1250E-05

    SolvedDof = np.hstack([SolvedDofF, SolvedDofA+fluid_ndof])

 
    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################
    #logger.info(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
    dK=list()
    dM=list()
    for itP in range(0,nbPara):
        logger.info(' Build gradient matrices for parameter '+NamePara[itP])
        #
        IIf, JJf, Vfak_gradient, Vfam_gradient = acousticsXFEM.getGradientMatrices(
            nodes=fluid_nodes,
            elements=fluid_elements1,
            levelset=LevelSet,
            levelsetGradient=LevelSetGradient[itP],
            material=[celerity, rho],
        )

        dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        #build full stiffness and mass gradient matrices
        dK.append(scipy.sparse.bmat( [
                    [None,fluid_damping*dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                    [fluid_damping*dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]] ))
        dM.append(scipy.sparse.bmat( [
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
        logger.info("Proc. {} / time at the beginning of the FRF: {}".format(rank, time.ctime()))

        if rank == 0:
            logger.info('nb of total dofs: {}'.format(len(SolvedDofF)+len(SolvedDofA)))

        press_save = []
        enrichment_save = []
        uncorrectedpress_save = []
        disp_save = []
        denrichment_save=[]
        duncorrectedpress_save=[]
        dpress_save=list()
        for it in range(0,nbPara):
            dpress_save.append([])
            denrichment_save.append([])
            duncorrectedpress_save.append([])

        #extract frequencies for the associated processors
        freqCompute=listFreqPerProc[:,rank]
        freqCompute=freqCompute[freqCompute>0]
        it=0
        itmax=len(freqCompute)         
        for freq in freqCompute:
            it=it+1
            #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega = 2*np.pi*freq

            logger.info("Freq. step  {}/{} - proc number {} - frequency= {}".format(it, itmax, rank, freq))

            tic = time.process_time()

            F = np.array(omega**2*UF[SolvedDof], dtype='c16')

            logger.info(rank)
            if rank>=0:
                 #logger.info(K)
                 #logger.info(M)
                 #logger.info(omega)
                #  sol = mumps.spsolve(K-(omega**2)*M, F,comm=None)#mycomm)
                 sol = scipy.sparse.linalg.spsolve(K-(omega**2)*M, F)
                 #sol = mumps.spsolve(scipy.sparse.coo_matrix( \
                 #   K-(omega**2)*M, dtype='complex'), F+0.j,comm=mycomm)
                 #sol
                 #sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(
                 #     K-(omega**2)*M, dtype='c16'), F)
            
            ## pressure field without enrichment
            press1 = np.zeros((fluid_ndof), dtype=complex)
            press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))].copy()
            ## enrichment field
            enrichment = np.zeros((fluid_nnodes), dtype=complex)
            enrichment[SolvedDofA] = sol[list(
                range(len(SolvedDofF), len(SolvedDofF)+len(SolvedDofA)))].copy()
            ## correction of the pressure field with enrichment
            CorrectedPressure =press1.copy() #np.zeros((fluid_ndof),dtype=complex) #press1.copy()
            CorrectedPressure[SolvedDofA] = press1[SolvedDofA] + \
                 enrichment[SolvedDofA]*np.sign(LevelSet[SolvedDofA])

            
            ## compute and store FRF on the test volume
            # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
            frf.append(
                acousticsXFEM.getQuadraticPressure(
                fluid_nodes,
                fluid_elements5,
                press1,
                enrichment,
                LevelSet,
                LevelSet * 0 - 1.0,
                )
            )

            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(CorrectedPressure.copy())
                enrichment_save.append(enrichment.copy())
                uncorrectedpress_save.append(press1.copy())
            
            #####################
            #####################
            ######################
            #####################
            Dpress_Dtheta = np.zeros([fluid_ndof,nbPara],dtype=complex)
            DCorrectedPressure_Dtheta=np.array(Dpress_Dtheta)
            Denrichment_Dtheta = np.zeros([fluid_ndof,nbPara],dtype=complex)
            #####################
            #####################
            ## compute gradients
            for itP in range(0,nbPara):
                ## solve gradient problem
                tmp=-(dK[itP]-(omega**2)*dM[itP])*sol
                # Dsol_Dtheta_RAW = mumps.spsolve(K-(omega**2)*M, tmp, comm=mycomm )
                Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve(K-(omega**2)*M, tmp)
                #Dsol_Dtheta_RAW = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
                # Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
                #####################
                #####################
                ## gradient of the pressure field without enrichment
                Dpress_Dtheta[SolvedDofF,itP] = Dsol_Dtheta_RAW[list(range(len(SolvedDofF)))].copy()
                #####################
                #####################
                ## gradient of the enrichment field
                Denrichment_Dtheta[SolvedDofA,itP]= Dsol_Dtheta_RAW[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))].copy()
                #####################
                #####################
                #compute the corrected gradient pressure field (via enrichment)
                DCorrectedPressure_Dtheta[:,itP]=np.array(Dpress_Dtheta[:,itP].copy())
                DCorrectedPressure_Dtheta[SolvedDofA,itP]=DCorrectedPressure_Dtheta[SolvedDofA,itP].T+ \
                    np.array(Denrichment_Dtheta[SolvedDofA,itP]*np.sign(LevelSet[SolvedDofA]).T)
                #####################
                #####################
                #store gradients
                frfgradient[itP].append(
                    acousticsXFEM.getGradientQuadraticPressure(
                        nodes=fluid_nodes,
                        elements=fluid_elements5,
                        pressureField=press1,
                        gradientPressureField=Dpress_Dtheta[:,itP],
                        levelset=LevelSet,
                    )
                )

                #####################
                #####################
                dpress_save[itP].append(DCorrectedPressure_Dtheta[:,itP].copy())
                denrichment_save[itP].append(Denrichment_Dtheta[:,itP].copy())
                duncorrectedpress_save[itP].append(Dpress_Dtheta[:,itP].copy())

        frfsave=[frequencies,frf,frfgradient]
        if rank!=0:
            comm.send(frfsave, dest=0, tag=11)

        logger.info("Proc. {} / time at the end of the FRF: {}".format(rank,time.ctime()))

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW=list()
            #prepare pressure field
            dataW.append({'data':np.real(press_save),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure (real)'})
            dataW.append({'data':np.imag(press_save),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure (imaginary)'})
            dataW.append({'data':np.absolute(press_save),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure (norm)'})
            #prepare grad{'data':ent pressure field
            itG=0
            for itP in NamePara:
                dataW.append({'data':np.real(dpress_save[itG]),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure gradient '+itP+' (real)'})
                dataW.append({'data':np.imag(dpress_save[itG]),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure gradient '+itP+' (imaginary)'})
                dataW.append({'data':np.absolute(dpress_save[itG]),'type':'nodal','nbsteps':len(freqCompute),'name':'pressure gradient '+itP+' (norm)'})
                itG=itG+1
            logger.info("Write pressure field and gradients in msh file")
            msh2.mshWriter(
                cwd / (results_file + str(rank) + "_results_fluid_frf.msh"),
                fluid_nodes,
                {"type": "TET4", "connectivity": fluid_elements1},
                fields=dataW,
                append=True)
            
            logger.info(">>> Done!!")

            #export results with discontinuities on .pos files
            varExport=np.vstack(uncorrectedpress_save).transpose()
            varExportC=np.vstack(press_save).transpose()
            varExportB=np.vstack(enrichment_save).transpose()
            logger.info("Write pressure field in pos file")
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=np.real(varExport),
                                        enrichedfield=np.real(varExportB),
                                        filename=cwd /'results/press_plus_real.pos',
                                        name='Pressure + Real',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=np.imag(varExport),
                                        enrichedfield=np.imag(varExportB),
                                        filename=cwd /'results/press_plus_imag.pos',
                                        name='Pressure + Imag',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=np.absolute(varExport),
                                        enrichedfield=np.absolute(varExportB),
                                        filename=cwd /'results/press_plus_abs.pos',
                                        name='Pressure + Abs',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=np.real(varExport),
                                        enrichedfield=-np.real(varExportB),
                                        filename=cwd /'results/press_moins_real.pos',
                                        name='Pressure - Real',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=np.imag(varExport),
                                        enrichedfield=-np.imag(varExportB),
                                        filename=cwd /'results/press_moins_imag.pos',
                                        name='Pressure - Imag',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=np.absolute(varExport),
                                        enrichedfield=-np.absolute(varExportB),
                                        filename=cwd /'results/press_moins_abs.pos',
                                        name='Pressure - Abs',
                                        nbFreq=len(freqCompute))
            os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
            

            # prepare final data
            objMesh = MeshField.MeshField(fluid_nodes, fluid_elements, LevelSet, LevelSetTangent)
            objMesh.addField(press,enrichment)
            # extract final data
            data = objMesh.getData(['nodes','mesh','levelset','levelset_tangent','fields'])



            #export data
            dataexport=list()
            dataexport.append(fluid_nodes)
            dataexport.append(fluid_elements1)
            dataexport.append(LevelSet)
            dataexport.append(varExport)
            dataexport.append(varExportB)
            with open('debug_export.pck','wb') as f:
                pickle.dump(dataexport, f)
            # logger.info(Allfrfsave)

            logger.info(">>> Done!!")
            #
            itG=0
            for key,itP in enumerate(NamePara):
                GvarExport=np.vstack(duncorrectedpress_save[key]).copy().transpose()
                GvarExportB=np.vstack(denrichment_save[key]).copy().transpose()
                logger.info("Write gradient of pressure field in pos file (",itP,")")
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=np.real(GvarExport),
                                        enrichedfield=np.real(GvarExportB),
                                        filename=cwd /('results/Gpress_plus_'+itP+'_real.pos'),
                                        name='Gpressure + '+itP+' Real',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=np.imag(GvarExport),
                                        enrichedfield=np.imag(GvarExportB),
                                        filename=cwd /('results/Gpress_plus_'+itP+'_imag.pos'),
                                        name='Gpressure + '+itP+' Imag',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=np.real(GvarExport),
                                        enrichedfield=-np.real(GvarExportB),
                                        filename=cwd /('results/Gpress_imag_'+itP+'_real.pos'),
                                        name='Gpressure - '+itP+' Real',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=np.imag(GvarExport),
                                        enrichedfield=-np.imag(GvarExportB),
                                        filename=cwd /('results/Gpress_imag_'+itP+'_imag.pos'),
                                        name='Gpressure - '+itP+' Imag',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                
                #     
                # gradPsquare=2*(np.real(GvarExport)*np.real(varExport)+np.imag(GvarExport)*np.imag(varExport))
                # gradPsquareB=2*(np.real(GvarExportB)*np.real(varExportB)+np.imag(GvarExportB)*np.imag(varExportB))
                # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,gradPsquare,gradPsquareB,'results/Gpress_plus_'+itP+'_dsquare.pos')                
                # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,gradPsquare,-gradPsquareB,'results/Gpress_moins_'+itP+'_dsquare.pos')
                #                          
                gradCalc=(np.real(GvarExport)*np.real(varExport)+np.imag(GvarExport)*np.real(varExport))/np.absolute(varExport)
                #deal with zeros values in enrichment field
                gradCalcB=np.zeros([fluid_nnodes,nbStep])
                gradCalcB[SolvedDofA,:]=(np.real(GvarExportB[SolvedDofA,:])*np.real(varExportB[SolvedDofA,:])+np.imag(GvarExportB[SolvedDofA,:])*np.imag(varExportB[SolvedDofA,:]))/np.absolute(varExportB[SolvedDofA,:])
                #remove inf value
                IX = np.absolute(varExport)==0.
                gradCalc[IX]=1
                gradCalcB[IX]=1
                #
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=LevelSet,
                                        field=gradCalc,
                                        enrichedfield=gradCalcB,
                                        filename=cwd /('Gpress_plus_'+itP+'_dabsolute.pos'),
                                        name='Gpressure + '+itP+' dAbs',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                acousticsXFEM.exportPosFile(nodes=fluid_nodes,
                                        elements=fluid_elements1,
                                        levelset=-LevelSet,
                                        field=gradCalc,
                                        enrichedfield=-gradCalcB,
                                        filename=cwd /('Gpress_moins_'+itP+'_dabsolute.pos'),
                                        name='Gpressure - '+itP+' dAbs',
                                        nbFreq=len(freqCompute))
                os.system('cd {}&&bzip2 -f *.pos'.format(cwd / 'results'))
                logger.info(">>> Done!!")
                itG=itG+1

        #####################
        #####################
        # save the FRF problem
        Allfrequencies=np.zeros(nbStep)
        Allfrf=np.zeros(nbStep)
        Allfrfgradient=np.zeros([nbStep,nbPara])
        k=0
        if rank==0:
            for i in range(nproc):
                if i==0:
                    data=frfsave
                    # logger.info(data)
                else:
                    # logger.info(i)
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
            IXsort=np.argsort(Allfrequencies)
            AllfreqSorted=np.zeros(nbStep)
            AllfrfSorted=np.zeros(nbStep)
            AllfrfgradientSorted=np.zeros([nbStep,nbPara])
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

            with open(cwd / (results_file + "_results.frf"), "wb") as f:
                pickle.dump(Allfrfsave, f)
            #####################
            #####################
            #save on mat file
            scipy.io.savemat(cwd / (results_file+'_results.mat'),mdict={'AllFRF': Allfrfsave})
            scipy.io.savemat(cwd / (results_file_ini+'results.mat'),mdict={'AllFRF': Allfrfsave})
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
    paraVal  = np.array(dV.paraVal)
    gradCompute = np.array(dV.gradCompute)
    #caseDefine = dV.caseDef
    
    #load info from MPI
    nbProc,rank,comm=mpiInfo()
    #load options
    opts,args = getopt.getopt(argv,"p:s:F:f:hp:c:g:")
    for opt,arg in opts:
        if opt == "-s":
            nbStep  = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-p":
            tmp = np.array(arg.split(','),dtype=scipy.float32)
            paraVal=tmp
        elif opt == "-c":
            caseDefine=str(arg)
        elif opt == "-g":
            tmp = np.array(arg.split(','),dtype=scipy.int32)
            gradCompute=tmp
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    logger.info("Number of processors: ",nbProc)
    logger.info("Parameters: ",paraVal)
    logger.info("Number of frequency steps: ",nbStep)
    logger.info("Maximum frequency: ",freqMax)
    logger.info("Minimum frequency: ",freqMin)
    logger.info("Components of grad: ",gradCompute)
    #logger.info("Case: ",caseDefine)
    it=0
    for itP in paraVal:
        logger.info('Parameter num '+str(it)+': '+str(itP))
        it=it+1
    logger.info("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,gradCompute,1)#,caseDefine)

#usage definition
def usage():
    dV=defaultV
    logger.info("Usage: ",sys.argv[0],"-psFfhg [+arg]")
    logger.info("\t -p : input parameters (default value ",dV.nbProc,")")
    logger.info("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    logger.info("\t -F : maximum frequency (default value ",dV.freqMax,")")
    logger.info("\t -f : minimum frequency (default value ",dV.freqMin,")")
    logger.info("\t -g : Components of grad (default value ",dV.gradCompute,")")

#default values
class defaultV:
    freqMin     = 10.0
    freqMax     = 150.0
    nbStep      = 5
    paraVal   = [1.,1.,0.5,0.8]
    gradCompute =  [0,1,2,3]
    nbProc=1
    #caseDef= 'thick_u'

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)

