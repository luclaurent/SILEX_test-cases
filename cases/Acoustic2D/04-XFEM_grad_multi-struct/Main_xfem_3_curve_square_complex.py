import string
import time
from pathlib import Path
import numpy as np
from loguru import logger
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle
import getopt

# import mumps

import sys
from meshRW import msh, msh2
from SILEXlib import silex_lib_fem, silex_lib_xfem
from SILEXlib import MeshField

# load classes
objFEM = silex_lib_fem.LinearAcousticsTRI3()
objXFEM = silex_lib_xfem.LinearAcousticsTRI3()

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
    nproc=comm.Get_size()
    rank = comm.Get_rank()
    return nproc,rank,comm


class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

def computeFreqPerProc(nbStep,nbProc,freqInit,freqEnd):
    #compute integer number of freq per proc and remaining steps
    nbFreqProc=nbStep // nbProc
    nbFreqProcRemain=nbStep % nbProc
    #compute frequencies steps
    varCase=1
    if nbFreqProcRemain==0:
        varCase=0
    listFreq=np.zeros((nbFreqProc+varCase,nbProc))
    listAllFreq=np.linspace(freqInit,freqEnd,nbStep)
    #logger.info(np.linspace(freqInit,freqEnd,nbStep))
    #build array of frequencies
    itF=0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP<nbStep:
                listFreq[itC,itP]=listAllFreq[itF]
                itF += 1

    #logger.info(listFreq)
    return listFreq

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py

# mpirun -np 2 python Main_xfem_2.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

def RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,positionStructX,positionStructY,radiusStruct):
    # parallepipedic cavity with plane structure
    mesh_file='geom/xfem3_square'
    results_file_ini='results_square/xfem_3_'
    cwd = Path(__file__).resolve().parent

    listFreqPerProc = computeFreqPerProc(nbStep,nbProc,freqMin,freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity=340.0
    rho=1.2
    fluid_damping=(1+0.01j)

    #freq_ini     = 100.0
    #freq_end     = 500.0
    #nb_freq_step_per_proc=400 # 50 pour 8 proc.

    nproc=comm.Get_size()
    rank = comm.Get_rank()

    #nb_freq_step = nb_freq_step_per_proc*nproc
    #deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

    flag_write_gmsh_results=1

    flag_edge_enrichment=0
    #flag_edge_enrichment=1

    #freq_comparaison = 210.0

    h_struc=2.0

    ValpStructX=positionStructX*10000.
    ValpStructY=positionStructY*10000.
    ValpStructR=radiusStruct*10000.

    file_extension=str(ValpStructX)[0:5]+'_'+str(ValpStructY)[0:5]+'_'+str(ValpStructR)[0:4]
    results_file=results_file_ini+file_extension
    logger.info(results_file)

    # center and radius of half-circle
    x_pos_struc=positionStructX
    y_pos_struc=positionStructY
    radius_hcircle=radiusStruct
    logger.info(x_pos_struc)
    logger.info(y_pos_struc)
    logger.info(radius_hcircle)

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

    # start reading mesh
    mesh = msh.mshReader(cwd / (mesh_file + "_fluid.msh"))

    fluid_nodes = mesh.getNodes()[:, 0:2]
    fluid_elements = mesh.getElements(tag=1)["TRI3"]
    fluid_elements5 = mesh.getElements(tag=5)["TRI3"]
    Idnodes = np.unique(fluid_elements.flatten())
    Idnodes5 = np.unique(fluid_elements5.flatten())

    fluid_nnodes = fluid_nodes.shape[0]
    fluid_nelem = fluid_elements.shape[0]
    fluid_ndof = len(np.unique(fluid_elements.flatten()))

    fluid_elements_boun = mesh.getElements(tag=2)["LIN2"]
    IdnodeS2 = np.unique(fluid_elements_boun.flatten())

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_fluid_mesh.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements},
        )
        msh2.mshWriter(
            cwd / (results_file + "_control_volume_fluid_mesh.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements5},
        )
        msh2.mshWriter(
            cwd / (results_file + "_fluid_boundary.msh"),
            fluid_nodes,
            {"type": "LIN2", "connectivity": fluid_elements_boun},
        )
    logger.info("nnodes for fluid= {}".format(fluid_nnodes))
    logger.info("nelem for fluid= {}".format(fluid_nelem))
    logger.info("nelem for control volume= {}".format(fluid_elements5.shape[0]))





    ##################################################################
    # compute level set
    ##################################################################

    tic = time.process_time()

    # create coordinates of nodes of the structure (half circle)
    nbNodesHC=50
    #thetaHC=np.linspace(0,np.pi,nbNodesHC)
    #thetaHC=np.linspace(-np.pi/2,np.pi/2,nbNodesHC)
    thetaHC=np.linspace(np.pi/2,3*np.pi/2,nbNodesHC)

    xNodesHC=x_pos_struc-radius_hcircle*np.sin(thetaHC)
    yNodesHC=y_pos_struc-radius_hcircle*np.cos(thetaHC)
    struc_nodes=np.vstack([xNodesHC,yNodesHC]).transpose()

    lA=np.linspace(1,nbNodesHC-1)
    lB=np.linspace(2,nbNodesHC)
    struc_elements=np.vstack([lB,lA]).transpose()



    LevelSet=np.sqrt((fluid_nodes[:,0]-x_pos_struc)**2+(fluid_nodes[:,1]-y_pos_struc)**2)-radius_hcircle
    LevelSetTangent=-(fluid_nodes[:,1]-y_pos_struc)

    # level set gradient with respect to parameters
    LevelSet_gradient_X=-(fluid_nodes[:,0]-x_pos_struc)/(LevelSet+radius_hcircle)
    LevelSet_gradient_Y=-(fluid_nodes[:,1]-y_pos_struc)/(LevelSet+radius_hcircle)
    LevelSet_gradient_R=-np.ones(fluid_nnodes)

    if (flag_write_gmsh_results==1) and (rank==0):
        msh2.mshWriter(
            cwd / (results_file + "_level_sets.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements},
            fields=[
                {
                    "data": LevelSet,
                    "type": "nodal",
                    "dim": 1,
                    "name": "levelset",
                },
                {
                    "data": LevelSetTangent,
                    "type": "nodal",
                    "dim": 1,
                    "name": "levelset_tangent",
                },
                {
                    "data": LevelSet_gradient_X,
                    "type": "nodal",
                    "dim": 1,
                    "name": "levelset_tangent /X",
                },
                {
                    "data": LevelSet_gradient_Y,
                    "type": "nodal",
                    "dim": 1,
                    "name": "levelset_tangent /Y",
                },
                {
                    "data": LevelSet_gradient_,
                    "type": "nodal",
                    "dim": 1,
                    "name": "levelset_tangent /",
                }
            ],
            append=True,
        )

    toc = time.process_time()
    logger.info("time to compute level set: {}".format(toc-tic))


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.process_time()

    
    struc_boun=np.array([3])

    msh2.mshWriter(
        cwd / (results_file + "_struc_mesh.msh"),
        struc_nodes,
        {"type": "LIN2", "connectivity": struc_elements},
    )

    EnrichedElements = objXFEM.getEnrichedElements(
        fluid_nodes, fluid_elements, struc_nodes, struc_elements
    )

    toc = time.process_time()
    logger.info("time to find surface enriched elements: {}".format(toc-tic))

    if (flag_write_gmsh_results==1) and (rank==0):
        msh2.mshWriter(
            cwd / (results_file + "_enriched_elements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[EnrichedElements]},
        )

    tic = time.process_time()

    EdgeEnrichedElements = objXFEM.getEdgeEnrichedElements(
        fluid_nodes, fluid_elements, struc_nodes, struc_boun
    )

    toc = time.process_time()
    logger.info("time to find edge enriched elements: {}".format(toc-tic))

    if (flag_write_gmsh_results==1) and (rank==0):
        msh2.mshWriter(
            cwd / (results_file + "_edge_enriched_elements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[EdgeEnrichedElements]},
        )

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################
    tic = time.process_time()

    IIf, JJf, Vffk, Vffm = objFEM.getMatrices(
        fluid_nodes, fluid_elements, [celerity, rho]
    )

    KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    SolvedDofF=np.setdiff1d(list(range(fluid_ndof)),IdnodeS2-1)
    #SolvedDofF=range(fluid_ndof)

    toc = time.process_time()
    logger.info("time to compute fluid matrices: {}".format(toc-tic))

    ##################################################################
    # Compute enrichment: Heaviside + Edge
    ##################################################################
    tic = time.process_time()
    (
        NegativeLSelements,
        PositiveLSelements,
        NegativeLStgtElements,
        PositiveLStgtElements,
    ) = objXFEM.getLocationElementsLS(fluid_elements, LevelSet, LevelSetTangent)

    EdgeEnrichedElementsInAllMesh = objXFEM.getEnrichedElements(
        fluidElements=fluid_elements, levelset=LevelSetTangent
    )
    IdElementTip = objXFEM.getElementContainingPoint(
        fluid_elements, fluid_nodes, [0.6, 0.65]
    )

    if (flag_write_gmsh_results == 1) and (rank == 0):
        msh2.mshWriter(
            cwd / (results_file + "_NegativeLSelements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[NegativeLSelements]},
        )
        msh2.mshWriter(
            cwd / (results_file + "_PositiveLSelements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[PositiveLSelements]},
        )
        msh2.mshWriter(
            cwd / (results_file + "_NegativeLStgtElements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[NegativeLStgtElements]},
        )
        msh2.mshWriter(
            cwd / (results_file + "_PositiveLStgtElements.msh"),
            fluid_nodes,
            {"type": "TRI3", "connectivity": fluid_elements[PositiveLStgtElements]},
        )



    IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = objXFEM.getMatrices(
        fluid_nodes,
        fluid_elements,
        LevelSet,
        LevelSetTangent,
        [celerity, rho],
        flag_edge_enrichment,
    )

    KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
    MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

    toc = time.process_time()
    logger.info("time to compute Heaviside enrichment: {}".format(toc-tic))

    #Enrichednodes = np.unique(fluid_elements[np.hstack(([HeavisideEnrichedElements,EdgeEnrichedElements]))])
    #Enrichednodes = np.unique(fluid_elements[np.hstack(([EnrichedElements,PositiveLStgtElements,EdgeEnrichedElementsInAllMesh]))])
    #Enrichednodes = np.unique(fluid_elements[np.hstack(([EnrichedElements,PositiveLStgtElements]))])
    #Enrichednodes = np.unique(fluid_elements[np.hstack(([NegativeLStgtElements]))])
    Enrichednodes = np.unique(fluid_elements[EnrichedElements])
    #Enrichednodes = np.unique(fluid_elements)
    SolvedDofA=Enrichednodes-1

    msh2.mshWriter(
        cwd / (results_file + "_EnrichedElements.msh"),
        fluid_nodes,
        {"type": "TRI3", "connectivity": fluid_elements[EnrichedElements.flatten()]},
    )

    #################################################################
    # Construct the whole system
    ##################################################################

    K=scipy.sparse.bmat( [[fluid_damping*KFF[SolvedDofF,:][:,SolvedDofF],fluid_damping*KAF[SolvedDofF,:][:,SolvedDofA]],
                                    [fluid_damping*KAF[SolvedDofA,:][:,SolvedDofF],fluid_damping*KAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )



    M=scipy.sparse.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                    [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )

    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################

    #gradient wrt Xc
    IIf_X, JJf_X, Vfak_gradient_X, Vfam_gradient_X, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_X,
        material=[celerity, rho],
    )
    
    dMFA_dtheta_X = scipy.sparse.csc_matrix( (Vfam_gradient_X,(IIf_X,JJf_X)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_X = scipy.sparse.csc_matrix( (Vfak_gradient_X,(IIf_X,JJf_X)), shape=(fluid_ndof,fluid_ndof) )

    #gradient wrt Yc
    IIf_Y, JJf_Y, Vfak_gradient_Y, Vfam_gradient_Y, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_Y,
        material=[celerity, rho],
    )

    
    dMFA_dtheta_Y = scipy.sparse.csc_matrix( (Vfam_gradient_Y,(IIf_Y,JJf_Y)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_Y = scipy.sparse.csc_matrix( (Vfak_gradient_Y,(IIf_Y,JJf_Y)), shape=(fluid_ndof,fluid_ndof) )

    #gradient wrt R
    IIf_R, JJf_R, Vfak_gradient_R, Vfam_gradient_R, _ = objXFEM.getGradientMatrices(
        fluid_nodes,
        fluid_elements[EnrichedElements],
        levelset=LevelSet,
        levelsetTg=None,
        levelsetGradient=LevelSet_gradient_R,
        material=[celerity, rho],
    )
    
    dMFA_dtheta_R = scipy.sparse.csc_matrix( (Vfam_gradient_R,(IIf_R,JJf_R)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_R = scipy.sparse.csc_matrix( (Vfak_gradient_R,(IIf_R,JJf_R)), shape=(fluid_ndof,fluid_ndof) )

    with open(cwd / (results_file + "_KAF_MAF.pck"), "wb") as f:
        pickle.dump(
            [
                dKFA_dtheta_X[SolvedDofF, :][:, SolvedDofA],
                dMFA_dtheta_X[SolvedDofF, :][:, SolvedDofA],
                KAF[SolvedDofF, :][:, SolvedDofA],
                MAF[SolvedDofF, :][:, SolvedDofA],
            ],
            f,
        )
    

    dK_X=scipy.sparse.bmat( [[None,fluid_damping*dKFA_dtheta_X[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_X[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_X=scipy.sparse.bmat( [[None,dMFA_dtheta_X[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_X[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )
    
    dK_Y=scipy.sparse.bmat( [[None,fluid_damping*dKFA_dtheta_Y[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_Y[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_Y=scipy.sparse.bmat( [[None,dMFA_dtheta_Y[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_Y[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )       

    dK_R=scipy.sparse.bmat( [[None,fluid_damping*dKFA_dtheta_R[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_R[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_R=scipy.sparse.bmat( [[None,dMFA_dtheta_R[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_R[SolvedDofA,:][:,SolvedDofF],None]
                                    ] ) 
    

    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis=1
    FF = np.zeros(fluid_ndof)
    frequencies=[]
    frf=[]
    frfgradient_X=[]
    frfgradient_Y=[]
    frfgradient_R=[]

    if (Flag_frf_analysis==1):
        logger.info("time at the beginning of the FRF: {}".format(time.ctime()))

        press_save=[]
        dpress_save_X=[]
        dpress_save_Y=[]
        dpress_save_R=[]
        disp_save=[]

        #extract frequencies for the associated processors
        freqCompute=listFreqPerProc[:,rank]
        freqCompute=freqCompute[freqCompute>0]

        for freq in freqCompute:

            #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega=2*np.pi*freq
            logger.info("proc number",rank,"frequency=",freq)

            FF[SolvedDofF]=-(KFF[SolvedDofF,:][:,IdnodeS2-1]-(omega**2)*MFF[SolvedDofF,:][:,IdnodeS2-1])*(np.ones((len(IdnodeS2))))
            FA = np.zeros(fluid_ndof)
            F  = FF[SolvedDofF]
            F  = np.concatenate((F,FA[SolvedDofA]))
            #F  = scipy.sparse.csc_matrix(F)

            #sol = mumps.spsolve(scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float'), F, comm=mycomm )
            sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=complex), F)
            
            #eigen_values_small,eigen_vectors= scipy.sparse.linalg.eigsh(K-(omega**2)*M,1,sigma=0,which='SM')
            #eigen_values_large,eigen_vectors= scipy.sparse.linalg.eigsh(K-(omega**2)*M,1,sigma=0,which='LM')


            #eigen_values_small,eigen_vectors= scipy.linalg.eig((K-(omega**2)*M).todense(),1,sigma=0,which='SM')
            #eigen_values_large,eigen_vectors= scipy.linalg.eig(K-(omega**2)*M,1,sigma=0,which='LM')

            
            #stop
            tmp_X=-(dK_X-(omega**2)*dM_X)*sol
            tmp_Y=-(dK_Y-(omega**2)*dM_Y)*sol
            tmp_R=-(dK_R-(omega**2)*dM_R)*sol

            #Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
            Dsol_Dtheta_X = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=complex)  , tmp_X )
            Dsol_Dtheta_Y = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=complex)  , tmp_Y )
            Dsol_Dtheta_R = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype=complex)  , tmp_R )

            press = np.zeros(fluid_ndof,dtype=complex)
            press[IdnodeS2-1] = np.ones(len(IdnodeS2))
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            enrichment=np.zeros(fluid_nnodes,dtype=complex)
            enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            CorrectedPressure=press
            #logger.info(SolvedDofA)
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*np.sign(LevelSet[SolvedDofA])
            #frf.append(silex_acou_lib_tri3.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
            #frf.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,press+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            frf.append(
                objXFEM.getQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    press,
                    enrichment,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )
            #frf[i]=xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent)
            #press_save.append(CorrectedPressure.copy())
            press_save.append(press.copy())

            Dpress_Dtheta = np.zeros([fluid_ndof,3],dtype=complex)
            Dpress_Dtheta[SolvedDofF,0] = Dsol_Dtheta_X[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF,1] = Dsol_Dtheta_Y[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF,2] = Dsol_Dtheta_R[list(range(len(SolvedDofF)))]
            Denrichment_Dtheta = np.zeros([fluid_ndof,3],dtype=complex)
            Denrichment_Dtheta[SolvedDofA,0]= Dsol_Dtheta_X[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            Denrichment_Dtheta[SolvedDofA,1]= Dsol_Dtheta_Y[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            Denrichment_Dtheta[SolvedDofA,2]= Dsol_Dtheta_R[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            #DCorrectedPressure_Dtheta=np.array(Dpress_Dtheta)
            #DCorrectedPressure_Dtheta[SolvedDofA]=DCorrectedPressure_Dtheta[SolvedDofA].T+np.array(Denrichment_Dtheta[SolvedDofA]*np.sign(LevelSet[SolvedDofA]).T)
            frfgradient_X.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    press,
                    enrichment*0.0,
                    Dpress_Dtheta[:,0],
                    Denrichment_Dtheta[:,0]*0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
            )                
            frfgradient_Y.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    press,
                    enrichment*0.0,
                    Dpress_Dtheta[:,1],
                    Denrichment_Dtheta[:,1]*0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
                )
            frfgradient_R.append(
                objXFEM.getGradientQuadraticPressure(
                    fluid_nodes,
                    fluid_elements5,
                    press,
                    enrichment*0.0,
                    Dpress_Dtheta[:,2],
                    Denrichment_Dtheta[:,2]*0.0,
                    LevelSet,
                    LevelSetTangent,
                    flag_edge_enrichment,
                )
                )
            dpress_save_X.append(Dpress_Dtheta[:,0].copy())
            dpress_save_Y.append(Dpress_Dtheta[:,1].copy())
            dpress_save_R.append(Dpress_Dtheta[:,2].copy())
        

        logger.info("time at the end of the FRF: {}".format(time.ctime()))
        frfsave=[frequencies,frf,frfgradient_X,frfgradient_Y,frfgradient_R]
        if rank!=0 :
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results==1) and (rank==0):
            msh2.mshWriter(
                cwd / (results_file + "_results_fluid_frf.msh"),
                fluid_nodes,
                [{"type": "TRI3", "connectivity": fluid_elements}],
                fields=[
                    {
                        "data": press_save,
                        "type": "nodal",
                        "dim": 1,
                        "nbsteps": len(press_save),
                        "name": "pressure",
                    },
                    {
                        "data": dpress_save_X,
                        "type": "nodal",
                        "dim": 1,
                        "nbsteps": len(press_save),
                        "name": "pressure gradient /X",
                    },
                    {
                        "data": dpress_save_Y,
                        "type": "nodal",
                        "dim": 1,
                        "nbsteps": len(press_save),
                        "name": "pressure gradient /Y",
                    },
                    {
                        "data": dpress_save_R,
                        "type": "nodal",
                        "dim": 1,
                        "nbsteps": len(press_save),
                        "name": "pressure gradient /R",
                    },
                ],
                append=True,
            )

        # save the FRF problem
        Allfrequencies=np.zeros(nbStep)
        Allfrf=np.zeros(nbStep)
        Allfrfgradient=np.zeros([nbStep,3])
        k=0
        if rank==0:
            for i in range(nproc):
                if i==0:
                    data=frfsave
                    logger.info(data)
                else:
                    logger.info(i)
                    data=comm.recv(source=i, tag=11)
                    #data=data_buffer

                for j in range(len(data[0])):
                    Allfrequencies[k]=data[0][j]
                    Allfrf[k]=data[1][j]
                    Allfrfgradient[k,0]=data[2][j]
                    Allfrfgradient[k,1]=data[3][j]
                    Allfrfgradient[k,2]=data[4][j]
                    k=k+1

            Allfrequencies, Allfrf,Allfrfgradient_X,Allfrfgradient_Y,Allfrfgradient_R = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient[:,0],Allfrfgradient[:,1],Allfrfgradient[:,2])))
            Allfrfsave=[np.array(list(Allfrequencies)),np.array(list(Allfrf)),np.array(list(Allfrfgradient_X)),np.array(list(Allfrfgradient_Y)),np.array(list(Allfrfgradient_R))]#,press_save,dpress_save]
            with open(cwd / (results_file + "_results.frf"), "wb") as f:
                pickle.dump(Allfrfsave, f)
            # save on mat file
            scipy.io.savemat(
                cwd / (results_file + "_results.mat"), mdict={"AllFRF": Allfrfsave}
            )
            scipy.io.savemat(
                cwd / (results_file_ini + "results.mat"), mdict={"AllFRF": Allfrfsave}
            )


#function for dealing with options
def manageOpt(argv,dV):
    #load default values
    freqMin     = dV.freqMin
    freqMax     = dV.freqMax
    nbStep      = dV.nbStep
    posStructX  = dV.posStructX
    posStructY  = dV.posStructY
    radiusStruct= dV.radiusStruct   
    
    #load info from MPI
    nbProc,rank,comm=mpiInfo()
    #load options
    opts,args = getopt.getopt(argv,"p:s:F:f:hp")
    for opt,arg in opts:
        if opt == "-s":
            nbStep  = int(arg)
        elif opt == "-F":
            freqMax = float(arg)
        elif opt == "-f":
            freqMin = float(arg)
        elif opt == "-p":
            tmp = np.array(arg.split(','),dtype=scipy.float32)
            posStructX=tmp[0]   
            if len(tmp)>1:
                posStructY=tmp[1]
            if len(tmp)>2:
                radiusStruct=tmp[2]
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    print ("Number of processors: ",nbProc)
    print ("Number of frequency steps: ",nbStep)
    print ("Maximum frequency: ",freqMax)
    print ("Minimum frequency: ",freqMin)
    print ("Structure position X: ",posStructX)
    print ("Structure position Y: ",posStructY)
    print ("Structure radius: ",radiusStruct)
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,posStructX,posStructY,radiusStruct)

#usage definition
def usage():
    dV=defaultV
    logger.info("Usage: ",sys.argv[0],"-psFfh [+arg]")
    logger.info("\t -p : number of processors (default value ",dV.nbProc,")")
    logger.info("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    logger.info("\t -F : maximum frequency (default value ",dV.freqMax,")")
    logger.info("\t -f : minimum frequency (default value ",dV.freqMin,")")

#default values
class defaultV:
    freqMin     = 35.0
    freqMax     = 80.0
    nbStep      = 22
    posStructX   = 1.5 
    posStructY   = 1
    radiusStruct   = 0.75

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)
