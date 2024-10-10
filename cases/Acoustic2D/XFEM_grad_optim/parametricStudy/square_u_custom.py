import string
import time
import numpy
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io
import getopt

import pylab as pl
import pickle

#import mumps

import sys
sys.path.append('../../../librairies')
import silex_lib_gmsh
import silex_lib_tri3_acou



import mumps

from mpi4py import MPI
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
    listFreq=scipy.zeros((nbFreqProc+varCase,nbProc))
    listAllFreq=scipy.linspace(freqInit,freqEnd,nbStep)
    #print(scipy.linspace(freqInit,freqEnd,nbStep))
    #build array of frequencies
    itF=0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP<nbStep:
                listFreq[itC,itP]=listAllFreq[itF]
                itF += 1

    #print(listFreq)
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

def RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,positionStructX,positionStructY,radiusStruct,angleStruct):
    # parallepipedic cavity with plane structure
    mesh_file='geom/xfem3_square'
    results_file_ini='results_square/xfem_3_'

    listFreqPerProc = computeFreqPerProc(nbStep,nbProc,freqMin,freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity=340.0
    rho=1.2
    fluid_damping=(1+0.01j)

    nproc=comm.Get_size()
    rank = comm.Get_rank()

    flag_write_gmsh_results=1

    flag_edge_enrichment=0

    h_struc=2.0

    ValpStructX=positionStructX*10000.
    ValpStructY=positionStructY*10000.
    ValpStructR=radiusStruct*10000.
    ValpStructT=angleStruct*10000.

    file_extension=str(ValpStructX)[0:5]+'_'+str(ValpStructY)[0:5]+'_'+str(ValpStructR)[0:4]+'_'+str(ValpStructT)[0:4]
    results_file=results_file_ini+file_extension
    print(results_file)

    # center and radius of half-circle
    x_pos_struc=positionStructX
    y_pos_struc=positionStructY
    radius_hcircle=radiusStruct
    print(x_pos_struc)
    print(y_pos_struc)
    print(radius_hcircle)
    print(angleStruct)
    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_fluid.msh',2)
    fluid_elements,Idnodes = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',2,1)
    fluid_elements5,Idnodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',2,5)

    fluid_nnodes   = fluid_nodes.shape[0]
    fluid_nelem    = fluid_elements.shape[0]
    fluid_ndof     = fluid_nnodes

    fluid_elements_boun,IdnodeS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',1,3)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,2)
        silex_lib_gmsh.WriteResults(results_file+'_control_volume_fluid_mesh',fluid_nodes,fluid_elements5,2)
        silex_lib_gmsh.WriteResults(results_file+'_fluid_boundary',fluid_nodes,fluid_elements_boun,1)
    print("nnodes for fluid=",fluid_nnodes)
    print("nelem for fluid=",fluid_nelem)



    ##################################################################
    # compute level set
    ##################################################################

    tic = time.clock()

    # create coordinates of nodes of the structure (half circle)
    nbNodesHC=50
    nbNodesSC=25
    nbNodesWall=1
    #
    angleU=angleStruct*scipy.pi/180
    thicknessU=1.
    # points
    xc2=x_pos_struc-radius_hcircle*scipy.cos(angleU)
    yc2=y_pos_struc-radius_hcircle*scipy.sin(angleU)
    xc4=x_pos_struc+radius_hcircle*scipy.cos(angleU)
    yc4=y_pos_struc+radius_hcircle*scipy.sin(angleU)
    #parameter for circles
    thetaHC=scipy.linspace(-angleU-scipy.pi/2,-angleU+scipy.pi/2,nbNodesHC)
    thetaSC=scipy.linspace(-angleU-scipy.pi/2,-angleU+scipy.pi/2,nbNodesSC)
    #inner large circle
    xNodesIHC=x_pos_struc-(radius_hcircle-thicknessU/2.)*scipy.sin(thetaHC)
    yNodesIHC=y_pos_struc-(radius_hcircle-thicknessU/2.)*scipy.cos(thetaHC)
    #outer larger circle
    xNodesOHC=x_pos_struc-(radius_hcircle+thicknessU/2.)*scipy.sin(thetaHC[::-1])
    yNodesOHC=y_pos_struc-(radius_hcircle+thicknessU/2.)*scipy.cos(thetaHC[::-1])
    #small circle 2
    xNodesC2=xc2-thicknessU/2*scipy.sin(thetaSC[::-1]+scipy.pi)
    yNodesC2=yc2-thicknessU/2*scipy.cos(thetaSC[::-1]+scipy.pi)
    #small circle 4
    xNodesC4=xc4+thicknessU/2*scipy.sin(thetaSC[::-1])
    yNodesC4=yc4+thicknessU/2*scipy.cos(thetaSC[::-1])

    struc_nodes=scipy.vstack([scipy.hstack([xNodesIHC,xNodesC2[1:],xNodesOHC[1:],xNodesC4[1:]]),scipy.hstack([yNodesIHC,yNodesC2[1:],yNodesOHC[1:],yNodesC4[1:]])]).transpose()
    print(struc_nodes)

    nbNodesAllStruct=2*nbNodesHC+2*nbNodesSC-4
    lCA=scipy.linspace(1,nbNodesAllStruct-1,nbNodesAllStruct-1)
    lCB=scipy.linspace(2,nbNodesAllStruct,nbNodesAllStruct-1)

    struc_elements=scipy.vstack([lCA,lCB]).transpose()
    struc_elements=scipy.vstack([struc_elements,[1,nbNodesAllStruct]])
    print(struc_elements)


    LevelSet=scipy.zeros(fluid_nnodes)
    LevelSet_gradient_X=scipy.zeros(fluid_nnodes)
    LevelSet_gradient_Y=scipy.zeros(fluid_nnodes)
    LevelSet_gradient_R=scipy.zeros(fluid_nnodes)
    LevelSet_gradient_T=scipy.zeros(fluid_nnodes)
    IndicZone=scipy.zeros(fluid_nnodes)
    # level set 3 cases
    for itN in range(fluid_nnodes):
        xNcurr=fluid_nodes[itN,0]
        yNcurr=fluid_nodes[itN,1]
        #check areas
        if angleU != scipy.pi/2 or angleU != 3*scipy.pi/2:
            zoneA=yNcurr-scipy.tan(angleU)*(xNcurr-x_pos_struc)-y_pos_struc<0
        else:
            zoneA=xNcurr-x_pos_struc<0

        zoneAA=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle>=0
        zoneAB=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle<0
        if angleU != scipy.pi and angleU != 0. and angleU != 2*scipy.pi:
            zoneB=yNcurr-scipy.tan(angleU+scipy.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc>0
            zoneC=yNcurr-scipy.tan(angleU+scipy.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc<0
        else:
            zoneB=xNcurr-x_pos_struc>0
            zoneC=xNcurr-x_pos_struc<0

        if zoneA:
            if zoneAA:
                radiusOC=radius_hcircle+thicknessU/2.
                LevelSet[itN]=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusOC
                LevelSet_gradient_X[itN]=-(xNcurr-x_pos_struc)/(LevelSet[itN]+radiusOC)
                LevelSet_gradient_Y[itN]=-(yNcurr-y_pos_struc)/(LevelSet[itN]+radiusOC)
                LevelSet_gradient_R[itN]=-1 
                LevelSet_gradient_T[itN]=0 
                IndicZone[itN]=1
            elif zoneAB: 
                radiusIC=radius_hcircle-thicknessU/2.
                LevelSet[itN]=-(scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusIC)
                LevelSet_gradient_X[itN]=(xNcurr-x_pos_struc)/(-LevelSet[itN]+radiusIC)
                LevelSet_gradient_Y[itN]=(yNcurr-y_pos_struc)/(-LevelSet[itN]+radiusIC)
                LevelSet_gradient_R[itN]=1 
                LevelSet_gradient_T[itN]=0
                IndicZone[itN]=2
        elif zoneC:
            radiusC2=thicknessU/2.
            LevelSet[itN]=scipy.sqrt((xNcurr-xc2)**2+(yNcurr-yc2)**2)-radiusC2
            LevelSet_gradient_X[itN]=-(xNcurr-xc2)/(LevelSet[itN]+radiusC2)
            LevelSet_gradient_Y[itN]=-(yNcurr-yc2)/(LevelSet[itN]+radiusC2)
            LevelSet_gradient_R[itN]=0
            LevelSet_gradient_T[itN]=-radius_hcircle*scipy.sin(angleU)*LevelSet_gradient_X[itN]+radius_hcircle*scipy.cos(angleU)*LevelSet_gradient_Y[itN]
            IndicZone[itN]=3
        elif zoneB:
            radiusC4=thicknessU/2.
            LevelSet[itN]=scipy.sqrt((xNcurr-xc4)**2+(yNcurr-yc4)**2)-radiusC4
            LevelSet_gradient_X[itN]=-(xNcurr-xc4)/(LevelSet[itN]+radiusC4)
            LevelSet_gradient_Y[itN]=-(yNcurr-yc4)/(LevelSet[itN]+radiusC4)
            LevelSet_gradient_R[itN]=0
            LevelSet_gradient_T[itN]=radius_hcircle*scipy.sin(angleU)*LevelSet_gradient_X[itN]-radius_hcircle*scipy.cos(angleU)*LevelSet_gradient_Y[itN]
            IndicZone[itN]=4
        
    LevelSetTangent=fluid_nodes[:,1]-max(fluid_nodes[:,1])

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,2,[[LevelSetTangent,'nodal',1,'Tangent level set']])
        silex_lib_gmsh.WriteResults(results_file+'_level_set',fluid_nodes,fluid_elements,2,[[LevelSet,'nodal',1,'Level set'],[LevelSet_gradient_X,'nodal',1,'Level set Grad X'],[LevelSet_gradient_Y,'nodal',1,'Level set Grad Y'],[LevelSet_gradient_R,'nodal',1,'Level set Grad R'],[LevelSet_gradient_T,'nodal',1,'Level set Grad T'],[IndicZone,'nodal',1,'Indic Zones']])

    toc = time.clock()
    print("time to compute level set:",toc-tic)


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.clock()

    
    struc_boun=scipy.array([3])

    silex_lib_gmsh.WriteResults(results_file+'_struc_mesh',struc_nodes,struc_elements,1)

    EnrichedElements,NbEnrichedElements=silex_lib_tri3_acou.getenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements)
    EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])-1

    toc = time.clock()
    print("time to find surface enriched elements:",toc-tic)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],2)

    tic = time.clock()

    EdgeEnrichedElements,nbenrelts = silex_lib_tri3_acou.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes,fluid_elements)
    EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1

    toc = time.clock()
    print("time to find edge enriched elements:",toc-tic)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_edge_enriched_elements',fluid_nodes,fluid_elements[EdgeEnrichedElements],2)

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################
    tic = time.clock()

    IIf,JJf,Vffk,Vffm=silex_lib_tri3_acou.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

    KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    SolvedDofF=scipy.setdiff1d(list(range(fluid_ndof)),IdnodeS2-1)

    toc = time.clock()
    print("time to compute fluid matrices:",toc-tic)

    ##################################################################
    # Compute enrichment: Heaviside + Edge
    ##################################################################
    tic = time.clock()

    IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_tri3_acou.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,LevelSetTangent,celerity,rho,flag_edge_enrichment)

    KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
    MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

    toc = time.clock()
    print("time to compute Heaviside enrichment:",toc-tic)

    Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

    SolvedDofA=Enrichednodes-1

    silex_lib_gmsh.WriteResults(results_file+'_EnrichedElements',fluid_nodes,fluid_elements[scipy.hstack(([EnrichedElements]))],2)

    #################################################################
    # Construct the whole system
    ##################################################################

    K=scipy.sparse.construct.bmat( [[fluid_damping*KFF[SolvedDofF,:][:,SolvedDofF],fluid_damping*KAF[SolvedDofF,:][:,SolvedDofA]],
                                    [fluid_damping*KAF[SolvedDofA,:][:,SolvedDofF],fluid_damping*KAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )



    M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                    [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )

    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################

    print(silex_lib_tri3_acou.globalacousticgradientmatrices.__doc__)

    #gradient wrt Xc
    IIf_X,JJf_X,Vfak_gradient_X,Vfam_gradient_X=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient_X,LevelSet)
    dMFA_dtheta_X = scipy.sparse.csc_matrix( (Vfam_gradient_X,(IIf_X,JJf_X)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_X = scipy.sparse.csc_matrix( (Vfak_gradient_X,(IIf_X,JJf_X)), shape=(fluid_ndof,fluid_ndof) )

    #gradient wrt Yc
    IIf_Y,JJf_Y,Vfak_gradient_Y,Vfam_gradient_Y=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient_Y,LevelSet)
    dMFA_dtheta_Y = scipy.sparse.csc_matrix( (Vfam_gradient_Y,(IIf_Y,JJf_Y)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_Y = scipy.sparse.csc_matrix( (Vfak_gradient_Y,(IIf_Y,JJf_Y)), shape=(fluid_ndof,fluid_ndof) )

    #gradient wrt R
    IIf_R,JJf_R,Vfak_gradient_R,Vfam_gradient_R=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient_R,LevelSet)
    dMFA_dtheta_R = scipy.sparse.csc_matrix( (Vfam_gradient_R,(IIf_R,JJf_R)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_R = scipy.sparse.csc_matrix( (Vfak_gradient_R,(IIf_R,JJf_R)), shape=(fluid_ndof,fluid_ndof) )
    
    #gradient wrt R
    IIf_T,JJf_T,Vfak_gradient_T,Vfam_gradient_T=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient_T,LevelSet)
    dMFA_dtheta_T = scipy.sparse.csc_matrix( (Vfam_gradient_T,(IIf_T,JJf_T)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta_T = scipy.sparse.csc_matrix( (Vfak_gradient_T,(IIf_T,JJf_T)), shape=(fluid_ndof,fluid_ndof) )

    
    f=open(results_file+'_KAF_MAF.pck','wb')
    pickle.dump([dKFA_dtheta_X[SolvedDofF,:][:,SolvedDofA],dMFA_dtheta_X[SolvedDofF,:][:,SolvedDofA],KAF[SolvedDofF,:][:,SolvedDofA],MAF[SolvedDofF,:][:,SolvedDofA]], f)
    f.close()

    dK_X=scipy.sparse.construct.bmat( [[None,fluid_damping*dKFA_dtheta_X[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_X[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_X=scipy.sparse.construct.bmat( [[None,dMFA_dtheta_X[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_X[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )
    
    dK_Y=scipy.sparse.construct.bmat( [[None,fluid_damping*dKFA_dtheta_Y[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_Y[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_Y=scipy.sparse.construct.bmat( [[None,dMFA_dtheta_Y[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_Y[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )       

    dK_R=scipy.sparse.construct.bmat( [[None,fluid_damping*dKFA_dtheta_R[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_R[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_R=scipy.sparse.construct.bmat( [[None,dMFA_dtheta_R[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_R[SolvedDofA,:][:,SolvedDofF],None]
                                    ] ) 

    dK_T=scipy.sparse.construct.bmat( [[None,fluid_damping*dKFA_dtheta_T[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta_T[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM_T=scipy.sparse.construct.bmat( [[None,dMFA_dtheta_T[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta_T[SolvedDofA,:][:,SolvedDofF],None]
                                    ] ) 
    

    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis=1
    FF = scipy.zeros(fluid_ndof)
    frequencies=[]
    frf=[]
    frfgradient_X=[]
    frfgradient_Y=[]
    frfgradient_R=[]
    frfgradient_T=[]

    if (Flag_frf_analysis==1):
        print("time at the beginning of the FRF:",time.ctime())

        press_save=[]
        dpress_save_X=[]
        dpress_save_Y=[]
        dpress_save_R=[]
        dpress_save_T=[]
        disp_save=[]

        #extract frequencies for the associated processors
        freqCompute=listFreqPerProc[:,rank]
        freqCompute=freqCompute[freqCompute>0]

        for freq in freqCompute:

            #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega=2*scipy.pi*freq
            print("proc number",rank,"frequency=",freq)

            FF[SolvedDofF]=-(KFF[SolvedDofF,:][:,IdnodeS2-1]-(omega**2)*MFF[SolvedDofF,:][:,IdnodeS2-1])*(scipy.ones((len(IdnodeS2))))
            FA = scipy.zeros(fluid_ndof)
            F  = FF[SolvedDofF]
            F  = scipy.append(F,FA[SolvedDofA])
            #F  = scipy.sparse.csc_matrix(F)

            ## 
            pbFreq=K-(omega**2)*M

            sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(pbFreq,dtype=complex), F)
       
           
            #stop
            tmp_X=-(dK_X-(omega**2)*dM_X)*sol
            tmp_Y=-(dK_Y-(omega**2)*dM_Y)*sol
            tmp_R=-(dK_R-(omega**2)*dM_R)*sol
            tmp_T=-(dK_T-(omega**2)*dM_T)*sol
           
            #Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
            Dsol_Dtheta_X = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(pbFreq,dtype=complex)  , tmp_X )
            Dsol_Dtheta_Y = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(pbFreq,dtype=complex)  , tmp_Y )
            Dsol_Dtheta_R = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(pbFreq,dtype=complex)  , tmp_R )
            Dsol_Dtheta_T = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(pbFreq,dtype=complex)  , tmp_T )

            press = scipy.zeros(fluid_ndof,dtype=complex)
            press[IdnodeS2-1] = scipy.ones(len(IdnodeS2))
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            enrichment=scipy.zeros(fluid_nnodes,dtype=complex)
            enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            CorrectedPressure=press
            #print(SolvedDofA)
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
            #####################
            #####################
            frf.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            #####################
            #####################
            press_save.append(CorrectedPressure.copy())
            #####################
            #####################
            Dpress_Dtheta = scipy.zeros([fluid_ndof,4],dtype=complex)
            Dpress_Dtheta[SolvedDofF,0] = Dsol_Dtheta_X[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF,1] = Dsol_Dtheta_Y[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF,2] = Dsol_Dtheta_R[list(range(len(SolvedDofF)))]
            Dpress_Dtheta[SolvedDofF,3] = Dsol_Dtheta_T[list(range(len(SolvedDofF)))]
            #####################
            #####################
            Denrichment_Dtheta = scipy.zeros([fluid_ndof,4],dtype=complex)
            Denrichment_Dtheta[SolvedDofA,0]= Dsol_Dtheta_X[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            Denrichment_Dtheta[SolvedDofA,1]= Dsol_Dtheta_Y[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            Denrichment_Dtheta[SolvedDofA,2]= Dsol_Dtheta_R[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            Denrichment_Dtheta[SolvedDofA,3]= Dsol_Dtheta_T[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            #####################
            #####################
            DCorrectedPressure_Dtheta=scipy.array(Dpress_Dtheta)
            DCorrectedPressure_Dtheta[SolvedDofA,0]=DCorrectedPressure_Dtheta[SolvedDofA,0].T+scipy.array(Denrichment_Dtheta[SolvedDofA,0]*scipy.sign(LevelSet[SolvedDofA]).T)
            DCorrectedPressure_Dtheta[SolvedDofA,1]=DCorrectedPressure_Dtheta[SolvedDofA,1].T+scipy.array(Denrichment_Dtheta[SolvedDofA,1]*scipy.sign(LevelSet[SolvedDofA]).T)
            DCorrectedPressure_Dtheta[SolvedDofA,2]=DCorrectedPressure_Dtheta[SolvedDofA,2].T+scipy.array(Denrichment_Dtheta[SolvedDofA,2]*scipy.sign(LevelSet[SolvedDofA]).T)
            DCorrectedPressure_Dtheta[SolvedDofA,3]=DCorrectedPressure_Dtheta[SolvedDofA,3].T+scipy.array(Denrichment_Dtheta[SolvedDofA,3]*scipy.sign(LevelSet[SolvedDofA]).T)
            #####################
            #####################
            frfgradient_X.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,CorrectedPressure+0j,DCorrectedPressure_Dtheta[:,0]+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta[:,0]+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            frfgradient_Y.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,CorrectedPressure+0j,DCorrectedPressure_Dtheta[:,1]+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta[:,1]+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            frfgradient_R.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,CorrectedPressure+0j,DCorrectedPressure_Dtheta[:,2]+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta[:,2]+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            frfgradient_T.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,CorrectedPressure+0j,DCorrectedPressure_Dtheta[:,3]+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta[:,3]+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            #####################
            #####################
            dpress_save_X.append(DCorrectedPressure_Dtheta[:,0].copy())
            dpress_save_Y.append(DCorrectedPressure_Dtheta[:,1].copy())
            dpress_save_R.append(DCorrectedPressure_Dtheta[:,2].copy())
            dpress_save_T.append(DCorrectedPressure_Dtheta[:,3].copy())
        
        #####################
        #####################
        print("time at the end of the FRF:",time.ctime())
        frfsave=[frequencies,frf,frfgradient_X,frfgradient_Y,frfgradient_R,frfgradient_T]
        if rank!=0 :
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results==1) and (rank==0):
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,2,[[numpy.real(press_save),'nodal',1,'pressure (real)'],[numpy.imag(press_save),'nodal',1,'pressure (imaginary)'],[numpy.absolute(press_save),'nodal',1,'pressure (norm)'],[numpy.real(dpress_save_X),'nodal',1,'pressure gradient X (real)'],[numpy.imag(dpress_save_X),'nodal',1,'pressure gradient X (imaginary)'],[numpy.absolute(dpress_save_X),'nodal',1,'pressure gradient X (norm)'],[numpy.real(dpress_save_Y),'nodal',1,'pressure gradient Y (real)'],[numpy.imag(dpress_save_Y),'nodal',1,'pressure gradient Y (imaginary)'],[numpy.absolute(dpress_save_Y),'nodal',1,'pressure gradient Y (norm)'],[numpy.real(dpress_save_R),'nodal',1,'pressure gradient R (real)'],[numpy.imag(dpress_save_R),'nodal',1,'pressure gradient R (imaginary)'],[numpy.absolute(dpress_save_R),'nodal',1,'pressure gradient R (norm)'],[numpy.real(dpress_save_T),'nodal',1,'pressure gradient T (real)'],[numpy.imag(dpress_save_T),'nodal',1,'pressure gradient T (imaginary)'],[numpy.absolute(dpress_save_T),'nodal',1,'pressure gradient T (norm)']])

        #####################
        #####################
        # save the FRF problem
        Allfrequencies=scipy.zeros(nbStep)
        Allfrf=scipy.zeros(nbStep)
        Allfrfgradient=scipy.zeros([nbStep,4])
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
                    Allfrfgradient[k,0]=data[2][j]
                    Allfrfgradient[k,1]=data[3][j]
                    Allfrfgradient[k,2]=data[4][j]
                    Allfrfgradient[k,3]=data[5][j]
                    k=k+1
            #####################
            #####################
            Allfrequencies, Allfrf,Allfrfgradient_X,Allfrfgradient_Y,Allfrfgradient_R,Allfrfgradient_T = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient[:,0],Allfrfgradient[:,1],Allfrfgradient[:,2],Allfrfgradient[:,3])))
            Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf)),scipy.array(list(Allfrfgradient_X)),scipy.array(list(Allfrfgradient_Y)),scipy.array(list(Allfrfgradient_R)),scipy.array(list(Allfrfgradient_T))]#,press_save,dpress_save]
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
    posStructX  = dV.posStructX
    posStructY  = dV.posStructY
    radiusStruct= dV.radiusStruct 
    angleStruct= dV.angleStruct 
    
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
            tmp = scipy.array(arg.split(','),dtype=scipy.float32)
            posStructX=tmp[0]   
            if len(tmp)>1:
                posStructY=tmp[1]
            if len(tmp)>2:
                radiusStruct=tmp[2]
            if len(tmp)>3:
                angleStruct=tmp[3]
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
    print ("Structure angle: ",angleStruct)
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,posStructX,posStructY,radiusStruct,angleStruct)

#usage definition
def usage():
    dV=defaultV
    print("Usage: ",sys.argv[0],"-psFfh [+arg]")
    print("\t -p : number of processors (default value ",dV.nbProc,")")
    print("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    print("\t -F : maximum frequency (default value ",dV.freqMax,")")
    print("\t -f : minimum frequency (default value ",dV.freqMin,")")

#default values
class defaultV:
    freqMin     = 35.0
    freqMax     = 80.0
    nbStep      = 22
    posStructX   = 1.5 
    posStructY   = 1
    radiusStruct   = 0.75
    angleStruct   = 0

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)
