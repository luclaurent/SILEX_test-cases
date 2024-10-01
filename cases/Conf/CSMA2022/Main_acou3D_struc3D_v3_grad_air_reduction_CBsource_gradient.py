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
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io

import pickle

import sys
import os
from shutil import copyfile
#sys.path.append('../../librairies')


import SILEX.silex_lib_xfem_acou_tet4 as silex_lib_xfem_acou_tet4
import SILEX.silex_lib_gmsh as silex_lib_gmsh
#import silex_lib_dkt_fortran as silex_lib_dkt

#import silex_lib_tet4_fortran_test as silex_lib_tet4

#import silex_lib_porous_tet4_fortran

# manage mumps librairie
import mumps
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#nproc = comm.Get_size()
#rank = comm.Get_rank()

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
    listFreq = np.zeros((nbFreqProc+varCase, nbProc))
    listAllFreq = np.linspace(freqInit, freqEnd, nbStep)
    # print(np.linspace(freqInit,freqEnd,nbStep))
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
    results_file_ini = 'results/cavity_acou3D_struc_3D_v3_air_reduction_CB_source_gradient'

    listFreqPerProc = computeFreqPerProc(nbStep, nbProc, freqMin, freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity = 340.0
    rho = 1.2
    fluid_damping = (1+0.01j)
    nb_mode_F = 1000

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
            file_extension = file_extension+'_'+"{:.{}E}".format(paraVal[i], 2)+'_NbModesFluid_'+str(nb_mode_F)

    results_file = results_file_ini+'_'+file_extension
    print(results_file)
    

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.process_time()

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

    fluid_ndof  = fluid_nnodes

    if rank == 0:
        print("Number of nodes:", fluid_nnodes)
        print("Number of elements in air:", fluid_nelem1)
        print("Number of nodes in air:", fluid_nnodes1)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults(
            results_file+'_air_cavity_Mesh1', fluid_nodes, fluid_elements1, 4)
        silex_lib_gmsh.WriteResults(
            results_file+'_air_controlled_volume_Mesh5', fluid_nodes, fluid_elements5, 4)

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
    #     print("nnodes for structure=", struc_nnodes)
    #     print("nelem for structure=", struc_nelem)

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
    print("Parameters values")
    print("Xc ",lx3," Yc ",ly3," Zc ",lz3," R ",R)
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
        print("time to compute level set:", toc-tic)

    if (flag_write_gmsh_results == 1) and (rank == 0):
        # silex_lib_gmsh.WriteResults2(
        #     results_file+'_struc_air_interface', struc_nodes, struc_elements, 2)
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
    tic = time.process_time()

    LSEnrichedElements, NbLSEnrichedElements = silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(
        fluid_elements1, LevelSet)
    LSEnrichedElements = LSEnrichedElements[list(range(NbLSEnrichedElements))]
    silex_lib_gmsh.WriteResults2(results_file+'_LS_enriched_elements',
                                fluid_nodes, fluid_elements1[LSEnrichedElements], 4)
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
        print("time to find enriched elements:", toc-tic)

    tic = time.process_time()

    if (flag_write_gmsh_results == 1) and (rank == 0):
        silex_lib_gmsh.WriteResults2(
            results_file+'_enriched_elements', fluid_nodes, fluid_elements1[EnrichedElements], 4)
        LS_moins_enriched = np.setdiff1d(LSEnrichedElements, EnrichedElements)
        enriched_moins_LS = np.setdiff1d(EnrichedElements, LSEnrichedElements)
        silex_lib_gmsh.WriteResults2(
            results_file+'_LS_moins_enriched', fluid_nodes, fluid_elements1[LS_moins_enriched], 4)
        silex_lib_gmsh.WriteResults2(
                results_file+'_enriched_moins_LS', fluid_nodes, fluid_elements1[enriched_moins_LS], 4)
    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################

    tic = time.process_time()

    IIf, JJf, Vffk, Vffm = silex_lib_xfem_acou_tet4.globalacousticmatrices(
        fluid_elements1, fluid_nodes, celerity, rho)

    KFF = scipy.sparse.csc_matrix(
        (Vffk, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))
    MFF = scipy.sparse.csc_matrix(
        (Vffm, (IIf, JJf)), shape=(fluid_ndof, fluid_ndof))

    SolvedDofF = list(range(fluid_ndof))
    SolvedDofB=np.hstack([9-1]) # 9 : node number where acoustic source is imposed
    SolvedDofI=np.setdiff1d(SolvedDofF,SolvedDofB)
    #SolvedDofI=SolvedDofF

    ##################################################################
    # Compute Heaviside enrichment
    ##################################################################
    tic = time.process_time()

    Enrichednodes = np.unique(fluid_elements1[EnrichedElements])

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

    toc = time.process_time()
    if rank == 0:
        print("time to compute Heaviside enrichment:", toc-tic)
        
    ##################################################################
    # Compute eigen modes of the fluid: internal dof I
    ##################################################################
    tic = time.process_time()

    eigen_values_I,eigen_vectors_I= scipy.sparse.linalg.eigsh(KFF[SolvedDofI,:][:,SolvedDofI],nb_mode_F,MFF[SolvedDofI,:][:,SolvedDofI],sigma=0,which='LM')

    freq_eigv_I=list(np.sqrt(eigen_values_I)/(2*scipy.pi))
    print(freq_eigv_I)
    eigen_vector_F_list=[]
    for i in range(nb_mode_F):
        tmp=np.zeros((fluid_ndof) , dtype='float')
        tmp[SolvedDofI]=eigen_vectors_I[:,i].real
        eigen_vector_F_list.append(tmp)

    if rank==0:
        print ("LAST fluid eigen frequencies : ",freq_eigv_I[-1])
     
    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes,fluid_elements1,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

    toc = time.process_time()
    if rank==0:
        print ("time for computing the fluid modes:",toc-tic)

    ##################################################################
    # Compute Psi_IA for the fluid: Psi_IA = - KII^{-1} * KIA
    ##################################################################
    tic = time.process_time()

    print ("Compute PSI_IA")
    #omega_cst=0.0*2.0*scipy.pi
    #MySolve = scipy.sparse.linalg.factorized( KFF[SolvedDofI,:][:,SolvedDofI]-(omega_cst**2)*MFF[SolvedDofI,:][:,SolvedDofI] ) # Makes LU decomposition.
    MySolve = scipy.sparse.linalg.factorized( KFF[SolvedDofI,:][:,SolvedDofI]) # Makes LU decomposition.
    print("LU decomposition has been made")
    Psi_IA=np.zeros((len(SolvedDofI),len(SolvedDofA)))

    j=0
    for One_dof in SolvedDofA:
        KIA_i_column=-KAF[SolvedDofI,:][:,One_dof]
        Xi=MySolve( KIA_i_column.todense() )
        Psi_IA[:,j]=np.array(Xi)[:,0]
        j=j+1

    Psi_IA=scipy.sparse.csc_matrix(Psi_IA)
    toc = time.process_time()
    print ("time to compute PSI_IA:",toc-tic)

    ##################################################################
    # Compute Psi_IB for the fluid: Psi_IB = - KII^{-1} * KIB
    ##################################################################
    tic = time.process_time()
    print ("Compute PSI_IB")
    Psi_IB=np.zeros((len(SolvedDofI),len(SolvedDofB)))

    j=0
    for One_dof in SolvedDofB:
        KIB_i_column=-KFF[SolvedDofI,:][:,One_dof]
        Xi=MySolve( KIB_i_column.todense() )
        Psi_IB[:,j]=np.array(Xi)[:,0]
        j=j+1


##    if rank==0:
##        print ("Compute PSI_IB")
##
##    if rank==0:
##        Psi_IB=np.zeros((len(SolvedDofI),len(SolvedDofB)))
##
##    i=0
##    j=0
##    while j<len(SolvedDofB):
##        j=i+rank
##        if j<len(SolvedDofB):
##            One_dof=[SolvedDofB[j]]
##            KIB_i_column=-KFF[SolvedDofI,:][:,One_dof]
##            Xi=MySolve( KIB_i_column.todense() )
##            if rank!=0:
##                comm.send([Xi,j], dest=0, tag=11)
##            if rank==0:
##                for k in range(nproc):
##                    if k!=0:
##                        if i+k<len(SolvedDofB):
##                            [Xi,j]=comm.recv(source=k, tag=11)
##                            Psi_IB[:,j]=np.array(Xi)[:,0]
##                    else:
##                        Psi_IB[:,j]=np.array(Xi)[:,0]
##            i=i+nproc

    if rank==0:
        print("End of PSI_IB computing, send to other proc.")

##    if rank==0:
##        for i in range(nproc):
##            if i!=0:
##                comm.send(Psi_IB, dest=i, tag=11)
##    if rank!=0:
##        Psi_IB=comm.recv(source=0, tag=11)

    Psi_IB=scipy.sparse.csc_matrix(Psi_IB)
    toc = time.process_time()

    if rank==0:
        print ("time to compute PSI_FB:",toc-tic)


    
    ##################################################################
    # Construct the whole system
    #################################################################

    # Fluid part

    VK_diag_mm = eigen_values_I
    VM_diag_mm = eigen_values_I/eigen_values_I
    IIDmm = list(range(nb_mode_F))
    JJDmm = list(range(nb_mode_F))

    K_diag_mm= scipy.sparse.csc_matrix( (VK_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
    M_diag_mm= scipy.sparse.csc_matrix( (VM_diag_mm,(IIDmm,JJDmm)), shape=(nb_mode_F,nb_mode_F) )
    Khat_AA = KAA[SolvedDofA,:][:,SolvedDofA]+Psi_IA.T*KAF[SolvedDofI,:][:,SolvedDofA]
    Khat_BB = KFF[SolvedDofB,:][:,SolvedDofB]+Psi_IB.T*KFF[SolvedDofI,:][:,SolvedDofB]

    Khat_BA = KFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA


    Mstar_IA = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IA+MAF[SolvedDofI,:][:,SolvedDofA]
    Mstar_IB = MFF[SolvedDofI,:][:,SolvedDofI]*Psi_IB+MFF[SolvedDofI,:][:,SolvedDofB]

    Mhat_AA = MAA[SolvedDofA,:][:,SolvedDofA]+scipy.sparse.csc_matrix((Psi_IA.T).todense()*Mstar_IA.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*Psi_IA
    Mhat_BB = MFF[SolvedDofB,:][:,SolvedDofB]+scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IB.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IB

    Mhat_mA = eigen_vectors_I.T*Mstar_IA
    Mhat_mB = eigen_vectors_I.T*Mstar_IB
    
    Mhat_BA = scipy.sparse.csc_matrix((Psi_IB.T).todense()*Mstar_IA.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*Psi_IA

    Kreduc=scipy.sparse.construct.bmat( [[fluid_damping*K_diag_mm,None,       None],
                                                       [None,     fluid_damping*Khat_BB,    fluid_damping*Khat_BA],
                                                       [None,     fluid_damping*Khat_BA.T,  fluid_damping*Khat_AA],
                                                       ]
                                                      )
        
    Mreduc=scipy.sparse.construct.bmat( [[M_diag_mm,    Mhat_mB,    Mhat_mA],
                                         [Mhat_mB.T,    Mhat_BB,    Mhat_BA],
                                         [Mhat_mA.T,    Mhat_BA.T,  Mhat_AA],
                                         ]
                                       )
     
    ##################################################################
    # Build Second member
    ##################################################################

    # To impose the load on the fluid:
    # fluid node number 9
    UF = np.zeros(2*fluid_ndof, dtype=float)
    UF[9-1] = 3.1250E-05

    #SolvedDof = np.hstack([SolvedDofF, SolvedDofA+fluid_ndof])

    UFreduc = np.zeros(fluid_ndof, dtype=float)
    UFreduc[9-1] = 3.1250E-05
    Freduced_F=np.hstack([np.zeros(nb_mode_F,dtype=float),UF[SolvedDofB]])

    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################
    #print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
    dKreduc=list()
    dMreduc=list()
    for itP in range(0,nbPara):
        print(' Build gradient matrices for parameter '+NamePara[itP])
        #
        IIf,JJf,Vfak_gradient,Vfam_gradient=\
            silex_lib_xfem_acou_tet4.globalacousticgradientmatrices(fluid_nodes,\
                fluid_elements1,LevelSet,celerity,rho,LevelSetGradient[itP])
        dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )


        tic = time.process_time()

        print ("Compute DPSI_IA_Dtheta")
        DPsi_IA_Dtheta=np.zeros((len(SolvedDofI),len(SolvedDofA)))

        j=0
        for One_dof in SolvedDofA:
            DKIA_Dtheta_i_column=-dKFA_dtheta[SolvedDofI,:][:,One_dof]
            Xi=MySolve( DKIA_Dtheta_i_column.todense() )
            DPsi_IA_Dtheta[:,j]=np.array(Xi)[:,0]
            j=j+1

        DPsi_IA_Dtheta=scipy.sparse.csc_matrix(DPsi_IA_Dtheta)
        toc = time.process_time()
        print ("time to compute DPSI_IA_Dtheta:",toc-tic)


        DKhat_AA_Dtheta = DPsi_IA_Dtheta.T*KAF[SolvedDofI,:][:,SolvedDofA]+Psi_IA.T*dKFA_dtheta[SolvedDofI,:][:,SolvedDofA]

        DKhat_BA_Dtheta = KFF[SolvedDofB,:][:,SolvedDofI]*DPsi_IA_Dtheta

        DMstar_IA_Dtheta= MFF[SolvedDofI,:][:,SolvedDofI]*DPsi_IA_Dtheta+dMFA_dtheta[SolvedDofI,:][:,SolvedDofA]

        DMhat_AA_Dtheta = scipy.sparse.csc_matrix((DPsi_IA_Dtheta.T).todense()*Mstar_IA.todense())+scipy.sparse.csc_matrix((Psi_IA.T).todense()*DMstar_IA_Dtheta.todense())+MAF[SolvedDofA,:][:,SolvedDofI]*DPsi_IA_Dtheta+dMFA_dtheta[SolvedDofA,:][:,SolvedDofI]*Psi_IA

        DMhat_mA_Dtheta = eigen_vectors_I.T*DMstar_IA_Dtheta

        DMhat_BA_Dtheta = scipy.sparse.csc_matrix((Psi_IB.T).todense()*DMstar_IA_Dtheta.todense())+MFF[SolvedDofB,:][:,SolvedDofI]*DPsi_IA_Dtheta

#build full stiffness and mass gradient matrices
##        dK.append(scipy.sparse.construct.bmat( [
##                    [None,fluid_damping*dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
##                    [fluid_damping*dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]] ))
##        dM.append(scipy.sparse.construct.bmat( [s
##                    [None,dMFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
##                    [dMFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]] ))
        DK_Dtheta=scipy.sparse.construct.bmat( [ [K_diag_mm*0.0,    None,               None],
                                                 [None,             None,               DKhat_BA_Dtheta],
                                                 [None,             DKhat_BA_Dtheta.T,  DKhat_AA_Dtheta]
                                         ]
                                       )

        DM_Dtheta=scipy.sparse.construct.bmat( [ [None,             None,               DMhat_mA_Dtheta],
                                                 [None,             None,               DMhat_BA_Dtheta],
                                                 [DMhat_mA_Dtheta.T,DMhat_BA_Dtheta.T,  DMhat_AA_Dtheta]
                                         ]
                                       )

        dKreduc.append(DK_Dtheta)
        dMreduc.append(DM_Dtheta)

    ##############################################################
    # FRF computation
    ##############################################################

    Flag_frf_analysis = 1
    frequencies = []
    frf = []
    frf_reduc = []
    frfgradient=list()
    for it in range(0,nbPara):
        frfgradient.append([])

    if (Flag_frf_analysis == 1):
        print("Proc. ", rank, " / time at the beginning of the FRF:", time.ctime())
        time0_frf=time.process_time()

        if rank == 0:
            print('nb of total dofs: ', len(SolvedDofF)+len(SolvedDofA))

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
            omega = 2*scipy.pi*freq

            print("Freq. step ",it,"/",itmax," proc number", rank, "frequency=", freq)

            tic = time.process_time()

            Freduc  = scipy.append(omega**2*Freduced_F,np.zeros( len(SolvedDofA) , dtype='c16' ))

            
            if rank>=0:
                 sol_reduc = mumps.spsolve(Kreduc-(omega**2)*Mreduc, Freduc,comm=None)#mycomm)
            
            ## Re-compute pressure with reduction strategy
            alpha_m = sol_reduc[list(range(nb_mode_F))].copy()
            P_B     = sol_reduc[list(range(nb_mode_F,nb_mode_F+len(SolvedDofB),1))].copy()
            P_A     = sol_reduc[list(range(nb_mode_F+len(SolvedDofB),nb_mode_F+len(SolvedDofB)+len(SolvedDofA),1))].copy()
            P_I     = eigen_vectors_I.dot(alpha_m)+Psi_IB.dot(P_B)+Psi_IA.dot(P_A)

            press_reduc      = np.zeros(fluid_ndof,dtype=complex)
            press_reduc[SolvedDofI] = P_I.copy()
            press_reduc[SolvedDofB] = P_B.copy()
            enrichment_reduc = np.zeros(fluid_ndof,dtype=complex)
            enrichment_reduc[SolvedDofA]= P_A.copy()
            CorrectedPressure=np.array(press_reduc.copy())
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA].T+np.array(enrichment_reduc[SolvedDofA]*np.sign(LevelSet[SolvedDofA]).T)

            ## compute and store FRF on the test volume
            ## frf_reduc.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(
            ##     fluid_elements5, fluid_nodes, press_reduc, enrichment_reduc, LevelSet, LevelSet*0-1.0))
            frf_reduc.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))

            ## correction of the pressure field with enrichment
            CorrectedPressure =press_reduc.copy() #np.zeros((fluid_ndof),dtype=complex) #press1.copy()
            CorrectedPressure[SolvedDofA] = press_reduc[SolvedDofA] + enrichment_reduc[SolvedDofA]*np.sign(LevelSet[SolvedDofA])
            
            if (flag_write_gmsh_results == 1) and (rank == 0):
                press_save.append(CorrectedPressure.copy())
                enrichment_save.append(enrichment_reduc.copy())
                uncorrectedpress_save.append(press_reduc.copy())
            
##            #####################
##            #####################
##            ######################
##            #####################
            Dpress_Dtheta = np.zeros([fluid_ndof,nbPara],dtype=complex)
            DCorrectedPressure_Dtheta=np.array(Dpress_Dtheta)
            Denrichment_Dtheta = np.zeros([fluid_ndof,nbPara],dtype=complex)
##            #####################
##            #####################
##            ## compute gradients
            for itP in range(0,nbPara):
                ## solve gradient problem
                tmp=-(dKreduc[itP]-(omega**2)*dMreduc[itP])*sol_reduc
                Dsol_Dtheta_RAW = mumps.spsolve(Kreduc-(omega**2)*Mreduc, tmp, comm=mycomm )
                Dalpha_m_Dtheta = Dsol_Dtheta_RAW[list(range(nb_mode_F))].copy()
                DP_B_Dtheta     = Dsol_Dtheta_RAW[list(range(nb_mode_F,nb_mode_F+len(SolvedDofB),1))].copy()
                DP_A_Dtheta     = Dsol_Dtheta_RAW[list(range(nb_mode_F+len(SolvedDofB),nb_mode_F+len(SolvedDofB)+len(SolvedDofA),1))].copy()
                DP_I_Dtheta     = eigen_vectors_I.dot(Dalpha_m_Dtheta)+Psi_IB.dot(DP_B_Dtheta)+Psi_IA.dot(DP_A_Dtheta)+DPsi_IA_Dtheta.dot(P_A)
                Dpress_Dtheta[SolvedDofI,itP] = DP_I_Dtheta
                Dpress_Dtheta[SolvedDofB,itP] = DP_B_Dtheta
                Denrichment_Dtheta[SolvedDofA,itP]= DP_A_Dtheta
                DCorrectedPressure_Dtheta[:,itP]=np.array(Dpress_Dtheta[:,itP].copy())
                DCorrectedPressure_Dtheta[SolvedDofA,itP]=DCorrectedPressure_Dtheta[SolvedDofA,itP].T+np.array(Denrichment_Dtheta[SolvedDofA,itP]*np.sign(LevelSet[SolvedDofA]).T)
##                #Dsol_Dtheta_RAW = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
##                # Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
##                #####################
##                #####################
##                ## gradient of the pressure field without enrichment
##                Dpress_Dtheta[SolvedDofF,itP] = Dsol_Dtheta_RAW[list(range(len(SolvedDofF)))].copy()
##                #####################
##                #####################
##                ## gradient of the enrichment field
##                Denrichment_Dtheta[SolvedDofA,itP]= Dsol_Dtheta_RAW[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))].copy()
##                #####################
##                #####################
##                #compute the corrected gradient pressure field (via enrichment)
##                DCorrectedPressure_Dtheta[:,itP]=np.array(Dpress_Dtheta[:,itP].copy())
##                DCorrectedPressure_Dtheta[SolvedDofA,itP]=DCorrectedPressure_Dtheta[SolvedDofA,itP].T+ \
##                    np.array(Denrichment_Dtheta[SolvedDofA,itP]*np.sign(LevelSet[SolvedDofA]).T)
                #####################
                #####################
                #store gradients
                frfgradient[itP].append(\
                    silex_lib_xfem_acou_tet4.computegradientcomplexquadratiquepressure(\
                        fluid_elements5,\
                            fluid_nodes,\
                                press_reduc+0j,\
                                    Dpress_Dtheta[:,itP]+0j,\
                                        LevelSet))                
##                #####################
##                #####################
                dpress_save[itP].append(DCorrectedPressure_Dtheta[:,itP].copy())
                denrichment_save[itP].append(Denrichment_Dtheta[:,itP].copy())
                duncorrectedpress_save[itP].append(Dpress_Dtheta[:,itP].copy())

        #print('nbPara',nbPara)
        #    for itP in range(0,nbPara):
        #        frfgradient[itP].append(1)# attention ici j'envoie frf dans frfgradient: pour tester
                
        frfsave=[frequencies,frf_reduc,frfgradient]
        if rank!=0:
            comm.send(frfsave, dest=0, tag=11)

        time1_frf=time.process_time()
        print("Proc. ", rank, " / time at the end of the FRF:", time.ctime())
        print("Time for FRF: ", time1_frf-time0_frf)
        print("Mean time for one freq. step: ", (time1_frf-time0_frf)/nbStep)
      
        # import pickle
        # data={'nodes':fluid_nodes,
        #       'elems':fluid_elements1,
        #       'press':press_save,
        #       'dpress':dpress_save}
        # output = open('out.db','wb')
        # p=pickle.Pickler(output)
        # p.dump(data)
        # output.close()
        

        if (flag_write_gmsh_results == 1) and (rank == 0):
            dataW=list()
            #prepare pressure field
            dataW.append([np.real(press_save),'nodal',1,'pressure (real)'])
            dataW.append([np.imag(press_save),'nodal',1,'pressure (imaginary)'])
            dataW.append([np.absolute(press_save),'nodal',1,'pressure (norm)'])
            #prepare gradient pressure field
            itG=0
            for itP in NamePara:
                dataW.append([np.real(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (real)'])
                dataW.append([np.imag(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (imaginary)'])
                dataW.append([np.absolute(dpress_save[itG]),'nodal',1,'pressure gradient '+itP+' (norm)'])
                itG=itG+1
            print("Write pressure field and gradients in msh file")
            silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',
                                        fluid_nodes, fluid_elements1, 4,dataW)
            print(">>> Done!!")

            #export results with discontinuities on .pos files
            varExport=np.vstack(uncorrectedpress_save).transpose()
            varExportC=np.vstack(press_save).transpose()
            varExportB=np.vstack(enrichment_save).transpose()
##            print("Write pressure field in pos file")
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.real(varExport),np.real(varExportB),'results/press_plus_real.pos','Pressure + Real')
##            os.system('cd results&&bzip2 *.pos')
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.imag(varExport),np.imag(varExportB),'results/press_plus_imag.pos','Pressure + Imag')
##            os.system('cd results&&bzip2 *.pos')
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.absolute(varExport),np.absolute(varExportB),'results/press_plus_abs.pos','Pressure + Abs')
##            os.system('cd results&&bzip2 *.pos')
##            # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.absolute(varExport)**2,np.absolute(varExportB)**2,'press_plus_square.pos')
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,np.real(varExport),-np.real(varExportB),'results/press_moins_real.pos','Pressure - Real')
##            os.system('cd results&&bzip2 *.pos')
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,np.imag(varExport),-np.imag(varExportB),'results/press_moins_imag.pos','Pressure - Imag')
##            os.system('cd results&&bzip2 *.pos')
##            silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,np.absolute(varExport),-np.absolute(varExportB),'results/press_moins_abs.pos','Pressure - Abs')
##            os.system('cd results&&bzip2 *.pos')


##            #export data
##            dataexport=list()
##            dataexport.append(fluid_nodes)
##            dataexport.append(fluid_elements1)
##            dataexport.append(LevelSet)
##            dataexport.append(varExport)
##            dataexport.append(varExportB)
##            f=open('debug_export.pck','wb')
##            pickle.dump(dataexport, f)
##            # print(Allfrfsave)
##            f.close()

            print(">>> Done!!")
##            #
##            f=open(results_file+'_results.frf','wb')
##            pickle.dump(frfsave, f)
##            # print(Allfrfsave)
##            f.close()

##            itG=0
##            for key,itP in enumerate(NamePara):
##                GvarExport=np.vstack(duncorrectedpress_save[key]).copy().transpose()
##                GvarExportB=np.vstack(denrichment_save[key]).copy().transpose()
##                print("Write gradient of pressure field in pos file (",itP,")")
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.real(GvarExport),np.real(GvarExportB),'results/Gpress_plus_'+itP+'_real.pos','Gpressure + '+itP+' Real')
##                os.system('cd results&&bzip2 *.pos')
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,np.imag(GvarExport),np.imag(GvarExportB),'results/Gpress_plus_'+itP+'_imag.pos','Gpressure + '+itP+' Imag')
##                os.system('cd results&&bzip2 *.pos')
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,np.real(GvarExport),-np.real(GvarExportB),'results/Gpress_moins_'+itP+'_real.pos','Gpressure - '+itP+' Real')
##                os.system('cd results&&bzip2 *.pos')
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,np.imag(GvarExport),-np.imag(GvarExportB),'results/Gpress_moins_'+itP+'_imag.pos','Gpressure - '+itP+' Imag')
##                os.system('cd results&&bzip2 *.pos')
##                #     
##                # gradPsquare=2*(np.real(GvarExport)*np.real(varExport)+np.imag(GvarExport)*np.imag(varExport))
##                # gradPsquareB=2*(np.real(GvarExportB)*np.real(varExportB)+np.imag(GvarExportB)*np.imag(varExportB))
##                # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,gradPsquare,gradPsquareB,'results/Gpress_plus_'+itP+'_dsquare.pos')                
##                # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,gradPsquare,-gradPsquareB,'results/Gpress_moins_'+itP+'_dsquare.pos')
##                #                          
##                gradCalc=(np.real(GvarExport)*np.real(varExport)+np.imag(GvarExport)*np.real(varExport))/np.absolute(varExport)
##                #deal with zeros values in enrichment field
##                gradCalcB=np.zeros([fluid_nnodes,nbStep])
##                gradCalcB[SolvedDofA,:]=(np.real(GvarExportB[SolvedDofA,:])*np.real(varExportB[SolvedDofA,:])+np.imag(GvarExportB[SolvedDofA,:])*np.imag(varExportB[SolvedDofA,:]))/np.absolute(varExportB[SolvedDofA,:])
##                #remove inf value
##                IX = np.absolute(varExport)==0.
##                gradCalc[IX]=1
##                gradCalcB[IX]=1
##                #
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,gradCalc,gradCalcB,'results/Gpress_plus_'+itP+'_dabsolute.pos','Gpressure + '+itP+' dAbs')       
##                os.system('cd results&&bzip2 *.pos')
##                silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,gradCalc,-gradCalcB,'results/Gpress_moins_'+itP+'_dabsolute.pos','Gpressure - '+itP+' dAbs')
##                os.system('cd results&&bzip2 *.pos')
##                print(">>> Done!!")
##                itG=itG+1

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
                    # print(data)
                else:
                    # print(i)
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

            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            # print(Allfrfsave)
            f.close()
##            #####################
##            #####################
##            #save on mat file
##            scipy.io.savemat(results_file+'_results.mat',mdict={'AllFRF': Allfrfsave})
##            scipy.io.savemat(results_file_ini+'results.mat',mdict={'AllFRF': Allfrfsave})
##            #####################
##            #####################
##            return Allfrfsave


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
    print ("Number of processors: ",nbProc)
    print ("Parameters: ",paraVal)
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
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,gradCompute,1)#,caseDefine)

#usage definition
def usage():
    dV=defaultV
    print("Usage: ",sys.argv[0],"-psFfhg [+arg]")
    print("\t -p : input parameters (default value ",dV.nbProc,")")
    print("\t -s : number of steps in the frequency range (default value ",dV.nbStep,")")
    print("\t -F : maximum frequency (default value ",dV.freqMax,")")
    print("\t -f : minimum frequency (default value ",dV.freqMin,")")
    print("\t -g : Components of grad (default value ",dV.gradCompute,")")

#default values
class defaultV:
    freqMin     = 10.0
    freqMax     = 150.0
    nbStep      = 149*2
    paraVal   = [1.,1.,0.5,0.8]
    gradCompute =  [0,1,2,3]
    nbProc=1
    #caseDef= 'thick_u'

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)

