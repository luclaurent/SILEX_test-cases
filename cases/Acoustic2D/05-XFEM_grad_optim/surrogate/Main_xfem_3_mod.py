import string
import time
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

def RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,positionStruct):
    # parallepipedic cavity with plane structure
    mesh_file='geom/xfem3_noarea_fin'
    results_file_ini='results/xfem_3_'

    listFreqPerProc = computeFreqPerProc(nbStep,nbProc,freqMin,freqMax)

    ##############################################################
    # Material, Boundary conditions
    ##############################################################

    # air
    celerity=340.0
    rho=1.2

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

    h_struc=1.0

    ValpStruct=positionStruct*10000.

    file_extension=str(ValpStruct)[0:5]
    results_file=results_file_ini+file_extension
    print(results_file)

    x_pos_struc=positionStruct
    print(x_pos_struc)

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

    fluid_elements_boun,IdnodeS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',1,2)

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

    LevelSet=fluid_nodes[:,0]-x_pos_struc
    LevelSetTangent=fluid_nodes[:,1]-h_struc

    # level set gradient with respect to parameters
    LevelSet_gradient=-scipy.ones(fluid_nnodes)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,2,[[LevelSetTangent,'nodal',1,'Tangent level set']])
        silex_lib_gmsh.WriteResults(results_file+'_level_set',fluid_nodes,fluid_elements,2,[[LevelSet,'nodal',1,'Level set']])

    toc = time.clock()
    print("time to compute level set:",toc-tic)


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.clock()

    struc_nodes=scipy.array([[x_pos_struc,0.0],[x_pos_struc,h_struc/2.0],[x_pos_struc,h_struc]])
    struc_elements=scipy.array([[1,2],[2,3]])
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
    #SolvedDofF=range(fluid_ndof)

    toc = time.clock()
    print("time to compute fluid matrices:",toc-tic)

    ##################################################################
    # Compute enrichment: Heaviside + Edge
    ##################################################################
    tic = time.clock()

    #HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

    #Enrichednodes = scipy.unique(fluid_elements[HeavisideEnrichedElements])
    #Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

    #print xvibacoufo.getpositivenegativeelts.__doc__

    NegativeLSelements,PositiveLSelements,NegativeLStgtElements,PositiveLStgtElements,nbNegLS,nbPosLS,nbNegLSt,nbPosLSt=silex_lib_tri3_acou.getpositivenegativeelts(fluid_elements,LevelSet,LevelSetTangent)

    NegativeLSelements=NegativeLSelements[list(range(nbNegLS))]
    PositiveLSelements=PositiveLSelements[list(range(nbPosLS))]
    NegativeLStgtElements=NegativeLStgtElements[list(range(nbNegLSt))]
    PositiveLStgtElements=PositiveLStgtElements[list(range(nbPosLSt))]

    EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_tri3_acou.getenrichedelementsfromlevelset(fluid_elements,LevelSetTangent)
    EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])

    IdElementTip=silex_lib_tri3_acou.getelementcontainingpoint(fluid_elements,fluid_nodes,[0.6,0.65])

    #if (flag_write_gmsh_results==1) and (rank==0):
    #    silex_lib_gmsh.WriteResults(results_file+'_NegativeLSelements',fluid_nodes,fluid_elements[NegativeLSelements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_PositiveLSelements',fluid_nodes,fluid_elements[PositiveLSelements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_NegativeLStgtElements',fluid_nodes,fluid_elements[NegativeLStgtElements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_PositiveLStgtElements',fluid_nodes,fluid_elements[PositiveLStgtElements],2)




    IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_tri3_acou.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,LevelSetTangent,celerity,rho,flag_edge_enrichment)

    KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
    MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

    toc = time.clock()
    print("time to compute Heaviside enrichment:",toc-tic)

    #Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([HeavisideEnrichedElements,EdgeEnrichedElements]))])
    #Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements,EdgeEnrichedElementsInAllMesh]))])
    #Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements]))])
    #Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([NegativeLStgtElements]))])
    Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])
    #Enrichednodes = scipy.unique(fluid_elements)
    SolvedDofA=Enrichednodes-1

    silex_lib_gmsh.WriteResults(results_file+'_EnrichedElements',fluid_nodes,fluid_elements[scipy.hstack(([EnrichedElements]))],2)

    #################################################################
    # Construct the whole system
    ##################################################################

    K=scipy.sparse.construct.bmat( [[KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
                                    [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )



    M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                    [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )

    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################

    print(silex_lib_tri3_acou.globalacousticgradientmatrices.__doc__)

    IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient,LevelSet)
    #IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements,fluid_nodes,celerity,rho,LevelSet_gradient,LevelSet)
    dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    f=open(results_file+'_KAF_MAF.pck','wb')
    pickle.dump([dKFA_dtheta[SolvedDofF,:][:,SolvedDofA],dMFA_dtheta[SolvedDofF,:][:,SolvedDofA],KAF[SolvedDofF,:][:,SolvedDofA],MAF[SolvedDofF,:][:,SolvedDofA]], f)
    f.close()

    dK=scipy.sparse.construct.bmat( [[None,dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM=scipy.sparse.construct.bmat( [[None,dMFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )


    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis=1
    FF = scipy.zeros(fluid_ndof)
    frequencies=[]
    frf=[]
    frfgradient=[]

    if (Flag_frf_analysis==1):
        print("time at the beginning of the FRF:",time.ctime())

        press_save=[]
        dpress_save=[]
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

            #sol = mumps.spsolve(scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float'), F, comm=mycomm )
            sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float'), F)
            
            #eigen_values_small,eigen_vectors= scipy.sparse.linalg.eigsh(K-(omega**2)*M,1,sigma=0,which='SM')
            #eigen_values_large,eigen_vectors= scipy.sparse.linalg.eigsh(K-(omega**2)*M,1,sigma=0,which='LM')


            #eigen_values_small,eigen_vectors= scipy.linalg.eig((K-(omega**2)*M).todense(),1,sigma=0,which='SM')
            #eigen_values_large,eigen_vectors= scipy.linalg.eig(K-(omega**2)*M,1,sigma=0,which='LM')

            
            #stop
            tmp=-(dK-(omega**2)*dM)*sol

            #Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
            Dsol_Dtheta = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='float')  , tmp )

            press = scipy.zeros(fluid_ndof)
            press[IdnodeS2-1] = scipy.ones(len(IdnodeS2))
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            enrichment=scipy.zeros(fluid_nnodes)
            enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            CorrectedPressure=press
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
            #frf.append(silex_acou_lib_tri3.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
            frf.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,press+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            #frf[i]=xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent)
            #press_save.append(CorrectedPressure.copy())
            press_save.append(press.copy())

            Dpress_Dtheta = scipy.zeros(fluid_ndof,dtype=float)
            Dpress_Dtheta[SolvedDofF] = Dsol_Dtheta[list(range(len(SolvedDofF)))]
            Denrichment_Dtheta = scipy.zeros(fluid_ndof,dtype=float)
            Denrichment_Dtheta[SolvedDofA]= Dsol_Dtheta[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            #DCorrectedPressure_Dtheta=scipy.array(Dpress_Dtheta)
            #DCorrectedPressure_Dtheta[SolvedDofA]=DCorrectedPressure_Dtheta[SolvedDofA].T+scipy.array(Denrichment_Dtheta[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)
            frfgradient.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,press+0j,Dpress_Dtheta+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            dpress_save.append(Dpress_Dtheta.copy())
        

        print("time at the end of the FRF:",time.ctime())
        frfsave=[frequencies,frf,frfgradient]
        if rank!=0 :
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results==1) and (rank==0):
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,2,[[press_save,'nodal',1,'pressure'],[dpress_save,'nodal',1,'pressure gradient']])

        # save the FRF problem
        Allfrequencies=scipy.zeros(nbStep)
        Allfrf=scipy.zeros(nbStep)
        Allfrfgradient=scipy.zeros(nbStep)
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
                    Allfrfgradient[k]=data[2][j]
                    k=k+1

            Allfrequencies, Allfrf,Allfrfgradient = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf)),scipy.array(list(Allfrfgradient))]#,press_save,dpress_save]
            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            print(Allfrfsave)
            f.close()
            #save on mat file
            scipy.io.savemat(results_file+'_results.mat',mdict={'AllFRF': Allfrfsave})
            scipy.io.savemat(results_file_ini+'results.mat',mdict={'AllFRF': Allfrfsave})


#function for dealing with options
def manageOpt(argv,dV):
    #load default values
    freqMin     = dV.freqMin
    freqMax     = dV.freqMax
    nbStep      = dV.nbStep
    posStruct   = dV.posStruct
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
            posStruct = float(arg)
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    print ("Number of processors: ",nbProc)
    print ("Number of frequency steps: ",nbStep)
    print ("Maximum frequency: ",freqMax)
    print ("Minimum frequency: ",freqMin)
    print ("Strucure position: ",posStruct)
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,posStruct)

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
    posStruct   = 0.0001

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)
