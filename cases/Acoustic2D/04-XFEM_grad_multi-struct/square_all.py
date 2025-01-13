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

import buildStruct as lvlB

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

def RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,caseDefine):
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


    # number of parameters
    nbPara=len(paraVal);
    #prepare save file
    file_extension="{:.{}E}".format(paraVal[0],2)
    if nbPara>1:
        for i in range(1,nbPara):
            file_extension=file_extension+'_'+"{:.{}E}".format(paraVal[i],2)

    results_file=results_file_ini+file_extension
    print(results_file)

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

    #build level set
    struc_nodes,struc_elements,LevelSet,LevelSetGradient,NamePara,LevelSetTangent,IndicZone=lvlB.buildStruct(paraVal,fluid_nodes,caseDefine)

    #number of parameters
    nbPara=len(NamePara)

    #export levelset data in a gmsh file
    if (flag_write_gmsh_results==1) and (rank==0):
        #export levelset tangent
        silex_lib_gmsh.WriteResults(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,2,[[LevelSetTangent,'nodal',1,'Tangent level set']])
        #export levelset and levelset gradient
        dataW=list()
        dataW.append([LevelSet,'nodal',1,'Level set'])
        itP=0
        for iN in NamePara:
            dataW.append([LevelSetGradient[itP],'nodal',1,'Level Set Grad '+iN])
            itP=itP+1

        dataW.append([IndicZone,'nodal',1,'Indic Zones'])
        silex_lib_gmsh.WriteResults(results_file+'_level_set',fluid_nodes,fluid_elements,2,dataW)

    toc = time.clock()
    print("time to compute level set:",toc-tic)


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.clock()

    
    struc_boun=scipy.array([3])

    silex_lib_gmsh.WriteResults(results_file+'_struc_mesh',struc_nodes,struc_elements,1)

    #EnrichedElements,NbEnrichedElements=silex_lib_tri3_acou.getenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements)
    EnrichedElements,NbEnrichedElements=silex_lib_tri3_acou.getenrichedelementsfromlevelset(fluid_elements,LevelSet)
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

    #build gradients matrices
    dK=list()
    dM=list()
    for itP in range(0,nbPara):
        print(' Build Matrix for parameter '+NamePara[itP])
        IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSetGradient[itP],LevelSet)
        dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
        #export data
        f=open(results_file+'_KAF_MAF.pck','wb')
        pickle.dump([dKFA_dtheta[SolvedDofF,:][:,SolvedDofA],dMFA_dtheta[SolvedDofF,:][:,SolvedDofA],KAF[SolvedDofF,:][:,SolvedDofA],MAF[SolvedDofF,:][:,SolvedDofA]], f)
        f.close()
        #build full stiffness and mass gradient matrices
        dK.append(scipy.sparse.construct.bmat( [[None,fluid_damping*dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [fluid_damping*dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )
                                    )

        dM.append(scipy.sparse.construct.bmat( [[None,dMFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )
                                    )


    

    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis=1
    FF = scipy.zeros(fluid_ndof)
    frequencies=[]
    frf=[]
    frfgradient=list()
    for it in range(0,nbPara):
        frfgradient.append([])

    if (Flag_frf_analysis==1):
        print("time at the beginning of the FRF:",time.ctime())

        press_save=[]
        dpress_save=list()
        for it in range(0,nbPara):
            dpress_save.append([])

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
            ## solve direct problem
            sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(pbFreq,dtype=complex), F)
            #store the pressure field
            press = scipy.zeros(fluid_ndof,dtype=complex)
            press[IdnodeS2-1] = scipy.ones(len(IdnodeS2))
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            enrichment=scipy.zeros(fluid_nnodes,dtype=complex)
            enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            #compute the corrected pressure field (via enrichment)
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
            ######################
            #####################
            Dpress_Dtheta = scipy.zeros([fluid_ndof,nbPara],dtype=complex)
            DCorrectedPressure_Dtheta=scipy.array(Dpress_Dtheta)
            Denrichment_Dtheta = scipy.zeros([fluid_ndof,nbPara],dtype=complex)
            ## compute gradients
            for itP in range(0,nbPara):
                tmpG=-(dK[itP]-(omega**2)*dM[itP])*sol
                #Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )
                Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve(  scipy.sparse.csc_matrix(pbFreq,dtype=complex)  , tmpG )
                #####################
                #####################
                Dpress_Dtheta[SolvedDofF,itP] = Dsol_Dtheta_RAW[list(range(len(SolvedDofF)))]
                #####################
                #####################
                Denrichment_Dtheta[SolvedDofA,itP]= Dsol_Dtheta_RAW[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
                #####################
                #####################
                #compute the corrected gradient pressure field (via enrichment)
                DCorrectedPressure_Dtheta[:,itP]=scipy.array(Dpress_Dtheta[:,itP])
                DCorrectedPressure_Dtheta[SolvedDofA,itP]=DCorrectedPressure_Dtheta[SolvedDofA,itP].T+scipy.array(Denrichment_Dtheta[SolvedDofA,itP]*scipy.sign(LevelSet[SolvedDofA]).T)
                #####################
                #####################
                #store gradients
                frfgradient[itP].append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,CorrectedPressure+0j,DCorrectedPressure_Dtheta[:,itP]+0j,0.0*enrichment+0j,0.0*Denrichment_Dtheta[:,itP]+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
                #####################
                #####################
                dpress_save[itP].append(DCorrectedPressure_Dtheta[:,itP].copy())
        
        #####################
        #####################
        print("time at the end of the FRF:",time.ctime())
        frfsave=[frequencies,frf,frfgradient]
        if rank!=0 :
            comm.send(frfsave, dest=0, tag=11)

        if (flag_write_gmsh_results==1) and (rank==0):
            dataW=list()
            #prepare pressure field
            pR=numpy.real(press_save)
            pI=numpy.imag(press_save)
            pA=numpy.absolute(press_save)
            dataW.append([pR,'nodal',1,'pressure (real)'])
            dataW.append([pI,'nodal',1,'pressure (imaginary)'])
            dataW.append([pA,'nodal',1,'pressure (norm)'])
            #prepare gradient pressure field
            itG=0
            for itP in NamePara:
                dpR=numpy.real(dpress_save[itG])
                dpI=numpy.imag(dpress_save[itG])
                dpA=numpy.absolute(dpress_save[itG])
                dpC=(dpI*pI+dpR*pR)/pA
                dataW.append([dpR,'nodal',1,'pressure gradient '+itP+' (real)'])
                dataW.append([dpI,'nodal',1,'pressure gradient '+itP+' (imaginary)'])
                dataW.append([dpA,'nodal',1,'pressure gradient '+itP+' (norm)'])
                dataW.append([dpC,'nodal',1,'pressure gradient '+itP+' (deriv. absolute)'])
                itG=itG+1
            #write fields
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,2,dataW)

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
    caseDefine = dV.caseDef
    
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
        elif opt == "-h":
            usage()
            sys.exit()
    #print chosen parameters
    print ("Number of processors: ",nbProc)
    print ("Number of frequency steps: ",nbStep)
    print ("Maximum frequency: ",freqMax)
    print ("Minimum frequency: ",freqMin)
    print ("Case: ",caseDefine)
    it=0
    for itP in paraVal:
        print ('Parameter num '+str(it)+': '+str(itP))
        it=it+1
    print ("\n\n")

    #run computation
    RunPb(freqMin,freqMax,nbStep,nbProc,rank,comm,paraVal,caseDefine)

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
    nbStep      = 11
    paraVal   = [1.5,1,0.75,0]
    caseDef= 'thick_u'

### Run autonomous
if __name__ == '__main__':
    #run with options
    dV=defaultV
    manageOpt(sys.argv[1:],dV)
