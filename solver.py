###########################################################
# Declaration of useful solver tools for SILEXclass
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import logging
import time
import scipy
import scipy.sparse
import numpy as np
#
import utils
import structTools
#
from SILEX import silex_lib_xfem_acou_tet4

class solverTools:

    def solveLinear(self, A, B):
        """
        Method used to solve linear problem (choose automatically the available approach)
        """
        if self.mumpsOk:
            sol = self.classMumps.spsolve(A, B, comm=self.commMPI)
            # mumps is globaly defined
        else:
            sol = scipy.sparse.linalg.spsolve(A, B)
        return sol

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def buildLevelSet(self,typeLS=None,typeGEO=None,paraVal=None):
        """
        ##################################################################
        # function used to build the LevelSet (signed) and the gradient of the level-set
        ##################################################################
        """
        logging.info("Build Level-set: type %s"%(typeLS))
        tic = time.process_time()
        if typeLS is "FromMesh":
            # the level-set is built using the structure mesh
            if paraVal is not None:
                fileOrigGeo = self.getDatafile(typeFile='geo')
                fileDestMesh = self.getResultFile(detPara=True,addTxt=True,ext=None)
                structTools.buildStructMesh(fileOrigGeo,fileDestMesh,paraVal)
            # load mesh from file
            self.loadMesh(type='nodesStruct')
            self.loadMesh(type='elemsStruct')
            # build Level-Set from a structure mesh
            LevelSet,LevelSetU = silex_lib_xfem_acou_tet4.computelevelset(
                self.fluidNodes,
                self.structNodes,
                self.structElems)
            #compute gradient of Level-Set
            if paraData['gradCompute']:
                pass

        if typeLS is "manual":
            #the level-set is built using an explicit function 
            LSobj=structTools.LSmanual(typeName=typeGEO,nodes=self.fluidNodes,paraVal=paraVal)
            #export values
            self.LevelSet,self.LevelSetU=LSobj.exportLS()
            self.loadParaU(namePara = LSobj.exportParaName())
            #compute gradient of Level-Set
            if self.paraData['gradCompute']:
                self.LevelSetGradient,nameGrad= LSobj.exportLSgrad(self.getNameGrad())
                self.loadParaU(nameParaGrad = nameGrad)
        #
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    
    def buildEnrichedPart(self):
        """
        ##################################################################
        # Find the enriched parts using the Level-set
        ##################################################################
        """
        logging.info("Find enriched elements and nodes using LS")
        tic = time.process_time()
        #
        self.EnrichedElems, self.EnrichedNbElems = silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(self.fluidElems, self.LevelSet)
        self.EnrichedElems = self.EnrichedElems[list(range(self.EnrichedNbElems))]
        #
        # self.LSEnrichedElems=self.LSEnrichedElems#[self.LSEnrichedElems-1]
        self.EnrichedNodes = np.unique(self.fluidElems[self.EnrichedElems])
        #
        # tmp = []
        # for i in LSEnrichednodes:
        #     for j in range(4):
        #         tmpp = scipy.where(fluid_elements1[:, j] == i)[0]
        #         for k in range(len(tmpp)):
        #             tmp.append(tmpp[k])
        # tmp.append(scipy.where(fluid_elements1[:,1]==i))
        # tmp.append(scipy.where(fluid_elements1[:,2]==i))
        # tmp.append(scipy.where(fluid_elements1[:,3]==i))
        #
        # tmp = scipy.unique(scipy.array(tmp))
        # tmp1,elttest0,tmp2=scipy.intersect1d(fluid_elements1[:,0],LSEnrichednodes,return_indices=True)
        # silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements_test0',fluid_nodes,fluid_elements1[tmp],4)
        #[75804, 97252, 97253,34973, 93135, 93137, 93248,83787, 93136,93525]
        # EnrichedElements0, NbEnrichedElements = silex_lib_xfem_acou_tet4.getsurfenrichedelements(
        #     struc_nodes, struc_elements, fluid_nodes, fluid_elements1[tmp])
        # EnrichedElements0 = scipy.unique(
        #     EnrichedElements0[list(range(NbEnrichedElements))])
        # EnrichedElements0 = EnrichedElements0-1
        # EnrichedElements = tmp[EnrichedElements0]
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def buildSecondMember(self):
        """
        ##############################################################
        # Build second member
        ##############################################################
        """    
        logging.info("Build second member")
        tic = time.process_time()
        #fluid vector second member
        self.UF = np.zeros(2*self.fluidNbDofs, dtype=float)
        # apply bc in displacement using data
        for bc in self.caseProp['bcdisp']:
            #get number of dofs and values to apply
            numDofs,valbc = self.getFormatBC(bc)
            self.UF[numDofs]=valbc
        #
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

        if self.paraData['gradCompute']:
            logging.info("Build gradient of second member")
            tic = time.process_time()
            for it in self.paraData['nameGrad']:
                self.dUF.append(None) 
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def buildFluidOperators(self):
        """
        ##############################################################
        # Compute Standard Fluid Matrices
        ##############################################################
        """
        logging.info("Build fluid operators")
        tic = time.process_time()   
        # build matrices using vector description
        IIf, JJf, Vffk, Vffm = silex_lib_xfem_acou_tet4.globalacousticmatrices(
            self.fluidElems, 
            self.fluidNodes, 
            self.mechaProp['celerity'], 
            self.mechaProp['rho'])

        self.KFF = scipy.sparse.csc_matrix(
            (Vffk, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        self.MFF = scipy.sparse.csc_matrix(
            (Vffm, (IIf, JJf)), 
            shape=(self.fluidNbDofs, self.fluidNbDofs))

        self.SolvedDofF = list(range(self.fluidNbDofs))
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
        #
        if self.paraData['gradCompute']:
            logging.info("Build fluid gradient operators")
            tic = time.process_time()
            for it in self.paraData['nameGrad']:
                self.dKFF.append(None) 
                self.dMFF.append(None) 
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    
    def buildEnrichedOperators(self):
        """
        ##################################################################
        # Compute Heaviside enrichment
        ##################################################################
        """
        logging.info("Build enriched operators")
        tic = time.process_time()  
        # build matrices using vector description
        IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(
            self.fluidElems, 
            self.fluidNodes, 
            self.LevelSet, 
            self.mechaProp['celerity'], 
            self.mechaProp['rho'])

        self.KAA = scipy.sparse.csc_matrix(
            (Vaak, (IIaa, JJaa)), 
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        self.MAA = scipy.sparse.csc_matrix(
            (Vaam, (IIaa, JJaa)), 
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        self.KAF = scipy.sparse.csc_matrix(
            (Vafk, (IIaf, JJaf)), 
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        self.MAF = scipy.sparse.csc_matrix(
            (Vafm, (IIaf, JJaf)), 
            shape=(self.fluidNbDofs, self.fluidNbDofs))

        self.SolvedDofA = self.EnrichedNodes-1
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

        if self.paraData['gradCompute']:
            logging.info("Build gradient of enriched operators")
            tic = time.process_time()
            for it in range(0,self.getNbGrad()):
                #compute gradients matrices for each selected parameter
                dKAAx,dKFAx,dMAAx,dMFAx=self.buildGoperators(self,it)
                #store
                self.dKAA.append(dKAAx) 
                self.dKFA.append(dKAAx) 
                self.dMAA.append(dKAAx) 
                self.dMFA.append(dKAAx) 
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def buildOperators(self):
        """
        #################################################################
        # Construct the whole system
        #################################################################
        """
        logging.info("Assembly of operators of the whole problem")
        tic = time.process_time()  
        #
        fd=self.mechaProp['fluid_damping']
        #
        self.K = scipy.sparse.construct.bmat([
            [fd*self.KFF[self.SolvedDofF, :][:, self.SolvedDofF], fd*self.KAF[self.SolvedDofF, :][:, self.SolvedDofA]],
            [fd*self.KAF[self.SolvedDofA, :][:, self.SolvedDofF], fd*self.KAA[self.SolvedDofA, :][:, self.SolvedDofA]]])
        #
        self.M = scipy.sparse.construct.bmat([
            [self.MFF[self.SolvedDofF, :][:, self.SolvedDofF], self.MAF[self.SolvedDofF, :][:, self.SolvedDofA]],
            [self.MAF[self.SolvedDofA, :][:, self.SolvedDofF], self.MAA[self.SolvedDofA, :][:, self.SolvedDofA]]])
        #
        #list of all dofs
        self.SolvedDof = np.hstack([self.SolvedDofF, self.SolvedDofA+self.fluidNbDofs])
        #
        logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
        #
        if self.paraData['gradCompute']:
            logging.info("Build gradient of the assembled operators")
            tic = time.process_time()            
            for it in range(0,self.getNbGrad()):
                # build full stiffness and mass gradient matrices
                dK.append(scipy.sparse.construct.bmat([
                    [fd*self.dKFF[it][self.SolvedDofF, :][:, self.SolvedDofF], fd*self.dKFA[it][self.SolvedDofF, :][:, self.SolvedDofA]],
                    [fd*self.dKFA[it][self.SolvedDofA, :][:, self.SolvedDofF], fd*self.dKAA[it][self.SolvedDofA, :][:, self.SolvedDofA]]]))
                dM.append(scipy.sparse.construct.bmat([
                    [self.dMFF[it][self.SolvedDofF, :][:, self.SolvedDofF], self.dMFA[it][self.SolvedDofF, :][:, self.SolvedDofA]],
                    [self.dMFA[it][self.SolvedDofA, :][:, self.SolvedDofF], self.dMAA[it][self.SolvedDofA, :][:, self.SolvedDofA]]]))
            #
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def buildGOperators(self,itG):
        """
        ##################################################################
        # Build gradients operators one by one
        ##################################################################
        """
        logging.info(' Build gradient matrices for parameter '+self.paraData['NameGrad'][itG])
        #compute gradients matrices and indices
        # print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
        IIf, JJf, Vfak_gradient, Vfam_gradient =\
            silex_lib_xfem_acou_tet4.globalacousticgradientmatrices(self.fluidNodes,
                                                                    self.fluidElems,
                                                                    self.LevelSet,
                                                                    self.paraData['celerity'],
                                                                    self.paraData['rho'],
                                                                    self.LevelSetGradient[itG])
        #build sparse matrices
        dKAA_dtheta=None
        dMAA_dtheta=None
        dKFA_dtheta = scipy.sparse.csc_matrix(
            (Vfak_gradient, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        dMFA_dtheta = scipy.sparse.csc_matrix(
            (Vfam_gradient, (IIf, JJf)),
             hape=(self.fluidNbDofs, self.fluidNbDofs))
        #
        return dKAA_dtheta,dKFA_dtheta,dMAA_dtheta,dMFA_dtheta

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def preProcessMaster(self):
        """
        ##################################################################
        # method to initialize data for master node (declare ie rank=0) and 
        # operators that will remain constant for all parameters
        ##################################################################
        """
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##    Start SILEX vibro-acoustics computation   ##")
        if self.paraData['gradCompute']:
            logging.info("##          (with gradients computations)        ##")            
        logging.info("##################################################")
        logging.info("##################################################")
        
        # check & create folders/files 
        self.createDatabase()
        # load MPI info
        self.loadMPI()
        # load MUMPS
        self.loadMUMPS()
        # load fluid mesh
        self.loadMesh(typeData='nodesFluid',dispFlag=True,filename=self.getDatafile('fluidmesh'))
        self.loadMesh(typeData='elemsFluid',dispFlag=True,filename=self.getDatafile('fluidmesh'))
        # load control volume
        self.loadMesh(typeData='elemsControlFluid',dispFlag=True,filename=self.getDatafile('fluidmesh'))
        #build fluid operators
        self.buildFluidOperators() 
        #build second member
        self.buildSecondMember()
        #generate the list of frequencies
        self.Frequencies=utils.computeFreqPerProc(\
            self.caseProp['nbSteps'],
            self.nbProcMPI,
            self.caseProp['freqMin'],
            self.caseProp['freqMax'])

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    
    def initRun(self,paraVal=None):
        """
        ##################################################################
        # method used to intialize the data before a run associated to a set of parameters
        ##################################################################
        """
        # initialize data
        self.initDataSolve()
        #
        self.showDataParaVal()
        #build levelset
        self.buildLevelSet(typeLS=self.caseProp['typeLS'],typeGEO=self.caseProp['typeGEOstruct'],paraVal=paraVal)
        #build enriched part
        self.buildEnrichedPart()
        #build operators
        self.buildEnrichedOperators()
        self.buildOperators()

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def solvePbOneStep(self,itF,itMax,freq):
        """
        Method used to solve the problem for one frequency step
        """
        tic = time.process_time()
        #
        logging.info("Rank: %i - Solve whole problem for frequency %g (step %i/%i)"%(self.rankMPI,freq,itF,itMax))
        # compute natural frequency
        omega=2*np.pi*freq
        #Build full second member
        F=(omega**2)*np.array(self.UF[self.SolvedDof],dtype=self.loadType())
        #solve the whole problem on pressure
        ticB = time.process_time()
        sol = self.solveLinear(self.K-(omega**2)*self.M,F)
        logging.info("Rank: %i - Solve linear problem - Done - %g s"%(self.rankMPI,time.process_time()-ticB))
        ##### build the fields
        ticB = time.process_time()
        # uncorrected pressure field        
        self.pressureUncorrect[itF][self.SolvedDofF] = sol[list(range(len(self.SolvedDofF)))].copy()
        # enrichment field
        self.pressureEnrichment[itF][self.SolvedDofA] = sol[list(range(len(self.SolvedDofF),len(self.SolvedDofF)+len(self.SolvedDofA)))].copy()
        # corrected pressure field
        self.pressure[itF] = self.pressureUncorrect[itF].copy()
        self.pressure[itF][self.SolvedDofA] += self.pressureEnrichment[itF][self.SolvedDofA]*np.sign(self.LevelSet[self.SolvedDofA])
        logging.info("Rank: %i - Fields computation - Done - %g s"%(self.rankMPI,time.process_time()-ticB))
        #compute and store FRF on the control volume
        if self.caseProp['computeFRF']:
            self.FRF[itF] = silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(
                self.fluidElemsControl,
                self.fluidNodes,
                self.pressureUncorrect[itF],
                self.pressureEnrichment[itF],
                self.LevelSet,
                self.LevelSet*0.-1.0)
        #Compute gradients of fields
        if self.paraData['gradCompute']:
            #initialize variables
            
            # along the parameters
            logging.info("Rank: %i - Start gradients computation"%(self.rankMPI))
            ticC = time.process_time()
            for itG in range(0,self.paraData['nbGrad']):
                ticB = time.process_time()
                #prepare data
                tmp = -(self.dK[itG]-(omega**2)*self.dM[itG])*sol
                Dsol = self.solveLinear(self.K-(omega**2)*self.M,tmp)
                logging.info("Rank: %i - Prepare and solve linear problem for grad (var: %s/num: %i) - Done - %g s"%(self.rankMPI,self.paraData['nameGrad'],itG,time.process_time()-ticB))
                # gradient of uncorrected pressure field
                self.pressureUncorrectGrad[itG][self.SolvedDofF,itF] = Dsol[list(range(len(self.SolvedDofF)))].copy()
                # gradient of enrichment field
                self.pressureUncorrectGrad[itG][self.SolvedDofA,itF] = Dsol[list(range(len(self.SolvedDofF),len(self.SolvedDofF)+len(self.SolvedDofA)))].copy()
                # gradient of corrected pressure field
                self.pressureGrad[itG][:,itF] = self.pressureUncorrectGrad.copy()
                self.pressureGrad[itG][self.SolvedDofA,itG] += self.pressureEnrichmentGrad[itF][self.SolvedDofA,itG]*np.sign(self.LevelSet[self.SolvedDofA].T)
                #compute and store gradient of FRF on the control volume
                if self.caseProp['computeFRF']:
                    self.FRFgrad[itG][itF] = silex_lib_xfem_acou_tet4.computegradientcomplexquadratiquepressure(
                        self.fluidElemsControl,
                        self.fluidNodes,
                        self.pressureUncorrect[itF],
                        self.pressureUncorrectGrad[itG][:,itF],
                        self.LevelSet)
            logging.info("Rank: %i - Gradients computation - Done - %g s"%(MPI,time.process_time()-ticC))

        logging.info("Rank: %i - Done - Solve whole problem for frequency %g (step %i/%i) in %g s"%(self.rankMPI,freq,itF,itMax,time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def initDataSolve(self):
        """
        ##################################################################
        # Initialize storage for computation
        ##################################################################
        """
        self.pressureUncorrect = [np.zeros((self.fluidNbDofs),dtype=self.loadType()) for _ in range(self.caseProp['nbSteps'])]
        self.pressureEnrichment = [np.zeros((self.fluidNbDofs),dtype=self.loadType()) for _ in range(self.caseProp['nbSteps'])]
        self.pressure = [None for _ in range(self.caseProp['nbSteps'])]
        #
        if self.caseProp['computeFRF']:
            self.FRF = np.zeros(self.caseProp['nbSteps'])
        #
        if self.paraData['gradCompute']:
            self.pressureUncorrectGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureEnrichmentGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            #
            self.FRFgrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
        # initialize class for exporting meshes
        self.classSave = None

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def solvePb(self,paraVal=None):
        """
        ##################################################################
        # Method used to solve the whole problem
        ##################################################################
        """
        ticS = time.process_time()
        if self.nbRuns == 0:
            self.preProcessMaster()
        if paraVal is not None:
            paraVal=np.array(paraVal)
        else:
            logging.error('>> Unable to start computation due to no given parameters values')
            raise
        #prepare parameters values to run
        paraValOk=self.prepPara(paraVal)
        #
        # along the parameters
        for valU in paraValOk:
            self.nbRuns += 1
            ticV = time.process_time()
            logging.info("##################################################")
            self.paraData['val']=valU            
            txtPara=self.formatPara(valIn=valU)
            logging.info('Start compute for parameters (nb %i): %s'%(self.nbRuns,txtPara))
            #initialization for run
            self.initRun(paraVal=valU)

            # along the frequencies
            for (itF,Freq) in enumerate(self.getFrequencies()):
                # solve the problem
                self.solvePbOneStep(itF,len(self.getFrequencies(total=True)),Freq)
            #export results (FRF, pressure fields and gradients)
            if self.flags['saveResults']:
                self.exportFieldsOnePara(paraName=True)
                self.saveFRF(paraName=True)
            # append dat to history
            self.paraData['oldval'].append(valU)
            self.allFRF.append(self.FRF)
            self.allFRFgrad.append(self.FRFgrad)
            #
            logging.info("Time to solve the whole problem for set of parameters nb %i - %g s"%(self.nbRuns,time.process_time()-ticV))
        #plot and export results (FRF, pressure fields and gradients)
        if self.flags['saveResults']:
            self.saveFRF(paraName=False,allData=True)
            self.plotFRF(allData=True,fileOut=self.getResultFile(detPara=False,addTxt='allData',ext='csv'))
            
        #
        logging.info("Time to solve the whole problem along sets of parameters - %g s"%(time.process_time()-ticS))
