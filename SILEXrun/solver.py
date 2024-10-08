###########################################################
# Declaration of useful solver tools for SILEXclass
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import time
import scipy
import scipy.sparse.linalg
import numpy as np
import logging
#
from .misc import utils
from . import structTools
from SILEXlib import silex_lib_xfem as libXFEM
#
# activate logger
Logger = logging.getLogger(__name__)


class solverTools(object):

    def solveLinear(self, A, B):
        """
        Method used to solve linear problem
        (choose automatically the available approach)
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

    def buildLevelSet(self, typeLS=None, typeGEO=None, paraVal=None):
        """
        ##################################################################
        # function used to build the LevelSet (signed)
        # and the gradient of the level-set
        ##################################################################
        """
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Build Level-set: type {}'.format(typeLS))
        tic = time.process_time()
        if typeLS == 'FromMesh':
            # the level-set is built using the structure mesh
            if paraVal is not None:
                fileOrigGeo = self.getDatafile(typeFile='geo')
                fileDestMesh = self.getResultFile(detPara=True,
                                                  addTxt=True,
                                                  ext=None)
                structTools.buildStructMesh(fileOrigGeo,
                                            fileDestMesh,
                                            paraVal)
            
            # load mesh from file
            self.loadMesh(type='nodesStruct')
            self.loadMesh(type='elemsStruct')
            # build Level-Set from a structure mesh
            self.LevelSet, self.LevelSetU = self.libXFEM.computeLevelSet(
                self.fluidNodes,
                self.structNodes,
                self.structElems)
            # compute gradient of Level-Set
            if self.paraData['gradCompute']:
                pass

        if typeLS == 'manual':
            # the level-set is built using an explicit function
            LSobj = structTools.LSmanual(typeName=typeGEO,
                                         nodes=self.fluidNodes,
                                         paraVal=paraVal)
            # export values
            self.LevelSet, self.LevelSetU = LSobj.exportLS()
            self.LevelSetTangent = LSobj.exportLST()
            self.loadParaU(namePara=LSobj.exportParaName())
            # compute gradient of Level-Set
            if self.paraData['gradCompute']:
                self.LevelSetGradient, nameGrad = LSobj.exportLSgrad(
                    self.getNameGrad())
                self.loadParaU(nameParaGrad=nameGrad)
        #
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))

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
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Find enriched elements and nodes using LS')
        tic = time.process_time()
        #
        self.EnrichedElems = self.libXFEM.getEnrichedElements(self.fluidNodes,
                                                              self.fluidElems,
                                                              self.structNodes,
                                                              self.structElems,
                                                              self.LevelSet)


        self.EdgeEnrichedElements = self.libXFEM.getEdgeEnrichedElements(self.fluidNodes,
                                                                        self.fluidElems,
                                                                        self.structNodes,
                                                                        self.structBounds)
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
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def analyseBC(self):
        """
        ##############################################################
        # Analyse BC to Build second member
        ##############################################################
        """
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Analyse boundary conditions')
        tic = time.process_time()
        # apply bc in displacement using data
        for bc in self.caseProp['bcdisp']:
            # get number of dofs and values to apply
            numDofs, valbc = self.getFormatBC(bc)
            bc['active'] = True
            bc['dofs'] = numDofs
            bc['load'] = valbc
            if len(numDofs) == 0:
                bc['active'] = False
                self.logger.warning('No nodes for BC in disp')
        # apply bc in pressure using data
        for bc in self.caseProp['bcpress']:
            # get number of dofs and values to apply
            numDofs, valbc = self.getFormatBC(bc)
            bc['active'] = True
            bc['dofs'] = numDofs
            bc['load'] = valbc
            # adpat list of internal fluid dofs
            self.SolvedDofF = np.setdiff1d(self.SolvedDofF, numDofs)
            if len(numDofs) == 0:
                bc['active'] = False
                self.logger.warning('No nodes for BC in press')
        #
        # fluid vector second member
        self.UF = np.zeros(2*self.fluidNbDofs, dtype=float)
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5),
                                 time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def buildSecondMember(self, omega):
        """
        ##############################################################
        # Build second member
        ##############################################################
        """
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Build second member')
        tic = time.process_time()
        # reset second member
        self.UF[:] = 0.0
        # apply bc in displacement using data
        for bc in self.caseProp['bcdisp']:
            if bc['active']:
                self.UF[bc['dofs']] = (omega**2)*bc['load']
        # apply bc in pressure using data
        for bc in self.caseProp['bcpress']:
            if bc['active']:
                self.UF[self.SolvedDofF] += -(self.KFF[self.SolvedDofF, :][:, bc['dofs']]
                                              -(omega**2)* self.MFF[self.SolvedDofF, :][:, bc['dofs']])*bc['load']
        #
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5),
                                 time.process_time()-tic))

        if self.paraData['gradCompute']:
            self.logger.info('Build gradient of second member')
            tic = time.process_time()
            for it in self.getNameGrad():
                self.dUF.append(None)
            self.logger.info('{} Done - {} s'
                             .format(self.deco.outPattern(lt=5),
                                     time.process_time()-tic))

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
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Build fluid operators')
        tic = time.process_time()
        # build matrices using vector description
        IIf, JJf, Vffk, Vffm = self.libXFEM.getGlobalAcousticMatrices(
            self.fluidNodes,
            self.fluidElems,
            self.mechaProp['celerity'],
            self.mechaProp['rho'])

        self.KFF = scipy.sparse.csc_matrix(
            (Vffk, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        self.MFF = scipy.sparse.csc_matrix(
            (Vffm, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))

        self.SolvedDofF = list(range(self.fluidNbDofs))
  
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))
        #
        if self.paraData['gradCompute']:
            self.logger.info('Build fluid gradient operators')
            tic = time.process_time()
            for it in self.getNameGrad():
                self.dKFF.append(None)
                self.dMFF.append(None)
            self.logger.info('{} Done - {} s'
                             .format(self.deco.outPattern(lt=5),
                                     time.process_time()-tic))

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
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Build enriched operators')
        tic = time.process_time()
        # build matrices using vector description
        IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = self.libXFEM.getEnrichedMatrices(
            self.fluidNodes,
            self.fluidElems,
            self.LevelSet,
            self.LevelSetTangent,
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
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))

        if self.paraData['gradCompute']:
            self.logger.info('Build gradient of enriched operators')
            tic = time.process_time()
            for it in range(self.getNbGrad()):
                # compute gradients matrices for each selected parameter
                dKAAx, dKFAx, dMAAx, dMFAx = self.buildGOperators(it)
                # store
                self.dKAA.append(dKAAx)
                self.dKFA.append(dKFAx)
                self.dMAA.append(dMAAx)
                self.dMFA.append(dMFAx)
            self.logger.info('{} Done - {} s'
                             .format(self.deco.outPattern(lt=5),
                                     time.process_time()-tic))

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
        self.logger.info(self.deco.equalPattern())
        self.logger.info('Assembly of operators of the whole problem')
        tic = time.process_time()
        #
        fd = self.mechaProp['fluid_damping']
        #
        self.K = scipy.sparse.construct.bmat([
            [fd*self.KFF[self.SolvedDofF, :][:, self.SolvedDofF],
                fd*self.KAF[self.SolvedDofF, :][:, self.SolvedDofA]],
            [fd*self.KAF[self.SolvedDofA, :][:, self.SolvedDofF], fd*self.KAA[self.SolvedDofA, :][:, self.SolvedDofA]]])
        #
        self.M = scipy.sparse.construct.bmat([
            [self.MFF[self.SolvedDofF, :][:, self.SolvedDofF],
                self.MAF[self.SolvedDofF, :][:, self.SolvedDofA]],
            [self.MAF[self.SolvedDofA, :][:, self.SolvedDofF], self.MAA[self.SolvedDofA, :][:, self.SolvedDofA]]])
        #
        # list of all dofs
        self.SolvedDof = np.hstack(
            [self.SolvedDofF, self.SolvedDofA+self.fluidNbDofs])
        #
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))
        #
        if self.paraData['gradCompute']:
            self.logger.info('Build gradient of the assembled operators')
            tic = time.process_time()
            #build full stiffness and mass matrices
            self.dK = list()
            self.dM = list()
            for it in range(self.getNbGrad()):
                # build submatrices
                K11 = None
                K12 = None
                K21 = None
                K22 = None
                M11 = None
                M12 = None
                M21 = None
                M22 = None
                #
                if self.dKFF[it] is not None:
                    K11 = fd*self.dKFF[it][self.SolvedDofF,
                                           :][:, self.SolvedDofF]
                if self.dKFA[it] is not None:
                    K12 = fd*self.dKFA[it][self.SolvedDofF,
                                           :][:, self.SolvedDofA]
                    K21 = fd*self.dKFA[it][self.SolvedDofA,
                                           :][:, self.SolvedDofF]
                if self.dKAA[it] is not None:
                    K22 = fd*self.dKAA[it][self.SolvedDofA,
                                           :][:, self.SolvedDofA]
                #
                if self.dMFF[it] is not None:
                    M11 = self.dMFF[it][self.SolvedDofF, :][:, self.SolvedDofF]
                if self.dMFA[it] is not None:
                    M12 = self.dMFA[it][self.SolvedDofF, :][:, self.SolvedDofA]
                    M21 = self.dMFA[it][self.SolvedDofA, :][:, self.SolvedDofF]
                if self.dMAA[it] is not None:
                    M22 = self.dMAA[it][self.SolvedDofA, :][:, self.SolvedDofA]
                # build full stiffness and mass gradient matrices
                self.dK.append(scipy.sparse.construct.bmat([
                    [K11, K12],
                    [K21, K22]]))
                self.dM.append(scipy.sparse.construct.bmat([
                    [M11, M12],
                    [M21, M22]]))
            #
            self.logger.info('{} Done - {} s'
                             .format(self.deco.outPattern(lt=5),
                                     time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def buildGOperators(self, itG):
        """
        ##################################################################
        # Build gradients operators one by one
        ##################################################################
        """
        self.logger.info(self.deco.equalPattern())
        self.logger.info(' Build gradient matrices for parameter {}'.format(
            self.paraData['nameGrad'][itG]))
        # compute gradients matrices and indices
        # print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)
        dataMat =\
            self.libXFEM.getGlobalAcousticGradientMatrices(
                self.fluidNodes,
                self.fluidElems,
                self.LevelSet,
                self.LevelSetTangent,
                self.LevelSetGradient[itG],
                self.mechaProp['celerity'],
                self.mechaProp['rho'])
        # manage obtained data
        IIf = dataMat[0]
        JJf = dataMat[1]
        Vfak_gradient = dataMat[2]
        Vfam_gradient = dataMat[3]
        # build sparse matrices
        dKAA_dtheta = None
        dMAA_dtheta = None
        dKFA_dtheta = scipy.sparse.csc_matrix(
            (Vfak_gradient, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        dMFA_dtheta = scipy.sparse.csc_matrix(
            (Vfam_gradient, (IIf, JJf)),
            shape=(self.fluidNbDofs, self.fluidNbDofs))
        #
        return dKAA_dtheta, dKFA_dtheta, dMAA_dtheta, dMFA_dtheta

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
        self.logger.info(self.deco.dashPattern())
        self.logger.info(self.deco.dashPattern())
        self.logger.info(self.deco.dashPattern())
        self.logger.info(self.deco.adaptTxtCenter(
            'Start SILEX vibro-acoustics computation'))
        if self.paraData['gradCompute']:
            self.logger.info(self.deco.adaptTxtCenter(
                'with gradients computations)'))
        self.logger.info(self.deco.dashPattern())
        self.logger.info(self.deco.dashPattern())

        # check & create folders/files
        self.createDatabase()
        # load MPI info
        self.loadMPI()
        # load MUMPS
        self.loadMUMPS()
        # load fluid mesh
        self.loadMesh(typeData='nodesFluid',
                      filename=self.getDatafile('fluidmesh'))
        self.loadMesh(typeData='elemsFluid',
                      filename=self.getDatafile('fluidmesh'))
        # load control volume
        self.loadMesh(typeData='elemsControlFluid',
                      filename=self.getDatafile('fluidmesh'))
        # initialize the builder class for Finite Element
        if self.libXFEM is None:
            self.libXFEM = libXFEM.SILEXlibXFEM(typeElem=self.getElemFluid())
        # build fluid operators
        self.buildFluidOperators()
        # analyse boundary conditions
        self.analyseBC()
        # generate the list of frequencies
        self.Frequencies = utils.computeFreqPerProc(
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

    def initRun(self, paraVal=None):
        """
        ##################################################################
        # method used to initialize the data before a run associated to a set of parameters
        ##################################################################
        """
        # initialize data
        self.initDataSolve()
        #
        self.showDataParaVal()
        # build levelset
        self.buildLevelSet(
            typeLS=self.caseProp['typeLS'],
            typeGEO=self.caseProp['typeGEOstruct'],
            paraVal=paraVal)
        # build enriched part
        self.buildEnrichedPart()
        # build operators
        self.buildEnrichedOperators()
        self.buildOperators()

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def solvePbOneStep(self, itF, itMax, freq):
        """
        Method used to solve the problem for one frequency step
        """
        self.logger.info(self.deco.equalPattern())
        tic = time.process_time()
        #
        self.logger.info(
            'Rank: {} - Solve whole problem for frequency {} (step {}/{})'
            .format(self.rankMPI, freq, itF+1, itMax))
        # compute natural frequency
        omega = 2*np.pi*freq
        # Build full second member
        self.buildSecondMember(omega)
        F = np.array(self.UF[self.SolvedDof], dtype=self.loadType())
        # solve the whole problem on pressure
        ticB = time.process_time()
        sol = self.solveLinear(self.K-(omega**2)*self.M, F)
        self.logger.info(
            'Rank: {} - Solve linear problem - Done - {} s'
            .format(self.rankMPI, time.process_time()-ticB))
        # build the fields
        ticB = time.process_time()
        # uncorrected pressure field
        self.pressureUncorrect[self.SolvedDofF, itF] = sol[list(
            range(len(self.SolvedDofF)))].copy()
        # enrichment field
        self.pressureEnrichment[self.SolvedDofA, itF] = sol[list(
            range(len(self.SolvedDofF), 
                  len(self.SolvedDofF)+len(self.SolvedDofA)))].copy()
        # corrected pressure field
        self.pressure[:, itF] = self.pressureUncorrect[:, itF].copy()
        self.pressure[self.SolvedDofA, itF] += self.pressureEnrichment[self.SolvedDofA,
                                                                       itF]*np.sign(self.LevelSet[self.SolvedDofA])
        self.logger.info(
            'Rank: {} - Fields computation - Done - {} s'
            .format(self.rankMPI, time.process_time()-ticB))
        # compute and store FRF on the control volume
        if self.caseProp['computeFRF']:
            self.FRF[itF] = self.libXFEM.computeQuadraticPressure(
                self.fluidNodes,
                self.fluidElemsControl,
                self.pressureUncorrect[:, itF],
                self.pressureEnrichment[:, itF],
                self.LevelSet)
        # Compute gradients of field
        if self.paraData['gradCompute']:
            # initialize variables

            # along the parameters
            self.logger.info(
                'Rank: {} - Start gradients computation'.format(self.rankMPI))
            ticC = time.process_time()
            for itG in range(self.getNbGrad()):
                ticB = time.process_time()
                # prepare data
                tmp = -(self.dK[itG]-(omega**2)*self.dM[itG])*sol
                Dsol = self.solveLinear(self.K-(omega**2)*self.M, tmp)
                self.logger.info('Rank: {} - Prepare and solve linear problem for grad (var: {}/num: {}) - Done - {} s'
                                 .format(self.rankMPI, self.getNameGrad()[itG], itG, time.process_time()-ticB))
                # gradient of uncorrected pressure field
                self.pressureUncorrectGrad[itG][self.SolvedDofF, itF] = Dsol[list(
                    range(len(self.SolvedDofF)))].copy()
                # gradient of enrichment field
                self.pressureEnrichmentGrad[itG][self.SolvedDofA, itF] = Dsol[list(
                    range(len(self.SolvedDofF), len(self.SolvedDofF)+len(self.SolvedDofA)))].copy()
                # gradient of corrected pressure field
                self.pressureGrad[itG][:,itF] = self.pressureUncorrectGrad[itG][:, itF].copy()
                self.pressureGrad[itG][self.SolvedDofA, itF] += self.pressureEnrichmentGrad[itG][self.SolvedDofA,
                                                                                                 itF]*np.sign(self.LevelSet[self.SolvedDofA].T)
                # compute and store gradient of FRF on the control volume
                if self.caseProp['computeFRF']:
                    self.FRFgrad[itG][itF] = self.libXFEM.computeGradQuadraticPressure(
                        nodesCavity=self.fluidNodes,
                        elemsControlVolume=self.fluidElemsControl,
                        Field=self.pressure[:, itF],
                        gradField=self.pressureGrad[itG][:, itF],
                        levelset=self.LevelSet)
            self.logger.info('Rank: {} - Gradients computation - Done - {} s'
                             .format(self.rankMPI, time.process_time()-ticC))

        self.logger.info('Rank: {} - Done - Solve whole problem for frequency {} (freq step {}/{}) in {} s'
                         .format(self.rankMPI, freq, itF+1, itMax, time.process_time()-tic))
        self.logger.info(self.deco.equalPattern())

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
        tic = time.process_time()
        self.logger.info('Initialize storage')
        self.pressureUncorrect = np.zeros(
            [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType())
        self.pressureEnrichment = np.zeros(
            [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType())
        self.pressure = np.zeros(
            [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType())
        #
        if self.caseProp['computeFRF']:
            self.FRF = np.zeros(self.caseProp['nbSteps'])
        #
        if self.paraData['gradCompute']:
            self.pressureUncorrectGrad = [np.zeros(
                [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureEnrichmentGrad = [np.zeros(
                [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureGrad = [np.zeros(
                [self.fluidNbDofs, self.caseProp['nbSteps']], dtype=self.loadType()) for _ in range(self.getNbGrad())]
            #
            if self.caseProp['computeFRF']:
                self.FRFgrad = np.zeros(
                    [self.getNbGrad(), self.caseProp['nbSteps']])
        # initialize class for exporting meshes
        self.classSave = None
        self.logger.info('{} Done - {} s'
                         .format(self.deco.outPattern(lt=5), 
                                 time.process_time()-tic))


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def solvePb(self, paraVal=None):
        """
        ##################################################################
        # Method used to solve the whole problem
        ##################################################################
        """
        self.logger.info(self.deco.equalPattern())
        self.logger.info(self.deco.equalPattern())
        self.logger.info(self.deco.equalPattern())
        self.logger.info(self.deco.equalPattern())
        ticS = time.process_time()
        if self.nbRuns == 0:
            self.preProcessMaster()
        if paraVal is not None:
            paraVal = np.array(paraVal)
        else:
            self.logger.error(
                '>> Unable to start computation due to no given parameters values')
            raise
        # prepare parameters values to run
        paraValOk = self.prepPara(paraVal)
        #
        # along the parameters
        for valU in paraValOk:
            self.nbRuns += 1
            ticV = time.process_time()
            self.logger.info(self.deco.dashPattern())
            self.paraData['val'] = valU
            txtPara = self.formatPara(valIn=valU)
            self.logger.info('Start compute for parameters (nb {}): {}'
                             .format(self.nbRuns, txtPara))
            # initialization for run
            self.initRun(paraVal=valU)

            # along the frequencies
            for (itF, Freq) in enumerate(self.getFrequencies()):
                # solve the problem
                self.solvePbOneStep(
                    itF, len(self.getFrequencies(total=True)), Freq)
            # export results (FRF, pressure fields and gradients)
            if self.flags['saveResults']:
                self.exportFieldsOnePara(paraName=True)
                self.saveFRF(paraName=True)
            # append data to history
            self.paraData['oldval'].append(valU)
            self.allFRF.append(self.FRF)
            self.allFRFgrad.append(self.FRFgrad)
            #
            self.logger.info('Time to solve the whole problem for set of parameters nb {} - {} s'
                             .format(self.nbRuns, time.process_time()-ticV))
        # plot and export results (FRF, pressure fields and gradients)
        if self.flags['saveResults']:
            self.saveFRF(paraName=False, allData=True)
            self.plotFRF(allData=True, fileOut=self.getResultFile(
                detPara=False, addTxt='allData', ext='csv'))

        #
        self.logger.info('Time to solve the whole problem along sets of parameters - {} s'
                         .format(time.process_time()-ticS))
        self.logger.info(self.deco.plusPattern())
        self.logger.info(self.deco.plusPattern())
        self.logger.info(self.deco.plusPattern())
        self.logger.info(self.deco.plusPattern())
