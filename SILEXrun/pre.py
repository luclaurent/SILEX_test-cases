###########################################################
# Declaration of useful pre-processing tools for SILEXclass
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################
import os
import numpy as np 
import logging
from datetime import datetime
import time
#
import structTools
import utils
#
from SILEXlib.mesh import msh

# activate logger
Logger = logging.getLogger(__name__)
class preProcess:


    def dataOk(self,dispFlag=None):
        """
        ##################################################################
        # method used to check if the data are available
        ##################################################################
        """
        statusData=True
        for key in ['geomFolder','resultsFolder','originalFluidMeshFile','currentStructMeshFile','resultsFile']:
            if not self.data[key]:
                statusData=False
        if dispFlag:
            Logger.info("Geometry folder: %s"%utils.prepareStr(self.data['geomFolder']))
            Logger.info("Results folder: %s"%utils.prepareStr(self.data['resultsFolder']))
            Logger.info("Original fluid mesh file: %s"%utils.prepareStr(self.data['originalFluidMeshFile']))
            Logger.info("Original structure mesh file: %s"%utils.prepareStr(self.data['originalStructMeshFile']))
            Logger.info("Current structure mesh file: %s"%utils.prepareStr(self.data['currentStructMeshFile']))
            Logger.info("Result mesh file: %s"%utils.prepareStr(self.data['resultsFile']))        
        return statusData

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def createDatabase(self):
        """
        ##################################################################
        # method used to create full path of folders
        ##################################################################
        """
        Logger.info('================================')
        #name of the symlink
        symlinkname='last'
        #create a specific folder for the results
        baseFolderName = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        if self.caseProp['name']:
            baseFolderName += '_'+self.caseProp['name']
        else:
            baseFolderName += '_'+self.data['originalFluidMeshFile']
        #initialize folders and files (could change along the runs)
        self.fullPathCurrentResultsFolder = os.path.join(self.data['resultsFolder'],baseFolderName)
        #create directory if not exists
        if not os.path.exists(self.fullPathCurrentResultsFolder):
            os.makedirs(self.fullPathCurrentResultsFolder)
        #remove symlink if exists
        fullpathSymLink = os.path.join(self.data['resultsFolder'],symlinkname)
        if os.path.islink(fullpathSymLink):
            os.unlink(fullpathSymLink)
        # create symlink to the folder
        os.symlink(os.path.relpath(self.fullPathCurrentResultsFolder,self.data['resultsFolder']),fullpathSymLink)

        #display
        Logger.info('Folder for results: %s'%(self.fullPathCurrentResultsFolder))
        Logger.info('================================')

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def loadMesh(self,typeData=None,dispFlag=True,filename=None):
        """
        ##################################################################
        # method used to load mesh files
        ##################################################################
        """
        Logger.info('================================')
        #dictionary of kind of data
        textDict = dict()
        textDict['nodesFluid'] = 'nodes of the fluid'
        textDict['elemsControlFluid'] = 'elements of the control volume in the fluid'
        textDict['elemsFluid'] = 'elements in the fluid'
        textDict['nodesStruct'] = 'nodes of the structure'
        textDict['elemsStruct'] = 'elements in the structure'
        #
        if filename is None:
            Logger.error("Filename of the mesh is missing")            

        if filename is not None:
            #check if file exist
            if utils.checkFile(filename,'file'):
                Logger.info('>> Read %s'.format(filename))
                tic = time.process_time()
                self.meshData = msh.mshReader(filename = filename, dim = self.getDim())
                Logger.info('++++++++++++++++++++ Done - {} s'.format(time.process_time()-tic))
                # deal with types
                if typeData=='nodesFluid':
                    self.fluidNodes = self.meshData.getnodes()
                    self.fluidNbNodes=self.fluidNodes.shape[0]
                    self.fluidNbDofs=self.fluidNbNodes
                    #
                    if dispFlag:
                        Logger.info("Fluid: {} nodes ({} dofs)".format(self.fluidNbNodes,self.fluidNbDofs))
                if typeData=='elemsControlFluid':
                    self.fluidElemsControl, self.idFluidNodesControl = self.meshData.getElements(tag = 5)   # air, ONLY controlled volume
                    self.fluidNbElemsControl=self.fluidElemsControl.shape[0]
                    self.fluidNbNodesControl=self.idFluidNodesControl.shape[0]
                    #
                    if dispFlag:
                        Logger.info("Fluid control volume: {} elems, {}} nodes".format(self.fluidNbElemsControl,self.fluidNbNodesControl))
                if typeData=='elemsFluid':
                    self.fluidElems, self.idFluidNodes = self.meshData.getElements(tag = 1)  # air, cavity + controlled volume
                    self.fluidNbElems=self.fluidElems.shape[0]
                    self.fluidNbNodes=self.idFluidNodes.shape[0]
                    #
                    if dispFlag:
                        Logger.info("Fluid whole volume: {} elems, {} nodes".format(self.fluidNbElems,self.fluidNbNodes))
                if typeData=='nodesStruct':
                    self.structNodes = msh.mshReader(filename = filename, dim = self.getDim())
                    self.structNbNodes = self.structNodes.shape[0]    
                    #
                    if dispFlag:
                        Logger.info("Structure: {} nodes".format(self.structNbNodes))
                if typeData=='elemsStruct':
                    self.structElems, self.idStructNodes = self.meshData.getElements(tag = 2)
                    self.structNbElems = self.structElems.shape[0]
                    self.structNbNodes = self.idStructNodes.shape[0]
                    #
                    if dispFlag:
                        Logger.info("Structure: {} elems, {} nodes".format(self.structNbElems,self.structNbNodes))
            else:
                Logger.error('Unable to read data: {}'.format(textDict[typeData]))
                Logger.error('>>> Unable to read file \'{}\' (file does not exist)'.format(filename))
                raise
        Logger.info('================================')


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def loadPara(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load parameters data with dictionary
        ##################################################################
        """
        Logger.info('================================')
        if dataIn is not None:
            Logger.info('>>> Load parameters properties <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.paraData[key]=dataIn[key]
                    Logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            Logger.info('>>> Available parameters properties and current values <<<')
            for key in self.paraData:
                Logger.info('>>>> {}: {}'.format(key,self.paraData[key]))
        Logger.info('================================')

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    
    def loadParaU(self,
                  namePara=None,
                  nameParaGrad=None,
                  valPara=None,
                  gradCompute=None,
                  nameGrad=None,
                  forceGradName = False):
        """
        ##################################################################
        # method used to load parameters data (require for gradient(s) computation(s))
        ##################################################################
        """
        if namePara is not None:
            self.paraData['name']=namePara
        if nameParaGrad is not None:
            self.paraData['nameGrad']=nameParaGrad
        if valPara is not None:
            self.paraData['val']=valPara
        if gradCompute is not None:
            self.paraData['gradCompute']=gradCompute
        if nameGrad is not None:
            self.paraData['nameGrad']=nameGrad
        #

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

        
    def loadMechaProperties(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load mechanical properties
        ##################################################################
        """
        Logger.info('================================')
        if dataIn is not None:
            Logger.info('>>> Load Mechanical properties <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.mechaProp[key]=dataIn[key]
                    Logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            Logger.info('>>> Available Mechanical properties and current values <<<')
            for key in self.mechaProp:
                Logger.info('>>>> {}: {}'.format(key,self.mechaProp[key]))
        Logger.info('================================')

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def loadData(self,dataIn=None,force=False):
        """
        ##################################################################
        # load data for the case (mesh file, directories,...)
        ##################################################################
        """
        Logger.info('================================')
        if dataIn is not None:
            Logger.info('>>> Load data  <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.data[key]=dataIn[key]
                    Logger.info('>> %s: %s'%(key,dataIn[key]))
        else:
            Logger.info('>>> Available data and current values <<<')
            for key in self.data:
                Logger.info('>>>> %s: %s'%(key,self.data[key]))
        Logger.info('================================')

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    
    def loadComputeProperties(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load computation properties
        ##################################################################
        """
        Logger.info('================================')
        if dataIn is not None:
            Logger.info('>>> Load properties for computation <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.caseProp[key]=dataIn[key]
                    Logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            Logger.info('>>> Available properties for computation and current values <<<')
            for key in self.caseProp:
                Logger.info('>>>> {}: {}'.format(key,self.mechcasePropaProp[key]))
        # try to prepare the case
        self.prepCase()
        Logger.info('================================')
        

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def loadBC(self,dataIn=None):
        """
        ##################################################################
        # method use to declare boundary conditions
        ##################################################################
        """
        if dataIn is not None:
            Logger.info('================================')
            Logger.info('>>> Load boundary condition(s) <<<')
            for key in dataIn:
                if key=='disp':
                    self.caseProp['bcdisp'].append(dataIn[key])
                    Logger.info('>> Add new bc: type displacement (nb {})'.format(len(self.caseProp['bcdisp']))
                if key=='press':
                    self.caseProp['bcpress'].append(dataIn[key])
                    Logger.info('>> Add new bc: type pressure (nb {})'.format(self.caseProp['bcpress']))
            Logger.info('================================')

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    
    def getFormatBC(self,bcIn):
        """
        ##############################################################
        # Prepare BC from data
        ##############################################################
        # bcIn list of dictionaries syntax:
        # - type: nodes, bbx (key)
        # - for nodes: value(s) is a list of nodes
        # - for bbx: value(s) are a declaration of boundary box (6 coordinates in 3D)
        # - value(s): values of the bc (many values could be given in the case of nodes (1 per nodes))
        # ex: {'type':'nodes','data':[12,3,14],'values':[xx,yy,zz]} or
        # ex: {'type':'bbx','data':[[xm,XM],[ym,YM],[zm,ZM]],'values':[uu]}
        # notice that 'bbx' could be also declared using 'bbx':[xm,XM,ym,YM,zm,ZM]
        """
        nodesList = []
        valuesOut = []
        nbNodes = 0
        if 'bbx' in bcIn['type']:
            nodesList=utils.getNodesBBX(self.fluidNodes,bcIn['data'])
            nbNodes=len(nodesList)
        if 'nodes' in bcIn['type']:
            # number of given nodes
            nbNodes=len(bcIn['data'])
            nodesList=bcIn['type']
        #deal with value(s) of bc: two cases
        val=np.array(bcIn['values'])
        if len(val.shape)==0 or len(val.shape)==1:
            valuesOut=np.zeros([nbNodes])
            valuesOut[:]=val
        elif len(bcIn['values'])==nbNodes:
            valuesOut=val
        else:
            Logger.error('Bad statement of boundary condition')
        #
        return nodesList,valuesOut

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def prepCase(self):
        """
        ##################################################################
        # method to initialize case's data
        ##################################################################
        """
        if self.caseProp['typeLS'] == 'manual':
            #load the object to build LS
            LSobj=structTools.LSmanual(self.caseProp['typeGEOstruct'])
            #export values
            self.loadParaU(namePara=LSobj.exportParaName())

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
