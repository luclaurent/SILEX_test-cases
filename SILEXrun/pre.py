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
from . import structTools
from .misc import utils
#
from meshRW import msh as IOmsh

# # activate logger
# Logger = logging.getLogger(__name__)

class preProcess(object):
    def dataOk(self,dispFlag=True):
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
            self.logger.info("Geometry folder: {}".format(utils.prepareStr(self.data['geomFolder'])))
            self.logger.info("Results folder: {}".format(utils.prepareStr(self.data['resultsFolder'])))
            self.logger.info("Original fluid mesh file: {}".format(utils.prepareStr(self.data['originalFluidMeshFile'])))
            self.logger.info("Original structure mesh file: {}".format(utils.prepareStr(self.data['originalStructMeshFile'])))
            self.logger.info("Current structure mesh file: {}".format(utils.prepareStr(self.data['currentStructMeshFile'])))
            self.logger.info("Result mesh file: {}".format(utils.prepareStr(self.data['resultsFile'])))
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
        self.logger.info(self.deco.equalPattern())
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
        self.logger.info('Folder for results: {}'.format(self.fullPathCurrentResultsFolder))
        self.logger.info(self.deco.equalPattern())

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def loadMesh(self,
                typeData = None,
                filename = None,
                force = False):
        """
        ##################################################################
        # method used to load mesh files
        ##################################################################
        """
        currMeshData = None # current mesh data
        dataOk = False      # flag for availability of the data
        #
        self.logger.info(self.deco.equalPattern())
        #dictionary of kind of data
        textDict = dict()
        textDict['nodesFluid'] = 'nodes of the fluid'
        textDict['elemsControlFluid'] = 'elements of the control volume in the fluid'
        textDict['elemsFluid'] = 'elements in the fluid'
        textDict['nodesStruct'] = 'nodes of the structure'
        textDict['elemsStruct'] = 'elements in the structure'
        #
        if filename is None and force:
            self.logger.error("Filename of the mesh is missing")            
        # type of area
        currArea = None
        if 'Fluid' in typeData:
            currArea = 'Fluid'
        if 'Struct' in typeData:
            currArea = 'Struct'
        # get the data
        if currArea in self.meshData.keys():
            currMeshData = self.meshData[currArea]
            dataOk = True
        #read file
        if filename is not None and (not dataOk or force):
            #check if file exist
            if utils.checkFile(filename,'file'):
                self.logger.info('>> Read {}'.format(filename))
                tic = time.process_time()
                self.meshData[currArea] = IOmsh.mshReader(filename = filename, dim = self.getDim())
                currMeshData = self.meshData[currArea]
                self.logger.info('++++++++++++++++++++ Done - {} s'.format(time.process_time()-tic))
            else:
                self.logger.error('>>> Unable to read file \'{}\' (file does not exist)'.format(filename))
                self.logger.error('Unable to read data: {}'.format(textDict[typeData]))                

        #store data
        # deal with types
        if typeData=='nodesFluid':
            self.fluidNodes = currMeshData.getNodes()
            self.fluidNbNodes = self.fluidNodes.shape[0]
            self.fluidNbDofs = self.fluidNbNodes
            #
            self.logger.debug("Fluid: {} nodes ({} dofs)".format(self.fluidNbNodes,self.fluidNbDofs))
        if typeData=='elemsControlFluid':
            self.fluidElemsControl, self.idFluidNodesControl = currMeshData.getElements(tag = 5,dictFormat = False)   # air, ONLY controlled volume
            self.fluidNbElemsControl = self.fluidElemsControl.shape[0]
            self.fluidNbNodesControl = self.idFluidNodesControl.shape[0]
            #
            self.logger.debug("Fluid control volume: {} elems, {} nodes".format(self.fluidNbElemsControl,self.fluidNbNodesControl))
        if typeData=='elemsFluid':
            self.fluidElems, self.idFluidNodes = currMeshData.getElements(tag = 1,dictFormat = False)  # air, cavity + controlled volume
            self.fluidNbElems = self.fluidElems.shape[0]
            self.fluidNbNodes = self.idFluidNodes.shape[0]
            #
            self.logger.debug("Fluid whole volume: {} elems, {} nodes".format(self.fluidNbElems,self.fluidNbNodes))
        if typeData=='nodesStruct':
            self.structNodes = currMeshData.getNodes()
            self.structNbNodes = self.structNodes.shape[0]
            #
            self.logger.debug("Structure: {} nodes".format(self.structNbNodes))
        if typeData=='elemsStruct':
            self.structElems, self.idStructNodes = currMeshData.getElements(tag = 2,dictFormat = False)
            self.structNbElems = self.structElems.shape[0]
            self.structNbNodes = self.idStructNodes.shape[0]
            #
            self.logger.debug("Structure: {} elems, {} nodes".format(self.structNbElems,self.structNbNodes))

        self.logger.info(self.deco.equalPattern())


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
        self.logger.info(self.deco.equalPattern())
        if dataIn is not None:
            self.logger.info(self.deco.adaptTxtCenter('Load parameters properties'))
            for key in dataIn:
                if dataIn[key] or force:
                    self.paraData[key]=dataIn[key]
                    self.logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            self.logger.info(self.deco.adaptTxtCenter('Available parameters properties and current values'))
            for key in self.paraData:
                self.logger.info('>>>> {}: {}'.format(key,self.paraData[key]))
        self.logger.info(self.deco.equalPattern())

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
            self.paraData['name'] = namePara
        if nameParaGrad is not None:
            self.paraData['nameGrad'] = nameParaGrad
        if valPara is not None:
            self.paraData['val'] = valPara
        if gradCompute is not None:
            self.paraData['gradCompute'] = gradCompute
        if nameGrad is not None:
            self.paraData['nameGrad'] = nameGrad
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
        self.logger.info(self.deco.equalPattern())
        if dataIn is not None:
            self.logger.info(self.deco.adaptTxtCenter('Load Mechanical properties'))
            for key in dataIn:
                self.mechaProp[key]=dataIn[key]
                self.logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            self.logger.info(self.deco.adaptTxtCenter('Available Mechanical properties and current values'))
            for key in self.mechaProp:
                self.logger.info('>>>> {}: {}'.format(key,self.mechaProp[key]))
        self.logger.info(self.deco.equalPattern())

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
        self.logger.info(self.deco.equalPattern())
        if dataIn is not None:
            self.logger.info(self.deco.adaptTxtCenter('Load data'))
            for key in dataIn:
                if dataIn[key] or force:
                    self.data[key]=dataIn[key]
                    self.logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            self.logger.info(self.deco.adaptTxtCenter('Available data and current values'))
            for key in self.data:
                self.logger.info('>>>> {}: {}'.format(key,self.data[key]))
        self.logger.info(self.deco.equalPattern())

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
        self.logger.info(self.deco.equalPattern())
        if dataIn is not None:
            self.logger.info(self.deco.adaptTxtCenter('Load properties for computation'))
            for key in dataIn:
                if dataIn[key] or force:
                    self.caseProp[key]=dataIn[key]
                    self.logger.info('>> {}: {}'.format(key,dataIn[key]))
        else:
            self.logger.info(self.deco.adaptTxtCenter('Available properties for computation and current values'))
            for key in self.caseProp:
                self.logger.info('>>>> {}: {}'.format(key,self.mechcasePropaProp[key]))
        # try to prepare the case
        self.prepCase()
        self.logger.info(self.deco.equalPattern())
        

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
            self.logger.info(self.deco.equalPattern())
            self.logger.info(self.deco.adaptTxtCenter('Load boundary condition(s)'))
            for key in dataIn:
                if key == 'disp':
                    self.caseProp['bcdisp'].append(dataIn[key])
                    self.logger.info('>> Add new bc: type displacement (nb {})'.format(len(self.caseProp['bcdisp'])))
                elif key == 'press':
                    self.caseProp['bcpress'].append(dataIn[key])
                    self.logger.info('>> Add new bc: type pressure (nb {})'.format(self.caseProp['bcpress']))
                else:
                    self.logger.warning('>> {} bc not accepted'.format(key))
            self.logger.info(self.deco.equalPattern())

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
            self.logger.error('Bad statement of boundary condition')
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
