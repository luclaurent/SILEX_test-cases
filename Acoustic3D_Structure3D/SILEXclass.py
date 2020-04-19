###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

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

###########################################################
# Libraries
###########################################################
import getopt
import importlib
import string
import time
import logging
from datetime import datetime 
#
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.io
#
import pickle
#
import sys
import os
from shutil import copyfile
#
import utils
import structTools

from SILEX import silex_lib_xfem_acou_tet4
from SILEX import silex_lib_gmsh

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# function to obtain info concerning MPI
def mpiInfo():
    # default values
    nproc=1
    rank=0
    comm=None
    #try to import MPI    
    mpi4py_loader=importlib.util.find_spec('mpi4py')
    if mpi4py_loader is not None:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nproc = comm.Get_size()
        rank = comm.Get_rank()
    else:
        print('No mpi4py.MPI module found')
        pass    
    return nproc, rank, comm

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# function to load MUMPS library
def loadMumps():
    mumps_loader=importlib.util.find_spec('mumps')
    foundMumps=mumps_loader is not None
    if foundMumps:
        globals()['mumps'] = __import__('mumps')
    return foundMumps




###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# computation class
class SILEX:
    #case properties
    caseProp = dict()
    caseProp['freqMax'] = []          # maximum frequency of the range
    caseProp['freqMin'] = []          # minimum frequency of the range
    caseProp['nbSteps'] = []          # number of frequency steps
    caseProp['modal'] = False         # modal version of the computation (building of modal basis and use it for gradients)
    caseProp['computeFRF'] = True     # computation of the FRF in control volume
    #
    caseProp['bcdisp'] = list()       # boundary conditions on displacement
    caseProp['bcpress'] = list()      # boundary conditions on pressure
    #
    caseProp['name'] = ''             # name of the case
    #
    caseProp['typeLS'] = ''           # type of Level-Set (FromMesh or manual)
    caseProp['typeGEOstruct'] = ''    # type of geometry of the structure (in the case of manual declaration (see structTools.py))

    #parameters values
    paraData = dict()
    paraData['oldval'] = list()       # previous values of parameters
    paraData['val'] = None            # current values of parameters
    paraData['name'] = None           # name of parameters
    paraData['nb'] = None             # number of parameters
    paraData['nameGrad'] = None       # name of parameters for gradients
    paraData['nbGrad'] = None         # number of gradients
    paraData['gradCompute'] = False   # compute gradients or not

    #material properties
    # air
    mechaProp = dict()
    mechaProp['celerity'] = 340.0
    mechaProp['rho'] = 1.2
    mechaProp['fluid_damping'] = (1+0.01j)

    #data properties
    data = dict()
    data['geomFolder'] = 'geom'             # folder of geometry and meshes
    data['resultsFolder'] = 'results'        # folder of results
    #
    data['originalFluidMeshFile'] = ''      # provided fluid mesh file
    data['originalStructGeoFile'] = ''      # provided structure geometry file (parametric, gmsh format...)
    data['originalStructMeshFile'] = ''     # provided structure mesh file
    data['currentStructMeshFile'] = ''      # build structure mesh file (generated from geometry)
    data['resultsFile'] = ''                # current results file basename
    data['prefixResults'] = 'SILEX'         # prefix use on results files
    #
    data['exportData'] = 'mat'              # default format to export data (as FRF, mat or pickle)
    data['exportMesh'] = 'msh'              # default format to export fields and mesh (msh, msh+, vtk)
    #
    fullPathCurrentResultsFolder = ''
    fullPathFluidMeshFile = ''
    fullPathStructMeshFile = ''
    fullPathMeshFile = ''
    fullPathResultsFile = ''
    #
    saveResults = True

    #architecture properties
    commMPI = None
    rankMPI = 0
    nbProcMPI = 0
    #
    mumpsOk = False

    #flags
    flags = dict()
    flags['saveResults'] = False  # flag to save results
    flags['edgeEnrichement'] = False  # flag to enrich the edge

    ###
    # storage variables
    LevelSet = []               # nodal values of LS (known at fluid nodes) signed
    LevelSetU = []              # nodal values of LS (known at fluid nodes) unsigned
    LevelSetGradient = list()   # nodal values of the parametric gradients of LS (known at fluid nodes)
    #
    fluidNodes = []             # coordinates of the fluid nodes
    fluidNbNodes = 0            # number of fluid nodes
    fluidNbDofs = 0             # number of dofs in fluid
    #
    fluidElems = []             # array of elements for the fluid volume
    idFluidNodes = []           # 
    fluidNbElems = 0            # number of elements in fluid volume
    fluidNbNodes = 0            # number of nodes in fluid volume
    #
    fluidElemsControl = []      # array of elements for the control volume
    idFluidNodesControl = []    # 
    fluidNbElemsControl = 0     # number of elements in control volume
    fluidNbNodesControl = 0     # number of nodes in control volume
    #
    structNodes = []            # coordinates of the structure nodes 
    structNbNodes = []          # number of structures nodes
    #
    structElems = []            # array of elements of structure
    idStructNodes = []          # 
    structNbElems = 0           # number of elements for structure
    structNbNodes = 0           # number of nodes of structure
    #enriched parts
    EnrichedNodes = []          # nodes associated to enriched elements
    EnrichedElems = []          # enriched elements
    EnrichedNbElems = 0        # number of enriched elements
    ##############
    #operators for fluid part
    KFF = []                    # fluid rigidity matrix
    MFF = []                    # fluid mass matrix
    dKFF = list()               # gradient matrices
    dMFF = list()               # gradient matrices
    SolvedDofF = []             # list of solved fluid dofs
    ##############
    #operators for enriched part
    KAA = []                    # rigidity matrix associated to full enriched dofs
    KFA = []                    # rigidity matrix associated to mixed enriched dofs
    MAA = []                    # mass matrix associated to full enriched dofs
    MFA = []                    # mass matrix associated to mixed enriched dofs
    dKAA = list()               # gradient matrices
    dKFA = list()               # gradient matrices
    dMAA = list()               # gradient matrices
    dMFA = list()               # gradient matrices
    SolvedDofA = []             # list of solved enriched dofs
    ##############
    #operators for the whole coupled problem
    K = []                      # rigidity matrix
    M = []                      # mass matrix
    dK = list()                 # gradient matrices
    dM = list()                 # gradients matrices
    UF = []                     # second member
    dUF = list()                # gradients of second member
    SolvedDof = []
    ##############
    #operators for results of the whole coupled problem
    pressure = []               # pressure field (corrected with enrichment field)
    pressureEnrichment = []    # enrichment part of the pressure field
    pressureUncorrect = []      # uncorrected part of the pressure field
    #
    pressureGrad = list()               # gradient of the pressure field (corrected with enrichment field)
    pressureEnrichmentGrad = list()    # gradient of the enrichment part of the pressure field
    pressureUncorrectGrad = list()      # gradient of the uncorrected part of the pressure field

    ##############
    #operators for FRF
    Frequencies = []            # frequencies for computation
    FRF = []                    # results of frequency response function
    FRFgrad = list()            # gradients of the FRF

    ##############
    #others data
    nbRuns = 0                   # number of runs

    ##################################################################
    ##################################################################
    ##################################################################
    ##################################################################
    def __init__(self):
        """
        ##################################################################
        # Constructor of the class
        ##################################################################
        """        
        #initialize logging
        loggingFile = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")+"_SILEX.log"
        logging.basicConfig(
            handlers=[
                logging.FileHandler(loggingFile),
                logging.StreamHandler(sys.stdout)
                ],
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
            #
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##               Load SILEX object              ##")
        
    def loadType(self):
        """
        ##################################################################
        # Load the type of variable (float or complex)
        ##################################################################
        """
        if np.imag(self.mechaProp['fluid_damping'])==0:
            typeRun=np.dtype('float')
        else:
            typeRun=np.dtype('complex')
        return typeRun

    def getResultFile(self,detPara=False,addTxt=None,ext=None):
        """
        ##################################################################
        # Build result file depending on status
        ##################################################################
        """
        fileName = data['prefixResults']
        #add parameters value in filename
        if detPara:
            fileName += '_'+self.formatPara(nospace=True)+'_'
        #add given text in filename
        if addTxt is not None:
            if addTxt[0] != '_':
                fileName += '_' + addTxt
            else:
                fileName += addTxt
        #add extension in filename
        if ext is not None:
            fileName += ext
        #store the full path of the file
        if self.fullPathCurrentResultsFolder == '':
            folder = data['resultsFolder']
        else:
            folder = self.fullPathCurrentResultsFolder
        self.fullPathResultsFile = os.path.join(folder,fileName)
        #
        return self.fullPathResultsFile

    def getDatafile(self,typeFile=None):
        """
        ##################################################################
        # Get data file (mesh, geometry,...)
        ##################################################################
        """
        filepath = None
        if typeFile is not None:
            if typeFile == 'fluidmesh':
                filepath = os.path.join(self.data['geomFolder'],self.data['originalFluidMeshFile'])
                self.fullPathFluidMeshFile = filepath
            if typeFile == 'structmesh':
                filepath = os.path.join(self.data['geomFolder'],self.data['originalStructMeshFile'])
                self.fullPathStructMeshFile = filepath
            if typeFile == 'geo':
                filepath = os.path.join(self.data['geomFolder'],self.data['originalStructGeoFile'])
        else:
            logging.info('no type of file given')
        #
        return filepath

    def loadMPI(self):
        """
        ##################################################################
        # method used to load openMPI information
        ##################################################################
        """
        self.nbProcMPI,self.rankMPI,self.commMPI=mpiInfo()

    def loadMUMPS(self):
        """
        ##################################################################
        # method used to load openMPI information
        ##################################################################
        """
        self.mumpsOk=loadMumps()

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
            logging.info("Geometry folder: %s"%utils.prepareStr(self.data['geomFolder']))
            logging.info("Results folder: %s"%utils.prepareStr(self.data['resultsFolder']))
            logging.info("Original fluid mesh file: %s"%utils.prepareStr(self.data['originalFluidMeshFile']))
            logging.info("Original structure mesh file: %s"%utils.prepareStr(self.data['originalStructMeshFile']))
            logging.info("Current structure mesh file: %s"%utils.prepareStr(self.data['currentStructMeshFile']))
            logging.info("Result mesh file: %s"%utils.prepareStr(self.data['resultsFile']))        
        return statusData


    def createDatabase(self):
        """
        ##################################################################
        # method used to create full path of folders
        ##################################################################
        """
        #create a specific folder for the results
        baseFolderName = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        if self.caseProp['name']:
            baseFolderName += '_'+self.caseProp['name']
        else:
            baseFolderName += '_'+data['originalFluidMeshFile']
        #initialize folders and files name (could change along the runs)
        createOk=False
        if self.dataOk(True):
            self.fullPathCurrentResultsFolder = os.path.join(self.data['resultsFolder'],baseFolderName)
            self.fullPathOriginalMeshFile = os.path.join(self.data['geomFolder'],self.data['originalMeshFile'])
            self.fullPathMeshFile = os.path.join(self.data['geomFolder'],self.data['currentMeshFile'])
            self.fullPathResultsFile = os.path.join(self.data['resultsFolder'],self.data['resultsFile'])
            #create directory if not exists
            if not os.path.exists(self.fullPathCurrentResultsFolder):
                os.makedirs(self.fullPathCurrentResultsFolder)
            #display
            logging.info('Folder for results: %s'%(self.fullPathCurrentResultsFolder))
        else:
            logging.info('Missing data to create database')


    def loadMesh(self,typeData=None,dispFlag=True,filename=None):
        """
        ##################################################################
        # method used to load mesh files
        ##################################################################
        """
        #dictionary of kind of data
        textDict = dict()
        textDict['nodesFluid'] = 'nodes of the fluid'
        textDict['elemsControlFluid'] = 'elements of the control volume in the fluid'
        textDict['elemsFluid'] = 'elements in the fluid'
        textDict['nodesStruct'] = 'nodes of the structure'
        textDict['elemsStruct'] = 'elements in the structure'
        #
        if filename is None:
            logging.error("Filename of the mesh is missing")            

        if filename is not None:
            #check if file exist
            if utils.checkFile(filename,'file'):
                # deal with types
                if typeData=='nodesFluid':
                    logging.info('>> Read %s'%textDict[typeData])
                    tic = time.process_time()
                    self.fluidNodes = silex_lib_gmsh.ReadGmshNodes(filename, 3)
                    self.fluidNbNodes=self.fluidNodes.shape[0]
                    self.fluidNbDofs=self.fluidNbNodes
                    logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
                    #
                    if dispFlag:
                        logging.info("Fluid: %i nodes (%i dofs)"%(self.fluidNbNodes,self.fluidNbDofs))
                if typeData=='elemsControlFluid':
                    logging.info('>> Read %s'%textDict[typeData])
                    tic = time.process_time()
                    self.fluidElemsControl, self.idFluidNodesControl = silex_lib_gmsh.ReadGmshElements(filename, 4, 5)  # air, ONLY controlled volume
                    self.fluidNbElemsControl=self.fluidElemsControl.shape[0]
                    self.fluidNbNodesControl=self.idFluidNodesControl.shape[0]
                    logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
                    #
                    if dispFlag:
                        logging.info("Fluid control volume: %i elems, %i nodes"%(self.fluidNbElemsControl,self.fluidNbNodesControl))
                if typeData=='elemsFluid':
                    logging.info('>> Read %s'%textDict[typeData])
                    tic = time.process_time()
                    self.fluidElems, self.idFluidNodes = silex_lib_gmsh.ReadGmshElements(filename, 4, 1)  # air, cavity + controlled volume
                    self.fluidNbElems=self.fluidElems.shape[0]
                    self.fluidNbNodes=self.idFluidNodes.shape[0]
                    logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
                    #
                    if dispFlag:
                        logging.info("Fluid whole volume: %i elems, %i nodes"%(self.fluidNbElems,self.fluidNbNodes))
                if typeData=='nodesStruct':
                    logging.info('>>Read %s'%textDict[typeData])
                    tic = time.process_time()
                    self.structNodes = silex_lib_gmsh.ReadGmshNodes(filename, 3)
                    self.structNbNodes = self.structNodes.shape[0]    
                    logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
                    #
                    if dispFlag:
                        logging.info("Structure: %i nodes"%(self.structNbNodes))
                if typeData=='elemsStruct':
                    logging.info('>> Read %s'%textDict[typeData])
                    tic = time.process_time()
                    self.structElems, self.idStructNodes = silex_lib_gmsh.ReadGmshElements(filename, 2, 2)
                    self.structNbElems = self.structElems.shape[0]
                    self.structNbNodes = self.idStructNodes.shape[0]
                    logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))
                    #
                    if dispFlag:
                        logging.info("Structure: %i elems, %i nodes"%(self.structNbElems,self.structNbNodes))
            else:
                logging.error('Unable to read data: %s'%textDict[typeData])
                logging.error('>>> Unable to read file \'%s\' (file does not exist)'%filename)
                raise

    def getFrequencies(self,total=False):
        """
        ##################################################################
        # obtain the list of frequencies for the rank or the whole list of frequencies
        ##################################################################
        """
        if total:
            listFreq = self.Frequencies.T.flatten()
        else:
            listFreq = self.Frequencies[:,self.rankMPI]
        return listFreq

    def getNbGrad(self):
        """
        ##################################################################
        # obtain the number of computed gradients
        ##################################################################
        """
        nbG = self.paraData['nbGrad']
        if not nbG:
            nbG = len(self.paraData['name'])
        return nbG
        

    def getNameGrad(self):
        """
        ##################################################################
        # obtain the name of computed gradients parameters names
        ##################################################################
        """
        nbG = self.paraData['nameGrad']
        if not nbG:
            nbG = self.paraData['name']
        return nbG


    def loadPara(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load parameters data with dictionary
        ##################################################################
        """
        if dataIn is not None:
            logging.info('>>> Load parameters properties <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.paraData[key]=dataIn[key]
                    logging.info('>> %s: %s'%(key,dataIn[key]))
        else:
            logging.info('>>> Available parameters properties and current values <<<')
            for key in self.paraData:
                logging.info('>>>> %s: %s'%(key,paraData[key]))
    
    def loadParaU(self,namePara=None,nameParaGrad=None,valPara=None,gradCompute=None,nameGrad=None,forceGradName = False):
        """
        ##################################################################
        # method used to load parameters data (require for gradient(s) computation(s))
        ##################################################################
        """
        if namePara is not None:
            self.paraData['name']=namePara
            self.paraData['nb']=len(namePara)
        if nameParaGrad is not None:
            if self.paraData['nameGrad'] is None or forceGradName:
                self.paraData['nameGrad']=namePara
                self.paraData['nbGrad']=len(namePara)
        if valPara is not None:
            self.paraData['val']=valPara
        if gradCompute is not None:
            self.paraData['gradCompute']=gradCompute
        if nameGrad is not None:
            self.paraData['nameGrad']=nameGrad
        #
        


    def showDataParaVal(self):
        """
        ##################################################################
        # method used to show the data concerning the design parameters along 
        # the gradients will be computed
        ##################################################################
        """
        #check number of names vs number of parameter values
        pName=self.paraData['name']
        pVal=self.paraData['val']
        if len(pName)<len(pVal):
            logging.warning('Bad number of parameters names')
            itP=0
            while len(pName)<len(pVal):
                pName.append("Temp_%i"%itP)
                itP=+1
        #prepare showing gradients
        pGrad=['No'] * len(pVal)
        if self.paraData['gradCompute']:            
            if self.getNameGrad():
                for vv in self.getNameGrad():
                    pGrad[vv]='Yes'
            else:
                pGrad=['Yes'] * len(pVal)

        #show information
        logging.info('>>> Parameters values <<<')
        itPara=1
        for (n,p,g) in zip(pName,pVal,pGrad):
            logging.info('>>> #%i| name: %s - value: %s - grad: %s',itPara,n,p,g)
            itPara=+1

    
    def loadMechaProperties(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load mechanical properties
        ##################################################################
        """
        if dataIn is not None:
            logging.info('>>> Load Mechanical properties <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.mechaProp[key]=dataIn[key]
                    logging.info('>> %s: %s'%(key,dataIn[key]))
        else:
            logging.info('>>> Available Mechanical properties and current values <<<')
            for key in self.mechaProp:
                logging.info('>>>> %s: %s'%(key,mechaProp[key]))

    def loadData(self,dataIn=None,force=False):
        """
        ##################################################################
        # load data for the case (mesh file, directories,...)
        ##################################################################
        """
        if dataIn is not None:
            logging.info('>>> Load data  <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.data[key]=dataIn[key]
                    logging.info('>> %s: %s'%(key,dataIn[key]))
        else:
            logging.info('>>> Available data and current values <<<')
            for key in self.data:
                logging.info('>>>> %s: %s'%(key,mechaProp[key]))
    
    def loadComputeProperties(self,dataIn=None,force=False):
        """
        ##################################################################
        # method used to load computation properties
        ##################################################################
        """
        if dataIn is not None:
            logging.info('>>> Load properties for computation <<<')
            for key in dataIn:
                if dataIn[key] or force:
                    self.caseProp[key]=dataIn[key]
                    logging.info('>> %s: %s'%(key,dataIn[key]))
        else:
            logging.info('>>> Available properties for computation and current values <<<')
            for key in self.caseProp:
                logging.info('>>>> %s: %s'%(key,mechaProp[key]))
        # try to prepare the case
        self.prepCase()

    def loadBC(self,dataIn=None):
        """
        ##################################################################
        # method use to declare boundary conditions
        ##################################################################
        """
        if dataIn is not None:
            logging.info('>>> Load boundary condition(s) <<<')
            for key in dataIn:
                if key=='disp':
                    self.caseProp['bcdisp'].append(dataIn[key])
                    logging.info('>> Add new bc: type displacement (nb %i)'%len(self.caseProp['bcdisp']))
                if key=='press':
                    self.caseProp['bcpress'].append(dataIn[key])
                    logging.info('>> Add new bc: type pressure (nb %i)'%len(self.caseProp['bcpress']))
    
    def createResultsName(self,txt=None):
        """
        ##################################################################
        # generate basename of the results file
        ##################################################################
        """
        # number of parameters
        nbPara = len(self.paraData['val'])
        # prepare save file
        file_basis = "{:.{}E}".format(self.paraData['val'][0], 2)
        if (nbPara > 1):
            for i in range(1, nbPara):
                file_basis = file_basis+'_'+"{:.{}E}".format(self.paraData['val'][i], 2)

        results_file = data['resultsFile'] +'_'+file_basis
        #
        logging.debug('Results base name %s'%results_file)
        return results_file

    def exportResults(self,typeExport=None,method='msh',dictFields = None,fileName = None):
        """
        ##################################################################
        # function used to export results and mesh
        ##################################################################
        """
        if typeExport is not None:
            logging.info("Export results: type %s, method %s"%(typeExport,method))
            tic = time.process_time()                
            if typeExport is "manuFields":
                #export many fields with details
                if dictFields is None:
                    logging.info("Fields are missing")
                else:
                    if method is "msh":
                        if dictFields is dict:
                            dictFields=[dictFields]
                        #prepare data to export
                        dataExport= list()
                        for fields in dictFields:
                            dataOneField = list()
                            dataOneField.append(fields['field'])
                            if fields['type'] is "nodal":
                                dataOneField.append('nodal')
                            dataOneField.append(1)
                            dataOneField.append(fields['name'])
                            #display
                            logging.info("Prepare to write field: %s"%fields['name'])
                            #fill list to export
                            dataExport.append(dataOneField)
                        #export data
                        if fileName is not None:
                            fileName = fileName+".msh"
                            logging.info("Write: %s"%fileName)
                            silex_lib_gmsh.WriteResults2(fileName,
                                         self.fluidNodes,
                                         self.fluidElems,
                                         4,
                                         dataExport)
                            logging.info('File size: %s'%utils.file_size(fileName))
                        else:
                            logging.info("Unable to write file: filename is missing")
                    #
                    if method is "pos":
                        #write a pos files (gmsh syntax) which include the discontinuities
                        # one file per field
                        if dictFields is dict:
                            dictFields=[dictFields]
                        #prepare data to export
                        for fields in dictFields:
                            fileName = fields['filename']+'.pos'
                            logging.info("Prepare to write field: %s in %s"%(fields['name'],fileName))
                            silex_lib_xfem_acou_tet4.makexfemposfilefreq(
                                self.fluidNodes,
                                self.fluidElems, 
                                fields['levelsetmod'],
                                fields['fielduncorrected'],
                                fields['fieldenrichment'],
                                fileName,
                                fields['name'])
                            logging.info('File size: %s'%utils.file_size(fileName))

            fileOk = None
            if fileName is not None:
                fileOk = fileName

            if typeExport is "cavitymesh":                
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_air_cavity_Mesh',ext=None)
                #export mesh of the cavity
                if method is "msh":
                    silex_lib_gmsh.WriteResults(fileOk, self.fluidNodes, self.fluidElems, 4)
                    
            if typeExport is "controlvolmesh":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_air_controlled_volume_Mesh',ext=None)
                #export mesh of the control volume
                if method is "msh":                    
                    silex_lib_gmsh.WriteResults(fileOk, self.fluidNodes, self.fluidElemsControl, 4)
                    
            if typeExport is "struct":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_struc_surface',ext=None)
                #export 2D mesh of the structur
                if method is "msh":                    
                    silex_lib_gmsh.WriteResults2(fileOk, self.structNodes, self.structElems, 2)
            if typeExport is "levelset":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_LS_data',ext=None)
                #export level-set and gradients of the level-set
                if method is "msh":
                    #create list of data
                    dataW = list()
                    #level-set
                    dataW.append([[self.LevelSet], 'nodal', 1, 'Level set'])
                    #gradients of level-set
                    itP = 0
                    for iN in self.paraData['nameGrad']:
                        dataW.append([[self.LevelSetGradient[itP]],'nodal', 1, 'Level Set Grad '+iN])
                        itP = itP+1
                    silex_lib_gmsh.WriteResults2(fileOk, self.fluidNodes, self.fluidElems, 4, dataW)
            if typeExport is "enrichedPart":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_LS_enriched_elements',ext=None)
                if method is "msh":
                    silex_lib_gmsh.WriteResults2(fileOk,self.fluidNodes, self.fluidElems[self.EnrichedElems], 4)
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

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
            logging.error('Bad statement of boundary condition')
        #
        return nodesList,valuesOut





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

    def solveLinear(self, A, B):
        """
        Method used to solve linear problem (choose automatically the available approach)
        """
        if self.mumpsOk:
            sol = mumps.spsolve(A, B, comm=self.commMPI)
            # mumps is globally defined
        else:
            sol = scipy.sparse.linalg.spsolve(A, B)
        return sol

    

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
        self.FRF = [None for _ in range(self.caseProp['nbSteps'])]
        #
        if self.paraData['gradCompute']:
            self.pressureUncorrectGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureEnrichmentGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            self.pressureGrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]
            #
            self.FRFgrad=[np.zeros([self.fluidNbDofs,self.caseProp['nbSteps']],dtype=self.loadType()) for _ in range(self.getNbGrad())]

    def formatPara(self,valIn=None,nospace=False):
        """
        ##################################################################
        # Format parameters as string
        ##################################################################
        """
        #if necessary load the current values of parameters
        if valIn is None:
            valU = self.paraData['val']
        else:
            valU=valIn
        dispPara = ''
        #condition valU
        if len(valU.shape) == 0:
            valU=[valU]
        #
        if self.paraData:
            for itV in range(0,len(valU)):
                dispPara += self.paraData['name'][itV]+' '+valU[itV].astype(str)+' '
        else:
            for itV in range(0,len(valU)):
                dispPara += valU[itV].astype(str)+' '
        #
        if nospace:
            dispPara = dispPara.replace(' ','_')
        return dispPara
    
    def prepPara(self,valIn):
        """
        ##################################################################
        # prepare parameters values to be ran
        ##################################################################
        """
        #if necessary load the current values of parameters
        if len(valIn.shape)==1:
            if self.paraData['nb']==1:
                paraArray = valIn
            else:
                paraArray = np.array([valIn])
        else:
            paraArray = valIn
        return paraArray

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
            #
            logging.info("Time to solve the whole problem for set of parameters nb %i - %g s"%(self.nbRuns,time.process_time()-ticV))
        #
        logging.info("Time to solve the whole problem along sets of parameters - %g s"%(time.process_time()-ticS))
        
    def exportFieldsOnePara(self,typeSave=None,paraName=False):
        """
        ##################################################################
        # Export fields in mesh file for one set of parameters
        ##################################################################
        """
        if typeSave is None:
            typeSave = self.data['exportMesh']
            logging.info('Use default format:  %s'%(typeSave))
        #
        # information concerning parameters
        txtPara=self.formatPara()
        #create database to export
        dataW = list()
        dataW.append({'field':np.real(self.pressure),'type':'nodal','name':'Pressure (real) ('+txtPara+')'})
        dataW.append({'field':np.imag(self.pressure),'type':'nodal','name':'Pressure (imaginary) ('+txtPara+')'})
        dataW.append({'field':np.absolute(self.pressure),'type':'nodal','name':'Pressure (norm) ('+txtPara+')'})
        #
        if self.paraData['gradCompute']:
            for itG,txtP in np.ndenumerate(self.paraData['nameGrad']):
                dataW.append({'field':np.real(self.pressureGrad),'type':'nodal','name':'Grad. '+txtP+' Pressure (real) ('+txtPara+')'})
                dataW.append({'field':np.imag(self.pressureGrad),'type':'nodal','name':'Grad. '+txtP+' Pressure (imaginary) ('+txtPara+')'})
                dataW.append({'field':np.absolute(self.pressureGrad),'type':'nodal','name':'Grad. '+txtP+' Pressure (norm) ('+txtPara+')'})
        #write the file
        self.exportResults(typeExport="manuFields",method=typeSave,dictFields = dataW,fileName = self.resultFile(detPara=paraName,addTxt='results_fluid',ext=None))
        #write the discontinuties of the field in file
        #new database to export
        dataW = list()
        #
        dataW.append({'fielduncorrected':np.real(self.pressureUncorrect),
        'fieldenrichment':np.real(self.pressureEnrichment),
        'levelsetmod':None,
        'filename':self.resultFile(detPara=paraName,addTxt='plus_real',ext=None),
        'name':'Pressure + (real) ('+txtPara+')'})
        dataW.append({'fielduncorrected':np.imag(self.pressureUncorrect),
        'fieldenrichment':np.imag(self.pressureEnrichment),
        'levelsetmod':None,
        'filename':self.resultFile(detPara=paraName,addTxt='plus_imag',ext=None),
        'name':'Pressure + (imaginary) ('+txtPara+')'})
        dataW.append({'fielduncorrected':np.absolute(self.pressureUncorrect),
        'fieldenrichment':np.absolute(self.pressureEnrichment),
        'levelsetmod':None,
        'filename':self.resultFile(detPara=paraName,addTxt='plus_abs',ext=None),
        'name':'Pressure + (norm) ('+txtPara+')'})
        #
        dataW.append({'fielduncorrected':np.real(self.pressureUncorrect),
        'fieldenrichment':-np.real(self.pressureEnrichment),
        'levelsetmod':-self.LevelSet,
        'filename':self.resultFile(detPara=paraName,addTxt='moins_real',ext=None),
        'name':'Pressure - (real) ('+txtPara+')'})
        dataW.append({'fielduncorrected':np.imag(self.pressureUncorrect),
        'fieldenrichment':-np.imag(self.pressureEnrichment),
        'levelsetmod':-self.LevelSet,
        'filename':self.resultFile(detPara=paraName,addTxt='moins_imag',ext=None),
        'name':'Pressure - (imaginary) ('+txtPara+')'})
        dataW.append({'fielduncorrected':np.absolute(self.pressureUncorrect),
        'fieldenrichment':-np.absolute(self.pressureEnrichment),
        'levelsetmod':-self.LevelSet,
        'filename':self.resultFile(detPara=paraName,addTxt='moins_abs',ext=None),
        'name':'Pressure - (norm) ('+txtPara+')'})
        #
        if self.paraData['gradCompute']:
            for itG,txtP in np.ndenumerate(self.paraData['nameGrad']):
                dataW.append({'fielduncorrected':np.real(self.pressureUncorrectGrad[itG]),
                'fieldenrichment':np.real(self.pressureEnrichmentGrad[itG]),
                'levelsetmod':self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_real',ext=None),
                'name':'Pressure + (real) ('+txtPara+')'})
                dataW.append({'fielduncorrected':np.imag(self.pressureUncorrectGrad[itG]),
                'fieldenrichment':np.imag(self.pressureEnrichmentGrad[itG]),
                'levelsetmod':self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_imag',ext=None),
                'name':'Pressure + (imaginary) ('+txtPara+')'})
                #
                dataW.append({'fielduncorrected':np.real(self.pressureUncorrectGrad[itG]),
                'fieldenrichment':-np.real(self.pressureEnrichmentGrad[itG]),
                'levelsetmod':-self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_real',ext=None),
                'name':'Pressure - (real) ('+txtPara+')'})
                dataW.append({'fielduncorrected':np.imag(self.pressureUncorrectGrad[itG]),
                'fieldenrichment':-np.imag(self.pressureEnrichmentGrad[itG]),
                'levelsetmod':-self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_imag',ext=None),
                'name':'Pressure - (imaginary) ('+txtPara+')'})
                #compute gradients of absolute values
                gradAbsUncorrect = (np.real(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect)+\
                    np.imag(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect))/np.absolute(self.pressureUncorrect)
                # deal with zeros values in enrichment field
                gradAbsEnrichment = np.zeros([self.fluidNodes, self.caseProp['nbSteps']])
                gradAbsEnrichment[self.SolvedDofA, :] = (np.real(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.real(self.pressureEnrichment[self.SolvedDofA, :])+\
                    np.imag(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.imag(self.pressureEnrichment[self.SolvedDofA, :]))/np.absolute(self.pressureEnrichment[self.SolvedDofA, :])
                # remove inf value
                IX = scipy.absolute(self.pressureUncorrect) == 0.
                gradAbsUncorrect[IX] = 1
                gradAbsEnrichment[IX] = 1
                #
                dataW.append({'fielduncorrected':self.gradAbsUncorrect[itG],
                'fieldenrichment':self.gradAbsEnrichment[itG],
                'levelsetmod':self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_abs',ext=None),
                'name':'Pressure + (norm) ('+txtPara+')'})
                dataW.append({'fielduncorrected':self.gradAbsUncorrect[itG],
                'fieldenrichment':-self.gradAbsEnrichment[itG],
                'levelsetmod':-self.LevelSet,
                'filename':self.resultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_abs',ext=None),
                'name':'Pressure - (norm) ('+txtPara+')'})

        #write the file
        self.exportResults(typeExport='manuFields',method='pos',dictFields = dataW,fileName = self.resultFile(detPara=paraName,addTxt='results_fluid',ext=None))


    def saveFRF(self,typeSave=None,paraName=False,allData=False):
        """
        ##################################################################
        # Function used to solve FRF in mat file or pickle
        ##################################################################
        """
        if typeSave is None:
            typesave = self.data['exportData']
            logging.info('Use default format:  %s'%(typeSave))
        if self.caseProp['computeFRF']:            
            #build dictionary
            dictOut=self.paraData
            dictOut['frequencies']=self.Frequencies
            dictOut['FRF']=self.FRF
            dictOut['FRFgrad']=self.FRFgrad
            if typeSave is 'mat':
                #export data
                scipy.io.savemat(self.resultFile(detPara=paraName,addTxt='results',ext='mat'),mdict=dictOut)
                logging.info('Export data in %s'%(self.resultFile(detPara=paraName,addTxt='results',ext='mat')))
            if typeSave is 'pickle':
                #export data
                f = open(self.resultFile(detPara=paraName,addTxt='results',ext='pck'))
                pickle.dump(dictOut,f)
                f.close()
                logging.info('Export data in %s'%(self.resultFile(detPara=paraName,addTxt='results',ext='pck')))
        else:
            logging.info('Nothing to export (FRF not computed)')

    def plotFRF(self):
        """
        ##################################################################
        # Function used to plot FRF
        ##################################################################
        """
        pass
    
    def saveResults(self,):
        """
        Method used to save results
        """
        pass

    



# ###########################################################
# ###########################################################
# ###########################################################
# ###########################################################
# ###########################################################
# ###########################################################




#     ##############################################################
#     # Material, Boundary conditions
#     ##############################################################

    

#     ##############################################################
#     # Load fluid mesh
#     ##############################################################
        

#     # ##############################################################
#     # # Load structure mesh
#     # ##############################################################

#     ##################################################################
#     # compute level set
#     ##################################################################

#     ##################################################################
#     # Get enriched nodes and elements
#     ##################################################################
   

    

#     ##################################################################
#     # Build Second member
#     ##################################################################

#     # To impose the load on the fluid:
#     # fluid node number 1
    

#     #################################################################
#     # Compute gradients with respect to parameters
#     ##################################################################
#     # print(silex_lib_xfem_acou_tet4.globalacousticgradientmatrices.__doc__)

        

#     ##############################################################
#     # FRF computation
#     ##############################################################

#     Flag_frf_analysis = 1
#     frequencies = []
#     frf = []
#     frfgradient = list()
#     for it in range(0, nbPara):
#         frfgradient.append([])

#     if (Flag_frf_analysis == 1):
#         print("Proc. ", rank, " / time at the beginning of the FRF:", time.ctime())

#         if rank == 0:
#             print('nb of total dofs: ', len(SolvedDofF)+len(SolvedDofA))

#         press_save = []
#         enrichment_save = []
#         uncorrectedpress_save = []
#         disp_save = []
#         denrichment_save = []
#         duncorrectedpress_save = []
#         dpress_save = list()
#         for it in range(0, nbPara):
#             dpress_save.append([])
#             denrichment_save.append([])
#             duncorrectedpress_save.append([])

#         # extract frequencies for the associated processors
#         freqCompute = listFreqPerProc[:, rank]
#         freqCompute = freqCompute[freqCompute > 0]
#         it = 0
#         itmax = len(freqCompute)
#         for freq in freqCompute:
#             it = it+1
#             #freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
#             frequencies.append(freq)
#             omega = 2*scipy.pi*freq

#             print("Freq. step ", it, " proc number", rank, "frequency=", freq)

#             tic = time.clock()

#             F = scipy.array(omega**2*UF[SolvedDof], dtype='c16')

#             if rank >= 0:
#                 # print(K)
#                 # print(M)
#                 # print(omega)
#                 sol = mumps.spsolve(K-(omega**2)*M, F+0.j, comm=mycomm)
#                 # sol = mumps.spsolve(scipy.sparse.coo_matrix( \
#                 #   K-(omega**2)*M, dtype='complex'), F+0.j,comm=mycomm)
#                 # sol
#                 # sol = scipy.sparse.linalg.spsolve(scipy.sparse.csc_matrix(
#                 #     K-(omega**2)*M, dtype='c16'), F)

#             # pressure field without enrichment
#             press1 = np.zeros((fluid_ndof), dtype=complex)
#             press1[SolvedDofF] = sol[list(range(len(SolvedDofF)))].copy()
#             # enrichment field
#             enrichment = np.zeros((fluid_nnodes), dtype=complex)
#             enrichment[SolvedDofA] = sol[list(
#                 range(len(SolvedDofF), len(SolvedDofF)+len(SolvedDofA)))].copy()
#             # correction of the pressure field with enrichment
#             # np.zeros((fluid_ndof),dtype=complex) #press1.copy()
#             CorrectedPressure = press1.copy()
#             CorrectedPressure[SolvedDofA] = press1[SolvedDofA] + \
#                 enrichment[SolvedDofA]*np.sign(LevelSet[SolvedDofA])

#             # compute and store FRF on the test volume
#             # frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements5,fluid_nodes,CorrectedPressure))
#             frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(
#                 fluid_elements5, fluid_nodes, press1, enrichment, LevelSet, LevelSet*0-1.0))

#             if (flag_write_gmsh_results == 1) and (rank == 0):
#                 press_save.append(CorrectedPressure.copy())
#                 enrichment_save.append(enrichment.copy())
#                 uncorrectedpress_save.append(press1.copy())

#             #####################
#             #####################
#             ######################
#             #####################
#             Dpress_Dtheta = np.zeros([fluid_ndof, nbPara], dtype=complex)
#             DCorrectedPressure_Dtheta = scipy.array(Dpress_Dtheta)
#             Denrichment_Dtheta = np.zeros(
#                 [fluid_ndof, nbPara], dtype=complex)
#             #####################
#             #####################
#             # compute gradients
#             for itP in range(0, nbPara):
#                 # solve gradient problem
#                 tmp = -(dK[itP]-(omega**2)*dM[itP])*sol
#                 Dsol_Dtheta_RAW = mumps.spsolve(
#                     K-(omega**2)*M, tmp, comm=mycomm)
#                 #Dsol_Dtheta_RAW = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
#                 # Dsol_Dtheta_RAW = scipy.sparse.linalg.spsolve( scipy.sparse.csc_matrix(K-(omega**2)*M,dtype='c16')  , tmp )
#                 #####################
#                 #####################
#                 # gradient of the pressure field without enrichment
#                 Dpress_Dtheta[SolvedDofF, itP] = Dsol_Dtheta_RAW[list(
#                     range(len(SolvedDofF)))].copy()
#                 #####################
#                 #####################
#                 # gradient of the enrichment field
#                 Denrichment_Dtheta[SolvedDofA, itP] = Dsol_Dtheta_RAW[list(
#                     range(len(SolvedDofF), len(SolvedDofF)+len(SolvedDofA)))].copy()
#                 #####################
#                 #####################
#                 # compute the corrected gradient pressure field (via enrichment)
#                 DCorrectedPressure_Dtheta[:, itP] = scipy.array(
#                     Dpress_Dtheta[:, itP].copy())
#                 DCorrectedPressure_Dtheta[SolvedDofA, itP] = DCorrectedPressure_Dtheta[SolvedDofA, itP].T + \
#                     scipy.array(
#                         Denrichment_Dtheta[SolvedDofA, itP]*np.sign(LevelSet[SolvedDofA]).T)
#                 #####################
#                 #####################
#                 # store gradients
#                 frfgradient[itP].append(
#                     silex_lib_xfem_acou_tet4.computegradientcomplexquadratiquepressure(
#                         fluid_elements5,
#                         fluid_nodes,
#                         press1+0j,
#                         Dpress_Dtheta[:, itP]+0j,
#                         LevelSet))
#                 #####################
#                 #####################
#                 dpress_save[itP].append(
#                     DCorrectedPressure_Dtheta[:, itP].copy())
#                 denrichment_save[itP].append(Denrichment_Dtheta[:, itP].copy())
#                 duncorrectedpress_save[itP].append(
#                     Dpress_Dtheta[:, itP].copy())

#         frfsave = [frequencies, frf, frfgradient]
#         if rank != 0:
#             comm.send(frfsave, dest=0, tag=11)

#         print("Proc. ", rank, " / time at the end of the FRF:", time.ctime())

#         if (flag_write_gmsh_results == 1) and (rank == 0):
#             dataW = list()
#             # prepare pressure field
#             dataW.append(
#                 [scipy.real(press_save), 'nodal', 1, 'pressure (real)'])
#             dataW.append([scipy.imag(press_save), 'nodal',
#                           1, 'pressure (imaginary)'])
#             dataW.append([scipy.absolute(press_save),
#                           'nodal', 1, 'pressure (norm)'])
#             # prepare gradient pressure field
#             itG = 0
#             for itP in NamePara:
#                 dataW.append([scipy.real(dpress_save[itG]), 'nodal',
#                               1, 'pressure gradient '+itP+' (real)'])
#                 dataW.append([scipy.imag(dpress_save[itG]), 'nodal',
#                               1, 'pressure gradient '+itP+' (imaginary)'])
#                 dataW.append([scipy.absolute(dpress_save[itG]),
#                               'nodal', 1, 'pressure gradient '+itP+' (norm)'])
#                 itG = itG+1
#             print("Write pressure field and gradients in msh file")
#             silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',
#                                          fluid_nodes, fluid_elements1, 4, dataW)
#             print(">>> Done!!")

#             # export results with discontinuities on .pos files
#             varExport = scipy.vstack(uncorrectedpress_save).transpose()
#             varExportC = scipy.vstack(press_save).transpose()
#             varExportB = scipy.vstack(enrichment_save).transpose()
#             print("Write pressure field in pos file")
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, LevelSet, scipy.real(
#                 varExport), scipy.real(varExportB), 'results/press_plus_real.pos', 'Pressure + Real')
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, LevelSet, scipy.imag(
#                 varExport), scipy.imag(varExportB), 'results/press_plus_imag.pos', 'Pressure + Imag')
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, LevelSet, scipy.absolute(
#                 varExport), scipy.absolute(varExportB), 'results/press_plus_abs.pos', 'Pressure + Abs')
#             # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,scipy.absolute(varExport)**2,scipy.absolute(varExportB)**2,'press_plus_square.pos')
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, -LevelSet, scipy.real(
#                 varExport), -scipy.real(varExportB), 'results/press_moins_real.pos', 'Pressure - Real')
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, -LevelSet, scipy.imag(
#                 varExport), -scipy.imag(varExportB), 'results/press_moins_imag.pos', 'Pressure - Imag')
#             silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, -LevelSet, scipy.absolute(
#                 varExport), -scipy.absolute(varExportB), 'results/press_moins_abs.pos', 'Pressure - Abs')
#             print(">>> Done!!")
#             #
#             itG = 0
#             for key, itP in enumerate(NamePara):
#                 GvarExport = scipy.vstack(
#                     duncorrectedpress_save[key]).copy().transpose()
#                 GvarExportB = scipy.vstack(
#                     denrichment_save[key]).copy().transpose()
#                 print("Write gradient of pressure field in pos file (", itP, ")")
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, LevelSet, scipy.real(
#                     GvarExport), scipy.real(GvarExportB), 'results/Gpress_plus_'+itP+'_real.pos', 'Gpressure + '+itP+' Real')
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, LevelSet, scipy.imag(
#                     GvarExport), scipy.imag(GvarExportB), 'results/Gpress_plus_'+itP+'_imag.pos', 'Gpressure + '+itP+' Imag')
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, -LevelSet, scipy.real(
#                     GvarExport), -scipy.real(GvarExportB), 'results/Gpress_moins_'+itP+'_real.pos', 'Gpressure - '+itP+' Real')
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes, fluid_elements1, -LevelSet, scipy.imag(
#                     GvarExport), -scipy.imag(GvarExportB), 'results/Gpress_moins_'+itP+'_imag.pos', 'Gpressure - '+itP+' Imag')
#                 #
#                 # gradPsquare=2*(scipy.real(GvarExport)*scipy.real(varExport)+scipy.imag(GvarExport)*scipy.imag(varExport))
#                 # gradPsquareB=2*(scipy.real(GvarExportB)*scipy.real(varExportB)+scipy.imag(GvarExportB)*scipy.imag(varExportB))
#                 # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,LevelSet,gradPsquare,gradPsquareB,'results/Gpress_plus_'+itP+'_dsquare.pos')
#                 # silex_lib_xfem_acou_tet4.makexfemposfilefreq(fluid_nodes,fluid_elements1,-LevelSet,gradPsquare,-gradPsquareB,'results/Gpress_moins_'+itP+'_dsquare.pos')
#                 #
#                 gradCalc = (scipy.real(GvarExport)*scipy.real(varExport)+scipy.imag(
#                     GvarExport)*scipy.real(varExport))/scipy.absolute(varExport)
#                 # deal with zeros values in enrichment field
#                 gradCalcB = np.zeros([fluid_nnodes, nbStep])
#                 gradCalcB[SolvedDofA, :] = (scipy.real(GvarExportB[SolvedDofA, :])*scipy.real(varExportB[SolvedDofA, :])+scipy.imag(
#                     GvarExportB[SolvedDofA, :])*scipy.imag(varExportB[SolvedDofA, :]))/scipy.absolute(varExportB[SolvedDofA, :])
#                 # remove inf value
#                 IX = scipy.absolute(varExport) == 0.
#                 gradCalc[IX] = 1
#                 gradCalcB[IX] = 1
#                 #
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(
#                     fluid_nodes, fluid_elements1, LevelSet, gradCalc, gradCalcB, 'results/Gpress_plus_'+itP+'_dabsolute.pos', 'Gpressure + '+itP+' dAbs')
#                 silex_lib_xfem_acou_tet4.makexfemposfilefreq(
#                     fluid_nodes, fluid_elements1, -LevelSet, gradCalc, -gradCalcB, 'results/Gpress_moins_'+itP+'_dabsolute.pos', 'Gpressure - '+itP+' dAbs')
#                 print(">>> Done!!")
#                 itG = itG+1

        

# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# # function for dealing with options
# def manageOpt(argv, dV):
#     # load default values
#     freqMin = dV.freqMin
#     freqMax = dV.freqMax
#     nbStep = dV.nbStep
#     paraVal = scipy.array(dV.paraVal)
#     gradCompute = scipy.array(dV.gradCompute)
#     #caseDefine = dV.caseDef

#     # load info from MPI
#     nbProc, rank, comm = mpiInfo()
#     # load options
#     opts, args = getopt.getopt(argv, "p:s:F:f:hp:c:g:")
#     for opt, arg in opts:
#         if opt == "-s":
#             nbStep = int(arg)
#         elif opt == "-F":
#             freqMax = float(arg)
#         elif opt == "-f":
#             freqMin = float(arg)
#         elif opt == "-p":
#             tmp = scipy.array(arg.split(','), dtype=scipy.float32)
#             paraVal = tmp
#         elif opt == "-c":
#             caseDefine = str(arg)
#         elif opt == "-g":
#             tmp = scipy.array(arg.split(','), dtype=scipy.int32)
#             gradCompute = tmp
#         elif opt == "-h":
#             usage()
#             sys.exit()
#     # print chosen parameters
#     print("Number of processors: ", nbProc)
#     print("Parameters: ", paraVal)
#     print("Number of frequency steps: ", nbStep)
#     print("Maximum frequency: ", freqMax)
#     print("Minimum frequency: ", freqMin)
#     print("Components of grad: ", gradCompute)
#     #print ("Case: ",caseDefine)
#     it = 0
#     for itP in paraVal:
#         print('Parameter num '+str(it)+': '+str(itP))
#         it = it+1
#     print("\n\n")

#     # run computation
#     RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm,
#           paraVal, gradCompute, 1)  # ,caseDefine)

# # usage definition


# def usage():
#     dV = defaultV
#     print("Usage: ", sys.argv[0], "-psFfhg [+arg]")
#     print("\t -p : input parameters (default value ", dV.nbProc, ")")
#     print("\t -s : number of steps in the frequency range (default value ", dV.nbStep, ")")
#     print("\t -F : maximum frequency (default value ", dV.freqMax, ")")
#     print("\t -f : minimum frequency (default value ", dV.freqMin, ")")
#     print("\t -g : Components of grad (default value ", dV.gradCompute, ")")

# # default values


# class defaultV:
#     freqMin = 10.0
#     freqMax = 150.0
#     nbStep = 1000
#     paraVal = [2., 2., 1., 1.]
#     gradCompute = [0, 1, 2, 3]
#     nbProc = 1
#     #caseDef= 'thick_u'


# # Run autonomous
# if __name__ == '__main__':
#     # run with options
#     dV = defaultV
#     manageOpt(sys.argv[1:], dV)
