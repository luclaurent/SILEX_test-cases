###########################################################
# Declaration of useful tools for SILEXclass
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import numpy as np
import os
import logging
from importlib import util as utilImport

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
    mpi4py_loader=utilImport.find_spec('mpi4py')
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
    mumps_loader=utilImport.find_spec('mumps')
    foundMumps = mumps_loader is not None
    clMumps = None
    if foundMumps:
        # globals()['mumps'] = __import__('mumps')
        import mumps as clMumps
    return foundMumps,clMumps

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

class tools:


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
    
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def loadMPI(self):
        """
        ##################################################################
        # method used to load openMPI information
        ##################################################################
        """
        self.nbProcMPI,self.rankMPI,self.commMPI=mpiInfo()

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def loadMUMPS(self):
        """
        ##################################################################
        # method used to load openMPI information
        ##################################################################
        """
        self.mumpsOk,self.classMumps=loadMumps()

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
    
    def getResultFile(self,detPara=False,addTxt=None,ext=None):
        """
        ##################################################################
        # Build result file depending on status
        ##################################################################
        """
        fileName = self.data['prefixResults']
        #add parameters value in filename
        if detPara:
            fileName += '_'+self.formatPara(nospace=True)
        #add given text in filename
        if addTxt is not None:
            if addTxt[0] != '_' and fileName[-1] != '_':
                fileName += '_' + addTxt
            else:
                fileName += addTxt
        #add extension in filename
        if ext is not None:
            fileName += '.'+ext
        #store the full path of the file
        if self.fullPathCurrentResultsFolder == '':
            folder = self.data['resultsFolder']
        else:
            folder = self.fullPathCurrentResultsFolder
        self.fullPathResultsFile = os.path.join(folder,fileName)
        #
        return self.fullPathResultsFile

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


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
        #remove last space
        dispPara=dispPara[:-1]
        #
        if nospace:
            dispPara = dispPara.replace(' ','_')
        return dispPara

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
