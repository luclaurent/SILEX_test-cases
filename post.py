###########################################################
# Declaration of useful post-processing tools for SILEXclass
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################
import logging
import time
import os
import re
import scipy
import scipy.io
import numpy as np
import pickle
#
import utils
#
from SILEX import classMeshField_mod as cl
from SILEX import silex_lib_xfem_acou_tet4

class postProcess:

    def writeMeshField(self,writerMethod=None,uncorrectedField=None,enrichmentField=None,fileName=None,fieldName=None):
        """
        ##################################################################
        # function used to export results and mesh via fortran classMeshField
        ##################################################################
        """
        #class used to catch the stdout of Fortran
        clSTD = utils.RedirectFortran('logfortran.tmp.log')
        funExport= lambda x: logging.info('%s%s'%('FORTRAN: ',x))
        
        #first initialization
        if not self.classSave:
            if self.LevelSet is not None and self.fluidNodes is not None and self.fluidElems is not None:
                
                #initialized class for results export                
                clSTD.start()
                self.classSave = cl.classmeshfield
                self.classSave.init(self.fluidNodes,self.fluidElems,self.LevelSet)
                clSTD.stop(funExport)
            else:
                logging.error('Unable to initialize the saving class due to a lack of data')
        #declare the type of output
        if writerMethod is not None and fileName is not None:
            if writerMethod is "msh" or writerMethod is "mshv2":
                typeExport="msh"
            if writerMethod is "vtk":                
                typeExport="vtk"
            #declare the writer
            clSTD.start()
            self.classSave.declarewriter(fileName,typeExport)
            clSTD.stop(funExport)
        else:
            if writerMethod is None:
                logging.warning('Writer method not declared (use %s)'%self.classSave.typeWriter)
            if fileName is None:
                logging.warning('fileName method not declared  (use %s)'%self.classSave.outputFileName)
        #append fields to the existing file(s)
        writeOk=True
        if uncorrectedField is not None and enrichmentField is not None:
            writeOk=True
            clSTD.start()
            if fieldName is not None:
                self.classSave.loadfields(uncorrectedField,enrichmentField,fieldName.replace(' ','_'))
            else:
                self.classSave.loadfields(uncorrectedField,enrichmentField)
            clSTD.stop(funExport)
        #append field without correction
        if uncorrectedField is not None and enrichmentField is None:
            writeOk=True
            clSTD.start()
            if fieldName is not None:
                self.classSave.loadfieldscorrected(uncorrectedField,fieldName.replace(' ','_'))
            else:
                self.classSave.loadfieldscorrected(uncorrectedField)
            clSTD.stop(funExport)
        #create list of written files
        listFiles=None
        if writeOk:            
            file = str(self.classSave.outputfilename.astype('str')).replace(' ','')
            basename = os.path.relpath(file,self.fullPathCurrentResultsFolder).replace('.msh','').replace('.vtk','')
            filelist= [f for f in os.listdir(self.fullPathCurrentResultsFolder) if re.match(basename+'.*\.(msh|vtk)', f)]
            return filelist

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


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
                            fileName = utils.addExt(fileName,".msh")
                            logging.info("Write: %s"%fileName)
                            #
                            if os.path.exists(fileName):
                                logging.warning(">>> File %s will OVERWRITTEN"%fileName)
                            #
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
                            fileName = utils.addExt(fields['filename'],'.pos')
                            logging.info("Prepare to write field: %s in %s"%(fields['name'],fileName))
                            #check if data are well shaped
                            dataLS = fields['levelsetmod']
                            dataUnc = fields['fielduncorrected']
                            dataEnr = fields['fieldenrichment']
                            #
                            if dataUnc.shape[1] == self.fluidNbNodes:
                                funU=lambda x: x.transpose()
                                logging.warning("Change shape of uncorrected fields")
                            elif dataUnc.shape[0] == self.fluidNbNodes:
                                funU=lambda x: x
                            else:
                                logging.error("Bad dimension of uncorrected fields to be exported")
                            #
                            if dataEnr.shape[1] == self.fluidNbNodes:
                                funE=lambda x: x.transpose()
                                logging.warning("Change shape of enrichment fields")
                            elif dataEnr.shape[0] == self.fluidNbNodes:
                                funE=lambda x: x
                            else:
                                logging.error("Bad dimension of enrichment fields to be exported")
                            #                            
                            if os.path.exists(fileName):
                                logging.warning(">>> File %s will OVERWRITTEN"%fileName)
                            #
                            silex_lib_xfem_acou_tet4.makexfemposfilefreq(
                                self.fluidNodes,
                                self.fluidElems, 
                                dataLS,
                                funU(dataUnc),
                                funE(dataEnr),
                                fileName,
                                fields['name'])
                            logging.info('File size: %s'%utils.file_size(fileName))
                    #
                    if method is "mshv2" or method is "vtk":
                        #write a pos files (gmsh syntax) which include the discontinuities
                        # one file per field
                        if dictFields is dict:
                            dictFields=[dictFields]
                        #prepare data to export
                        for fields in dictFields:
                            #adapt extension
                            if method is "mshv2":
                                fileName = utils.addExt(fileName,'.msh')
                            if method is "vtk":
                                fileName = utils.addExt(fileName,'.vtk')
                            logging.info("Prepare to write field: %s in %s"%(fields['name'],fileName))
                            #check if data are well shaped
                            dataUnc = fields['fielduncorrected']
                            dataEnr = fields['fieldenrichment']
                            #
                            if dataUnc.shape[1] == self.fluidNbNodes:
                                funU=lambda x: x.transpose()
                                logging.warning("Change shape of uncorrected fields")
                            elif dataUnc.shape[0] == self.fluidNbNodes:
                                funU=lambda x: x
                            else:
                                logging.error("Bad dimension of uncorrected fields to be exported")
                            #
                            if dataEnr.shape[1] == self.fluidNbNodes:
                                funE=lambda x: x.transpose()
                                logging.warning("Change shape of enrichment fields")
                            elif dataEnr.shape[0] == self.fluidNbNodes:
                                funE=lambda x: x
                            else:
                                logging.error("Bad dimension of enrichment fields to be exported")
                            #                            
                            if os.path.exists(fileName):
                                logging.warning(">>> File %s will OVERWRITTEN"%fileName)
                            #
                            fileList = self.writeMeshField(
                                method,
                                funU(dataUnc),
                                funE(dataEnr),
                                fileName,
                                fields['name'])
                            for file in fileList:
                                logging.info('File size: %s (%s)'%(utils.file_size(os.path.join(self.fullPathCurrentResultsFolder,file)),file))

            fileOk = None
            if fileName is not None:
                fileOk = fileName

            if typeExport is "cavitymesh":                
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_air_cavity_Mesh',ext='msh')
                #export mesh of the cavity
                if method is "msh":
                    if os.path.exists(fileOk):
                        logging.warning(">>> File %s will OVERWRITTEN"%fileOk)
                    #
                    silex_lib_gmsh.WriteResults(fileOk, self.fluidNodes, self.fluidElems, 4)
                    
            if typeExport is "controlvolmesh":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_air_controlled_volume_Mesh',ext='msh')
                #export mesh of the control volume
                if method is "msh":
                    if os.path.exists(fileOk):
                        logging.warning(">>> File %s will OVERWRITTEN"%fileOk)
                    #
                    silex_lib_gmsh.WriteResults(fileOk, self.fluidNodes, self.fluidElemsControl, 4)
                    #
            if typeExport is "struct":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_struc_surface',ext='msh')
                #export 2D mesh of the structur
                if method is "msh":                    
                    if os.path.exists(fileOk):
                        logging.warning(">>> File %s will OVERWRITTEN"%fileOk)
                    #
                    silex_lib_gmsh.WriteResults2(fileOk, self.structNodes, self.structElems, 2)
                    #
            if typeExport is "levelset":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_LS_data',ext='msh')
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
                    #
                    if os.path.exists(fileOk):
                        logging.warning(">>> File %s will OVERWRITTEN"%fileOk)
                    #
                    silex_lib_gmsh.WriteResults2(fileOk, self.fluidNodes, self.fluidElems, 4, dataW)
                    #
            if typeExport is "enrichedPart":
                if not fileOk:
                    fileOk = self.getResultFile(detPara=True,addTxt='_LS_enriched_elements',ext='msh')
                if method is "msh":
                    #
                    if os.path.exists(fileOk):
                        logging.warning(">>> File %s will OVERWRITTEN"%fileOk)
                    #
                    silex_lib_gmsh.WriteResults2(fileOk,self.fluidNodes, self.fluidElems[self.EnrichedElems], 4)
            logging.info("++++++++++++++++++++ Done - %g s"%(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

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
        dataW.append({
            'field':np.real(self.pressure),
            'fielduncorrected':np.real(self.pressureUncorrect),
            'fieldenrichment':np.real(self.pressureEnrichment),
            'type':'nodal',
            'name':'Pressure (real) ('+txtPara+')'})
        dataW.append({
            'field':np.imag(self.pressure),
            'fielduncorrected':np.imag(self.pressureUncorrect),
            'fieldenrichment':np.imag(self.pressureEnrichment),
            'type':'nodal',
            'name':'Pressure (imaginary) ('+txtPara+')'})
        dataW.append({
            'field':np.absolute(self.pressure),
            'fielduncorrected':np.absolute(self.pressureUncorrect),
            'fieldenrichment':np.absolute(self.pressureEnrichment),
            'type':'nodal',
            'name':'Pressure (norm) ('+txtPara+')'})
        #
        if self.paraData['gradCompute']:
            for itG,txtP in np.ndenumerate(self.paraData['nameGrad']):
                dataW.append({
                    'field':np.real(self.pressureGrad),
                    'fielduncorrected':np.real(self.pressureUncorrectGrad),
                    'fieldenrichment':np.real(self.pressureEnrichmentGrad),
                    'type':'nodal',
                    'name':'Grad. '+txtP+' Pressure (real) ('+txtPara+')'})
                dataW.append({
                    'field':np.imag(self.pressureGrad),
                    'fielduncorrected':np.imag(self.pressureUncorrectGrad),
                    'fieldenrichment':np.imag(self.pressureEnrichmentGrad),
                    'type':'nodal',
                    'name':'Grad. '+txtP+' Pressure (imaginary) ('+txtPara+')'})
                dataW.append({
                    'field':np.absolute(self.pressureGrad),
                    'fielduncorrected':np.absolute(self.pressureUncorrectGrad),
                    'fieldenrichment':np.absolute(self.pressureEnrichmentGrad),
                    'type':'nodal',
                    'name':'Grad. '+txtP+' Pressure (norm) ('+txtPara+')'})
        #write the file
        self.exportResults(typeExport="manuFields",method=typeSave,dictFields = dataW,fileName = self.getResultFile(detPara=paraName,addTxt='results_fluid',ext=None))
        
        if self.debug:
            #write the discontinuties of the field in file
            #new database to export
            dataW = list()
            #
            dataW.append({'fielduncorrected':np.real(self.pressureUncorrect),
            'fieldenrichment':np.real(self.pressureEnrichment),
            'levelsetmod':self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='plus_real',ext=None),
            'name':'Pressure + (real) ('+txtPara+')'})
            dataW.append({'fielduncorrected':np.imag(self.pressureUncorrect),
            'fieldenrichment':np.imag(self.pressureEnrichment),
            'levelsetmod':self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='plus_imag',ext=None),
            'name':'Pressure + (imaginary) ('+txtPara+')'})
            dataW.append({'fielduncorrected':np.absolute(self.pressureUncorrect),
            'fieldenrichment':np.absolute(self.pressureEnrichment),
            'levelsetmod':self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='plus_abs',ext=None),
            'name':'Pressure + (norm) ('+txtPara+')'})
            #
            dataW.append({'fielduncorrected':np.real(self.pressureUncorrect),
            'fieldenrichment':-np.real(self.pressureEnrichment),
            'levelsetmod':-self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='moins_real',ext=None),
            'name':'Pressure - (real) ('+txtPara+')'})
            dataW.append({'fielduncorrected':np.imag(self.pressureUncorrect),
            'fieldenrichment':-np.imag(self.pressureEnrichment),
            'levelsetmod':-self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='moins_imag',ext=None),
            'name':'Pressure - (imaginary) ('+txtPara+')'})
            dataW.append({'fielduncorrected':np.absolute(self.pressureUncorrect),
            'fieldenrichment':-np.absolute(self.pressureEnrichment),
            'levelsetmod':-self.LevelSet,
            'filename':self.getResultFile(detPara=paraName,addTxt='moins_abs',ext=None),
            'name':'Pressure - (norm) ('+txtPara+')'})
            #
            if self.paraData['gradCompute']:
                for itG,txtP in np.ndenumerate(self.paraData['nameGrad']):
                    dataW.append({'fielduncorrected':np.real(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment':np.real(self.pressureEnrichmentGrad[itG]),
                    'levelsetmod':self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_real',ext=None),
                    'name':'Pressure + (real) ('+txtPara+')'})
                    dataW.append({'fielduncorrected':np.imag(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment':np.imag(self.pressureEnrichmentGrad[itG]),
                    'levelsetmod':self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_imag',ext=None),
                    'name':'Pressure + (imaginary) ('+txtPara+')'})
                    #
                    dataW.append({'fielduncorrected':np.real(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment':-np.real(self.pressureEnrichmentGrad[itG]),
                    'levelsetmod':-self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_real',ext=None),
                    'name':'Pressure - (real) ('+txtPara+')'})
                    dataW.append({'fielduncorrected':np.imag(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment':-np.imag(self.pressureEnrichmentGrad[itG]),
                    'levelsetmod':-self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_imag',ext=None),
                    'name':'Pressure - (imaginary) ('+txtPara+')'})
                    #compute gradients of absolute values
                    gradAbsUncorrect = (np.real(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect)+\
                        np.imag(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect))/np.absolute(self.pressureUncorrect)
                    # deal with zeros values in enrichment field
                    gradAbsEnrichment = np.zeros([self.fluidNodes, self.caseProp['nbSteps']])
                    gradAbsEnrichment[self.SolvedDofA, :] = (np.real(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.real(self.pressureEnrichment[self.SolvedDofA, :])+\
                        np.imag(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.imag(self.pressureEnrichment[self.SolvedDofA, :]))/np.absolute(self.pressureEnrichment[self.SolvedDofA, :])
                    # remove inf value
                    IX = np.absolute(self.pressureUncorrect) == 0.
                    gradAbsUncorrect[IX] = 1
                    gradAbsEnrichment[IX] = 1
                    #
                    dataW.append({'fielduncorrected':self.gradAbsUncorrect[itG],
                    'fieldenrichment':self.gradAbsEnrichment[itG],
                    'levelsetmod':self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'plus_abs',ext=None),
                    'name':'Pressure + (norm) ('+txtPara+')'})
                    dataW.append({'fielduncorrected':self.gradAbsUncorrect[itG],
                    'fieldenrichment':-self.gradAbsEnrichment[itG],
                    'levelsetmod':-self.LevelSet,
                    'filename':self.getResultFile(detPara=paraName,addTxt='_grad_'+str(itG)+'moins_abs',ext=None),
                    'name':'Pressure - (norm) ('+txtPara+')'})

            #write the file
            self.exportResults(typeExport='manuFields',method='pos',dictFields = dataW,fileName = self.getResultFile(detPara=paraName,addTxt='results_fluid',ext=None))


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def saveFRF(self,typeSave=None,paraName=False,allData=False):
        """
        ##################################################################
        # Function used to solve FRF in mat file or pickle
        ##################################################################
        """
        if typeSave is None:
            typeSave = self.data['exportData']
            logging.info('Use default format:  %s'%(typeSave))
        if self.caseProp['computeFRF']:            
            dictOut=self.getOutputDict(allData)
            #
            txt=''
            if allData:
                txt='_all'
            if typeSave is 'mat':
                #clean dictionary
                utils.replace_none(dictOut)
                #export data
                scipy.io.savemat(self.getResultFile(detPara=paraName,addTxt='results'+txt,ext='mat'),mdict=dictOut)
                logging.info('Export data in %s'%(self.getResultFile(detPara=paraName,addTxt='results'+txt,ext='mat')))
            if typeSave is 'pickle':
                #export data
                f = open(self.getResultFile(detPara=paraName,addTxt='results'+txt,ext='pck'))
                pickle.dump(dictOut,f)
                f.close()
                logging.info('Export data in %s'%(self.getResultFile(detPara=paraName,addTxt='results'+txt,ext='pck')))
        else:
            logging.info('Nothing to export (FRF not computed)')


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def getOutputDict(self,allData=False):
        """
        ##################################################################
        # Function used to create an output dictionary
        ##################################################################
        """
        #build dictionary
        dictOut=self.paraData
        if allData: 
            dictOut['frequencies']=self.Frequencies
            dictOut['FRF']=self.FRF
            dictOut['FRFgrad']=self.FRFgrad
            dictOut['paraVal']=self.paraData['val']
        else:
            dictOut['frequencies']=self.Frequencies
            dictOut['FRF']=self.allFRF
            dictOut['FRFgrad']=self.allFRFgrad
            dictOut['paraVal']=self.paraData['oldval']
        return dictOut


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def plotFRF(self,allData=False,fileOut=None):
        """
        ##################################################################
        # Function used to plot FRF
        ##################################################################
        """
        if self.caseProp['computeFRF']:            
            from plotFRF import plotFRF
            dictOut=self.getOutputDict(allData)
            plotFRF(dictIn=dictOut,fileOut=fileOut)            
        else:
            logging.info('Nothing to plot (FRF not computed)')


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
    
    def saveResults(self,):
        """
        Method used to save results
        """
        pass

    

