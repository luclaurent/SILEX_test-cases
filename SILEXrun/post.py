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
from .misc import utils
#
from SILEXlib import silex_lib_xfem,MeshField
from meshRW import msh, vtk

# activate logger
Logger = logging.getLogger(__name__)

class postProcess(object):

    def writeMeshField(self,
                       writerMethod=None,
                       uncorrectedField=None,
                       enrichmentField=None,
                       fileName=None,
                       fieldName=None):
        """
        ##################################################################
        # function used to export results and mesh via fortran classMeshField
        ##################################################################
        """
        data = [self.fluidNodes, self.fluidElems,
                self.LevelSet, uncorrectedField, enrichmentField]
        import pickle
        with open('debug_export_2d.pck', 'wb') as f1:
            pickle.dump(data, f1)
        # first initialization
        if not self.classSave:
            # initialized class for results export
            self.logger.info('Initialization of the writing class')
            self.classSave = MeshField.MeshField(
                                self.fluidNodes,
                                self.fluidElems,
                                self.LevelSet,
                                self.LevelSetTangent)
        # declare the writer
        self.classSave.setWriter(fileName, writerMethod)
        # append fields to the existing file(s)
        writeOk = True
        if uncorrectedField is not None and enrichmentField is not None:
            writeOk = True
            self.classSave.addField(
                uncorrectedField, enrichmentField, fieldName)
        # append field without correction
        if uncorrectedField is not None and enrichmentField is None:
            writeOk = True
            self.classSave.addFieldcorrected(uncorrectedField, fieldName)
        # create list of written files
        listFiles = None
        if writeOk:
            file = str(self.classSave.getOutputFilename().astype(
                'str')).replace(' ', '')
            basename = os.path.relpath(file, self.fullPathCurrentResultsFolder).replace(
                '.msh', '').replace('.vtk', '')
            filelist = [f for f in os.listdir(
                self.fullPathCurrentResultsFolder) if re.match(basename+'.*\.(msh|vtk)', f)]
            filelist.sort()
            return filelist

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def exportResults(self,
                      typeExport=None,
                      method='msh',
                      dictFields=None,
                      fileName=None):
        """
        ##################################################################
        # function used to export results and mesh
        ##################################################################
        """
        if typeExport is not None:
            self.logger.info(
                'Export results: type {}, method {}'.format(typeExport, method))
            tic = time.process_time()
            if typeExport == 'manuFields':
                # export many fields with details
                if dictFields is None:
                    self.logger.info('Fields are missing')
                else:
                    if method == 'msh':
                        if dictFields is dict:
                            dictFields = [dictFields]
                        # prepare data to export
                        dataExport = list()
                        for fields in dictFields:
                            dataOneField = list()
                            dataOneField.append(fields['field'])
                            if fields['type'] == "nodal":
                                dataOneField.append('nodal')
                            dataOneField.append(1)
                            dataOneField.append(fields['name'])
                            # display
                            self.logger.info(
                                'Prepare to write field: {}'.format(fields['name']))
                            # fill list to export
                            dataExport.append(dataOneField)
                        # export data
                        if fileName is not None:
                            fileName = utils.addExt(fileName, ".msh")
                            self.logger.info('Write: {}'.format(fileName))
                            #
                            if os.path.exists(fileName):
                                self.logger.warning(
                                    '>>> File {} will be OVERWRITTEN'.format(fileName))
                            #
                            msh.mshWriter(filename=fileName,
                                          nodes=self.fluidNodes,
                                          elems=self.fluidElems,
                                          fields=dataExport)
                            self.logger.info('File size: {}'.format(
                                utils.file_size(fileName)))
                        else:
                            self.logger.info(
                                'Unable to write file: filename is missing')
                    #
                    if method == 'pos':
                        # write a pos files (gmsh syntax) which include the discontinuities
                        # one file per field
                        if dictFields is dict:
                            dictFields = [dictFields]
                        # prepare data to export
                        for fields in dictFields:
                            fileName = utils.addExt(fields['filename'], '.pos')
                            self.logger.info('Prepare to write field: {} in {}'.format(
                                fields['name'], fileName))
                            # check if data are well shaped
                            dataLS = fields['levelsetmod']
                            dataUnc = fields['fielduncorrected']
                            dataEnr = fields['fieldenrichment']
                            #
                            silex_lib_xfem.makexfemposfilefreq(
                                self.fluidNodes,
                                self.fluidElems,
                                dataLS,
                                dataUnc,
                                dataEnr,
                                fileName,
                                fields['name'])
                    #
                    if method == 'mshv2' or method == 'vtk':
                        # write a pos files (gmsh syntax) which include the discontinuities
                        # one file per field
                        if dictFields is dict:
                            dictFields = [dictFields]
                        # prepare data to export
                        for fields in dictFields:
                            # adapt extension
                            if method == 'mshv2':
                                fileName = utils.addExt(fileName, '.msh')
                            if method == 'vtk':
                                fileName = utils.addExt(fileName, '.vtk')
                            self.logger.info('Prepare to write field: {} in {}'.format(
                                fields['name'], fileName))
                            # check if data are well shaped
                            dataUnc = fields['fielduncorrected']
                            dataEnr = fields['fieldenrichment']
                            #
                            funU = utils.fixShapeArray(
                                dataUnc, self.fluidNbNodes, 'uncorrected fields')
                            funE = utils.fixShapeArray(
                                dataEnr, self.fluidNbNodes, 'enrichment fields')
                            #
                            if os.path.exists(fileName):
                                self.logger.warning(
                                    '>>> File {} will be OVERWRITTEN'.format(fileName))
                            #
                            fileList = self.writeMeshField(
                                method,
                                funU(dataUnc),
                                funE(dataEnr),
                                fileName,
                                fields['name'])
                            #
                            for file in fileList:
                                self.logger.info('File size: {} ({})'.format(utils.file_size(
                                    os.path.join(self.fullPathCurrentResultsFolder, file)), file))

            fileOk = None
            if fileName is not None:
                fileOk = fileName

            if typeExport == 'cavitymesh':
                if not fileOk:
                    fileOk = self.getResultFile(
                        detPara=True, addTxt='_air_cavity_Mesh', ext='msh')
                # export mesh of the cavity
                if method == 'msh':
                    if os.path.exists(fileOk):
                        self.logger.warning(
                            '>>> File {} will be OVERWRITTEN'.format(fileOk))
                    #
                    silex_lib_gmsh.WriteResults(
                        fileOk, self.fluidNodes, self.fluidElems, 4)

            if typeExport == 'controlvolmesh':
                if not fileOk:
                    fileOk = self.getResultFile(
                        detPara=True, addTxt='_air_controlled_volume_Mesh', ext='msh')
                # export mesh of the control volume
                if method == 'msh':
                    if os.path.exists(fileOk):
                        self.logger.warning(
                            '>>> File {} will be OVERWRITTEN'.format(fileOk))
                    #
                    silex_lib_gmsh.WriteResults(
                        fileOk, self.fluidNodes, self.fluidElemsControl, 4)
                    #
            if typeExport == 'struct':
                if not fileOk:
                    fileOk = self.getResultFile(
                        detPara=True, addTxt='_struc_surface', ext='msh')
                # export 2D mesh of the structur
                if method == 'msh':
                    if os.path.exists(fileOk):
                        self.logger.warning(
                            '>>> File {} will be OVERWRITTEN'.format(fileOk))
                    #
                    silex_lib_gmsh.WriteResults2(
                        fileOk, self.structNodes, self.structElems, 2)
                    #
            if typeExport == 'levelset':
                if not fileOk:
                    fileOk = self.getResultFile(
                        detPara=True, addTxt='_LS_data', ext='msh')
                # export level-set and gradients of the level-set
                if method == 'msh':
                    # create list of data
                    dataW = list()
                    # level-set
                    dataW.append([[self.LevelSet], 'nodal', 1, 'Level set'])
                    # gradients of level-set
                    itP = 0
                    for iN in self.paraData['nameGrad']:
                        dataW.append([[self.LevelSetGradient[itP]],
                                     'nodal', 1, 'Level Set Grad '+iN])
                        itP = itP+1
                    #
                    if os.path.exists(fileOk):
                        self.logger.warning(
                            '>>> File {} will be OVERWRITTEN'.format(fileOk))
                    #
                    silex_lib_gmsh.WriteResults2(
                        fileOk, self.fluidNodes, self.fluidElems, 4, dataW)
                    #
            if typeExport == 'enrichedPart':
                if not fileOk:
                    fileOk = self.getResultFile(
                        detPara=True, addTxt='_LS_enriched_elements', ext='msh')
                if method == 'msh':
                    #
                    if os.path.exists(fileOk):
                        self.logger.warning(
                            '>>> File {} will be OVERWRITTEN'.format(fileOk))
                    #
                    silex_lib_gmsh.WriteResults2(
                        fileOk, self.fluidNodes, self.fluidElems[self.EnrichedElems], 4)
            self.logger.info(
                '++++++++++++++++++++ Done - {} s'.format(time.process_time()-tic))

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def exportFieldsOnePara(self, typeSave=None, paraName=False):
        """
        ##################################################################
        # Export fields in mesh file for one set of parameters
        ##################################################################
        """
        if typeSave is None:
            typeSave = self.data['exportMesh']
            self.logger.info('Use default format: {}'.format(typeSave))
        #
        # information concerning parameters
        txtPara = self.formatPara()
        # create database to export
        dataW = list()
        dataW.append({
            'field': np.real(self.pressure),
            'fielduncorrected': np.real(self.pressureUncorrect),
            'fieldenrichment': np.real(self.pressureEnrichment),
            'type': 'nodal',
            'name': 'Pressure (real) ('+txtPara+')'})
        dataW.append({
            'field': np.imag(self.pressure),
            'fielduncorrected': np.imag(self.pressureUncorrect),
            'fieldenrichment': np.imag(self.pressureEnrichment),
            'type': 'nodal',
            'name': 'Pressure (imaginary) ('+txtPara+')'})
        dataW.append({
            'field': np.absolute(self.pressure),
            'fielduncorrected': np.absolute(self.pressureUncorrect),
            'fieldenrichment': np.absolute(self.pressureEnrichment),
            'type': 'nodal',
            'name': 'Pressure (norm) ('+txtPara+')'})
        #
        if self.paraData['gradCompute']:
            for itG, txtP in enumerate(self.getNameGrad()):
                dataW.append({
                    'field': np.real(self.pressureGrad[itG]),
                    'fielduncorrected': np.real(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment': np.real(self.pressureEnrichmentGrad[itG]),
                    'type': 'nodal',
                    'name': 'Grad. '+txtP+' Pressure (real part) ('+txtPara+')'})
                dataW.append({
                    'field': np.imag(self.pressureGrad[itG]),
                    'fielduncorrected': np.imag(self.pressureUncorrectGrad[itG]),
                    'fieldenrichment': np.imag(self.pressureEnrichmentGrad[itG]),
                    'type': 'nodal',
                    'name': 'Grad. '+txtP+' Pressure (imaginary part) ('+txtPara+')'})
                dataW.append({
                    'field': utils.computeNormGradComplexField(self.pressure, self.pressureGrad[itG]),
                    'fielduncorrected': utils.computeNormGradComplexField(self.pressureUncorrect, self.pressureUncorrectGrad[itG]),
                    'fieldenrichment': utils.computeNormGradComplexField(self.pressureEnrichment, self.pressureEnrichmentGrad[itG]),
                    'type': 'nodal',
                    'name': 'Grad. '+txtP+' Pressure (gradient of norm) ('+txtPara+')'})
        # write the file
        self.exportResults(
            typeExport="manuFields",
            method=typeSave,
            dictFields=dataW,
            fileName=self.getResultFile(detPara=paraName, addTxt='results_fluid', ext=None))

        if self.debug:
            # write the discontinuties of the field in file
            # new database to export
            dataW = list()
            #
            dataW.append({'fielduncorrected': np.real(self.pressureUncorrect),
                          'fieldenrichment': np.real(self.pressureEnrichment),
                          'levelsetmod': self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='plus_real', ext=None),
                          'name': 'Pressure + (real) ('+txtPara+')'})
            dataW.append({'fielduncorrected': np.imag(self.pressureUncorrect),
                          'fieldenrichment': np.imag(self.pressureEnrichment),
                          'levelsetmod': self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='plus_imag', ext=None),
                          'name': 'Pressure + (imaginary) ('+txtPara+')'})
            dataW.append({'fielduncorrected': np.absolute(self.pressureUncorrect),
                          'fieldenrichment': np.absolute(self.pressureEnrichment),
                          'levelsetmod': self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='plus_abs', ext=None),
                          'name': 'Pressure + (norm) ('+txtPara+')'})
            #
            dataW.append({'fielduncorrected': np.real(self.pressureUncorrect),
                          'fieldenrichment': -np.real(self.pressureEnrichment),
                          'levelsetmod': -self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='moins_real', ext=None),
                          'name': 'Pressure - (real) ('+txtPara+')'})
            dataW.append({'fielduncorrected': np.imag(self.pressureUncorrect),
                          'fieldenrichment': -np.imag(self.pressureEnrichment),
                          'levelsetmod': -self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='moins_imag', ext=None),
                          'name': 'Pressure - (imaginary) ('+txtPara+')'})
            dataW.append({'fielduncorrected': np.absolute(self.pressureUncorrect),
                          'fieldenrichment': -np.absolute(self.pressureEnrichment),
                          'levelsetmod': -self.LevelSet,
                          'filename': self.getResultFile(detPara=paraName, addTxt='moins_abs', ext=None),
                          'name': 'Pressure - (norm) ('+txtPara+')'})
            #
            if self.paraData['gradCompute']:
                for itG, txtP in np.ndenumerate(self.paraData['nameGrad']):
                    dataW.append({'fielduncorrected': np.real(self.pressureUncorrectGrad[itG]),
                                  'fieldenrichment': np.real(self.pressureEnrichmentGrad[itG]),
                                  'levelsetmod': self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'plus_real', ext=None),
                                  'name': 'Pressure + (real) ('+txtPara+')'})
                    dataW.append({'fielduncorrected': np.imag(self.pressureUncorrectGrad[itG]),
                                  'fieldenrichment': np.imag(self.pressureEnrichmentGrad[itG]),
                                  'levelsetmod': self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'plus_imag', ext=None),
                                  'name': 'Pressure + (imaginary) ('+txtPara+')'})
                    #
                    dataW.append({'fielduncorrected': np.real(self.pressureUncorrectGrad[itG]),
                                  'fieldenrichment': -np.real(self.pressureEnrichmentGrad[itG]),
                                  'levelsetmod': -self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'moins_real', ext=None),
                                  'name': 'Pressure - (real) ('+txtPara+')'})
                    dataW.append({'fielduncorrected': np.imag(self.pressureUncorrectGrad[itG]),
                                  'fieldenrichment': -np.imag(self.pressureEnrichmentGrad[itG]),
                                  'levelsetmod': -self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'moins_imag', ext=None),
                                  'name': 'Pressure - (imaginary) ('+txtPara+')'})
                    # compute gradients of absolute values
                    gradAbsUncorrect = (np.real(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect) +
                                        np.imag(self.pressureUncorrectGrad[itG])*np.real(self.pressureUncorrect))/np.absolute(self.pressureUncorrect)
                    # deal with zeros values in enrichment field
                    gradAbsEnrichment = np.zeros(
                        [self.fluidNodes, self.caseProp['nbSteps']])
                    gradAbsEnrichment[self.SolvedDofA, :] = (np.real(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.real(self.pressureEnrichment[self.SolvedDofA, :]) +
                                                             np.imag(self.pressureEnrichmentGrad[self.SolvedDofA, :])*np.imag(self.pressureEnrichment[self.SolvedDofA, :]))/np.absolute(self.pressureEnrichment[self.SolvedDofA, :])
                    # remove inf value
                    IX = np.absolute(self.pressureUncorrect) == 0.
                    gradAbsUncorrect[IX] = 1
                    gradAbsEnrichment[IX] = 1
                    #
                    dataW.append({'fielduncorrected': self.gradAbsUncorrect[itG],
                                  'fieldenrichment': self.gradAbsEnrichment[itG],
                                  'levelsetmod': self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'plus_abs', ext=None),
                                  'name': 'Pressure + (norm) ('+txtPara+')'})
                    dataW.append({'fielduncorrected': self.gradAbsUncorrect[itG],
                                  'fieldenrichment': -self.gradAbsEnrichment[itG],
                                  'levelsetmod': -self.LevelSet,
                                  'filename': self.getResultFile(detPara=paraName, addTxt='_grad_'+str(itG)+'moins_abs', ext=None),
                                  'name': 'Pressure - (norm) ('+txtPara+')'})

            # write the file
            self.exportResults(typeExport='manuFields',
                               method='pos',
                               dictFields=dataW,
                               fileName=self.getResultFile(
                                        detPara=paraName,
                                        addTxt='results_fluid',
                                        ext=None))


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def saveFRF(self, typeSave=None, paraName=False, allData=False):
        """
        ##################################################################
        # Function used to solve FRF in mat file or pickle
        ##################################################################
        """
        if typeSave is None:
            typeSave = self.data['exportData']
            self.logger.info('Use default format: {}'.format(typeSave))
        if self.caseProp['computeFRF']:
            dictOut = self.getOutputDict(allData)
            #
            txt = ''
            if allData:
                txt = '_all'
            if typeSave == 'mat':
                # clean dictionary
                utils.replace_none(dictOut)
                # export data
                scipy.io.savemat(self.getResultFile(
                    detPara=paraName, addTxt='results'+txt, ext='mat'), mdict=dictOut)
                self.logger.info('Export data in %s' % (self.getResultFile(
                    detPara=paraName, addTxt='results'+txt, ext='mat')))
            if typeSave == 'pickle':
                # export data
                f = open(self.getResultFile(detPara=paraName,
                         addTxt='results'+txt, ext='pck'))
                pickle.dump(dictOut, f)
                f.close()
                self.logger.info('Export data in {}'.format(self.getResultFile(
                    detPara=paraName, addTxt='results'+txt, ext='pck')))
        else:
            self.logger.info('Nothing to export (FRF not computed)')


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


    def getOutputDict(self, allData=False):
        """
        ##################################################################
        # Function used to create an output dictionary
        ##################################################################
        """
        # build dictionary
        dictOut = self.paraData
        if allData:
            dictOut['frequencies'] = self.Frequencies
            dictOut['FRF'] = np.vstack(self.allFRF)
            dictOut['FRFgrad'] = self.allFRFgrad
            dictOut['paraVal'] = self.paraData['oldval']
        else:
            dictOut['frequencies'] = self.Frequencies
            dictOut['FRF'] = self.FRF
            dictOut['FRFgrad'] = self.FRFgrad
            dictOut['paraVal'] = self.paraData['val']

        return dictOut.copy()


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

    def plotFRF(self, allData=False, fileOut=None):
        """
        ##################################################################
        # Function used to plot FRF
        ##################################################################
        """
        if self.caseProp['computeFRF']:
            from plotFRF import plotFRF
            dictOut = self.getOutputDict(allData)
            plotFRF(dictIn=dictOut, fileOut=fileOut)
        else:
            self.logger.info('Nothing to plot (FRF not computed)')


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
