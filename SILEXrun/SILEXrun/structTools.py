###########################################################
# Various tools to manage stuctures and levelset
# L. Laurent - 2020 - luc.laurent@lecnam.net
###########################################################


import os
import logging
import numpy as np
#
from .buildStruct3D import struct3D
from .buildStruct2D import struct2D
#
Logger = logging.getLogger(__name__)

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


# change parameters in geometry file (gmsh) and build the mesh
def buildStructMesh(fileOrig, destFile, paraVal, build=False):
    # copy original file to the used one
    copyfile(fileOrig+'.geo', destFile+'.geo')
    # change value of parameters in the new file
    for key, value in enumerate(paraVal):
        oldText = "<val##"+str(key)+">"
        newText = '%g' % value
        # print(oldText)
        # print(newText)
        cmdSed = "sed -i 's/"+oldText+"/"+newText+"/g' "+destFile+'.geo'
        # print(cmdSed)
        os.system(cmdSed)

    # run gmsh to build the mesh
    if build:
        os.system('gmsh -3 -format msh2 '+destFile+'.geo')


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
#class to build a Level-set and associated gradients manually (simple topologies)
class LSmanual:    

    def __init__(self,typeName,nodes=None,paraVal=None):
        """
        constructor of the object
        """
        #initialize data in the class
        self.clearData()
        #build the LS
        self.buildLS(typeName,nodes,paraVal)

    def buildLS(self,typeName = None,nodes = None,paraVal = None):
        """
        build the Level-set and associated gradients
        """
        if typeName is not None:
            self.updateType(typeName)
        if nodes is not None:
            self.nodes = nodes
        if paraVal is not None:
            self.paraVal = paraVal
        #
        if typeName[0] == '2':
            self.objBuild = struct2D(typeName)
            self.dim = 2
        elif typeName[0] == '3':
            self.objBuild = struct3D(typeName)
            self.dim = 3
        else:
            Logger.error('Bad type: {} (must start with 2 or 3)'.format(typeName))

        # get the parameters names
        self.paraName = self.objBuild.getParaName()
        #
        if self.nodes is not None and self.paraVal is not None:
            # check data
            if self.nodes.shape[1] != self.dim:
                Logger.error('(manualLS) >> bad dimension of nodes array ({})'.format(self.nodes.shape[1]))
            if len(paraVal)!=len(self.paraName):
                Logger.error('(manualLS) >> bad number of parameters ({}) - expected: {}'.format(len(self.paraVal),len(self.paraName)))
            #build LS
            self.LevelSet,self.LevelSetTangent,self.LevelSetU,self.LevelSetGrad,self.structNodes,self.structElems,self.structBounds = self.objBuild.getLS(self.nodes,paraVal)
    
    def clearData(self):
        """
        clear all data in the class
        """
        self.LevelSet = None         # signed distance (Level-set)
        self.LevelSetTangent = None  # tangent Level-set
        self.LevelSetU = None        # unsigned distance
        self.LevelSetGrad =list()    # gradients of the Level-set
        self.paraName = None         # name of the parameters
        self.typeName = None         # name of the case
        self.nodes = None            # nodes
        self.paraVal = None          # value of the parameters
    
    def updateType(self,typeName):
        """
        update the name of the type
        """
        self.clearData()
        self.typeName = 'typeName'

    def exportParaName(self):
        """
        export the names of the parameters
        """
        return self.paraName

    def exportLS(self):
        """
        export the level-set nodal values
        """
        if self.LevelSet is None:
            self.buildLS()
        return self.LevelSet,self.LevelSetU
    
    def exportLST(self):
        """
        export the tangent level-set nodal values
        """
        if self.LevelSet is None:
            self.buildLS()
        return self.LevelSetTangent

    def exportLSgrad(self,listExport=None):
        """
        export the required level-set gradients
        """
        if self.LevelSet is None:
            self.buildLS()
        if listExport is None:
            #export all data
            LevelSetgradExport=self.LevelSetGrad
            paraNameExport=self.paraName
        else:
            #export required data
            if len(listExport)>0:
                LevelSetgradExport=list()
                paraNameExport=list()
                for it,v in enumerate(listExport):
                    LevelSetgradExport.append(self.LevelSetGrad[it])
                    paraNameExport.append(self.paraName[it])
            else:
                #export all data
                LevelSetgradExport=self.LevelSetGrad
                paraNameExport=self.paraName
        return LevelSetgradExport,paraNameExport

