###########################################################
# Various tools to manage stuctures and levelset
# L. Laurent - 2020 - luc.laurent@lecnam.net
###########################################################


import os
import logging
import numpy as np


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# change parameters in geometry file (gmsh) and build the mesh
def buildStructMesh(fileOrig, destFile, paraVal,build=False):
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

    LevelSet=[]         # signed distance (Level-set)
    LevelSetU=[]        # unsigned distance
    LevelSetGrad=list()     # gradients of the Level-set
    paraName=[]         # name of the parameters


    def __init__(self,typeName,nodes,paraVal):
        """
        constructor of the object
        """
        #build the LS
        self.buildLS(typeName,nodes,paraVal)

    def buildLS(self,typeName,nodes,paraVal):
        """
        build the Level-set and associated gradients
        """
        if  typeName is '3D_sphere':
            # check data
            if nodes.shape[1]<3:
                print('ERROR in manualLS >> bad dimension of nodes array (%i)'%(nodes.shape[1]))
            if len(paraVal)!=4:
                print('ERROR in manualLS >> bad number of parameter (%i)'%(len(paraVal)))
            # LS from a simple analytic shape (sphere)
            lx3 = paraVal[0]  # 2.0 # Xc
            ly3 = paraVal[1]  # 2.0 # Yc
            lz3 = paraVal[2]  # 0.0 # YZ
            R = paraVal[3]  # 1.0 # sphere radius
            #
            self.namePara = ['X', 'Y', 'Z', 'R']
            #
            logging.debug("Parameters values")
            logging.debug("X ", lx3, " Y ", ly3, " Z ", lz3, " R ", R)
            # analytic LS
            self.LevelSet = np.sqrt((nodes[:, 0]-lx3)**2+(nodes[:, 1]-ly3)**2+(nodes[:, 2]-lz3)**2)-R
            self.LevelSetU=np.abs(LevelSet)
            # levelset gradients
            # Compute LS gradient according to Xc
            self.LevelSetGrad.append((lx3-nodes[:, 0])/(self.LevelSet+R))
            # Compute LS gradient according to Yc
            self.LevelSetGrad.append((ly3-nodes[:, 1])/(self.LevelSet+R))
            # Compute LS gradient according to Zc
            self.LevelSetGrad.append((lz3-nodes[:, 2])/(self.LevelSet+R))
            # Compute LS gradient according to R
            self.LevelSetGrad.append(nodes[:, 0]*0.-1.)
        else:
            print('Type %s does not exist'%(typeName))

    def exportParaName(self):
        """
        export the names of the parameters
        """
        return self.paraName

    def exportLS(self):
        """
        export the level-set nodal values
        """
        return self.LevelSet,self.LevelSetU

    def exportLSgrad(self,listExport=None):
        """
        export the required level-set gradients
        """
        if listExport is None:
            #export all data
            LevelSetgradExport=self.LevelSetGrad
            paraNameExport=self.paraName
        else:
            #export required data
            if len(listExport)>0:
                LevelSetgradExport=list()
                paraNameExport=list()
                for it in listExport:
                    LevelSetgradExport.append(self.LevelSetGrad[it])
                    paraNameExport.append(self.paraName[it])
            else:
                #export all data
                LevelSetgradExport=self.LevelSetGrad
                paraNameExport=self.paraName
        return LevelSetgradExport,paraNameExport

