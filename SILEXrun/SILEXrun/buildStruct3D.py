import logging
import numpy as np

class struct3D :

    def __init__(self,name):
        self.name = name #name of the geometric case        
        #
        self.BDparaname = {'3D_sphere':['X', 'Y', 'Z', 'R']}
        #
        self.paraName = self.getParaName()

    def getParaName(self):
        return self.BDparaname[self.name]
    
    def showParaVal(self,paraval):
        #
        logging.debug("Parameters values")
        logging.debug(' '.join('{}={}'.format(x,n) for x,n in zip(self.getParaName(),paraval)))
    
    def getLS(self,nodes,paraval):

        self.showParaVal(paraval)
        #
        LevelSet = []
        LevelSetTangent = None
        LevelSetU = []
        LevelSetGrad = list()
        
        if self.name == '3D_sphere':
            lx3 = paraval[0]  # 2.0 # Xc
            ly3 = paraval[1]  # 2.0 # Yc
            lz3 = paraval[2]  # 0.0 # YZ
            R = paraval[3]  # 1.0 # sphere radius
            
            # analytic LS
            LevelSet = np.sqrt((nodes[:, 0]-lx3)**2+(nodes[:, 1]-ly3)**2+(nodes[:, 2]-lz3)**2)-R
            LevelSetU = np.abs(LevelSet)
            LevelSetTangent = 0.0*LevelSetU-1.0
            # levelset gradients
            # Compute LS gradient according to Xc
            LevelSetGrad.append((lx3-nodes[:, 0])/(LevelSet+R))
            # Compute LS gradient according to Yc
            LevelSetGrad.append((ly3-nodes[:, 1])/(LevelSet+R))
            # Compute LS gradient according to Zc
            LevelSetGrad.append((lz3-nodes[:, 2])/(LevelSet+R))
            # Compute LS gradient according to R
            LevelSetGrad.append(nodes[:, 0]*0.-1.) 
        else:
            logging.error('Undefined geometry')
        return LevelSet,LevelSetTangent,LevelSetU,LevelSetGrad
