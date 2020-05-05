import logging
import numpy as np

class struct2D :

    def __init__(self,name):
        self.name = name #name of the geometric case
        self.paraName = self.getParaName()
        #
        self.BDparaname = {
            '2D_thick_u':['X','Y','R','T'],
            '2D_thin_x_low_wall':None,
            '2D_thin_x_up_wall':None,
            '2D_thick_x_up_wall':['X','Y','R']
            }

    def getParaName(self):
        return self.BDparaname[self.name]
    
    def showParaVal(self,paraval):
        #
        logging.debug("Parameters values")
        logging.debug(' '.join('{}={}'.format(x,n) for x,n in zip(self.getParaName(paraval))))
    
    def getLS(self,nodes,paraval):

        self.showParaVal(paraval)
        #
        LevelSet = []
        LevelSetU = []
        LevelSetGrad = list()
        
        if self.name == '2D_thick_u':
            #load values of parameters
            x_pos_struc=paraval[0]
            y_pos_struc=paraval[1]
            radius_hcircle=paraval[2]
            angleStruct=paraval[3]

            # create coordinates of nodes of the structure (half circle)
            nbNodesHC=50
            nbNodesSC=25
            nbNodesWall=1
            #
            angleU=angleStruct*np.pi/180
            thicknessU=1.
            # points
            xc2=x_pos_struc-radius_hcircle*np.cos(angleU)
            yc2=y_pos_struc-radius_hcircle*np.sin(angleU)
            xc4=x_pos_struc+radius_hcircle*np.cos(angleU)
            yc4=y_pos_struc+radius_hcircle*np.sin(angleU)
            #parameter for circles
            thetaHC=np.linspace(-angleU-np.pi/2,-angleU+np.pi/2,nbNodesHC)
            thetaSC=np.linspace(-angleU-np.pi/2,-angleU+np.pi/2,nbNodesSC)
            #inner large circle
            xNodesIHC=x_pos_struc-(radius_hcircle-thicknessU/2.)*np.sin(thetaHC)
            yNodesIHC=y_pos_struc-(radius_hcircle-thicknessU/2.)*np.cos(thetaHC)
            #outer larger circle
            xNodesOHC=x_pos_struc-(radius_hcircle+thicknessU/2.)*np.sin(thetaHC[::-1])
            yNodesOHC=y_pos_struc-(radius_hcircle+thicknessU/2.)*np.cos(thetaHC[::-1])
            #small circle 2
            xNodesC2=xc2-thicknessU/2*np.sin(thetaSC[::-1]+np.pi)
            yNodesC2=yc2-thicknessU/2*np.cos(thetaSC[::-1]+np.pi)
            #small circle 4
            xNodesC4=xc4+thicknessU/2*np.sin(thetaSC[::-1])
            yNodesC4=yc4+thicknessU/2*np.cos(thetaSC[::-1])
            #nodes list
            strucNodes=np.vstack([np.hstack([xNodesIHC,xNodesC2[1:],xNodesOHC[1:],xNodesC4[1:]]),np.hstack([yNodesIHC,yNodesC2[1:],yNodesOHC[1:],yNodesC4[1:]])]).transpose()
            # print(strucNodes)
            # total number of nodes
            nbNodesAllStruct=2*nbNodesHC+2*nbNodesSC-4
            #build list of elements
            lCA=np.linspace(1,nbNodesAllStruct-1,nbNodesAllStruct-1)
            lCB=np.linspace(2,nbNodesAllStruct,nbNodesAllStruct-1)
            #
            strucElem=np.vstack([lCA,lCB]).transpose()
            strucElem=np.vstack([strucElem,[1,nbNodesAllStruct]])
            # print(strucElem)
            #build level-sets
            NbNodesFluid=nodes.shape[0]
            LevelSet=np.zeros(NbNodesFluid)
            LevelSet_gradient_X=np.zeros(NbNodesFluid)
            LevelSet_gradient_Y=np.zeros(NbNodesFluid)
            LevelSet_gradient_R=np.zeros(NbNodesFluid)
            LevelSet_gradient_T=np.zeros(NbNodesFluid)
            IndicZone=np.zeros(NbNodesFluid)
            # level set 3 cases
            for itN in range(NbNodesFluid):
                xNcurr=nodes[itN,0]
                yNcurr=nodes[itN,1]
                #check areas
                if angleU != np.pi/2 or angleU != 3*np.pi/2:
                    zoneA=yNcurr-np.tan(angleU)*(xNcurr-x_pos_struc)-y_pos_struc<0
                else:
                    zoneA=xNcurr-x_pos_struc<0

                zoneAA=np.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle>=0
                zoneAB=np.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle<0
                if angleU != np.pi and angleU != 0. and angleU != 2*np.pi:
                    zoneB=yNcurr-np.tan(angleU+np.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc>0
                    zoneC=yNcurr-np.tan(angleU+np.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc<0
                else:
                    zoneB=xNcurr-x_pos_struc>0
                    zoneC=xNcurr-x_pos_struc<0

                if zoneA:
                    if zoneAA:
                        radiusOC=radius_hcircle+thicknessU/2.
                        LevelSet[itN]=np.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusOC
                        LevelSet_gradient_X[itN]=-(xNcurr-x_pos_struc)/(LevelSet[itN]+radiusOC)
                        LevelSet_gradient_Y[itN]=-(yNcurr-y_pos_struc)/(LevelSet[itN]+radiusOC)
                        LevelSet_gradient_R[itN]=-1 
                        LevelSet_gradient_T[itN]=0 
                        IndicZone[itN]=1
                    elif zoneAB: 
                        radiusIC=radius_hcircle-thicknessU/2.
                        LevelSet[itN]=-(np.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusIC)
                        LevelSet_gradient_X[itN]=(xNcurr-x_pos_struc)/(-LevelSet[itN]+radiusIC)
                        LevelSet_gradient_Y[itN]=(yNcurr-y_pos_struc)/(-LevelSet[itN]+radiusIC)
                        LevelSet_gradient_R[itN]=1 
                        LevelSet_gradient_T[itN]=0
                        IndicZone[itN]=2
                elif zoneC:
                    radiusC2=thicknessU/2.
                    LevelSet[itN]=np.sqrt((xNcurr-xc2)**2+(yNcurr-yc2)**2)-radiusC2
                    LevelSet_gradient_X[itN]=-(xNcurr-xc2)/(LevelSet[itN]+radiusC2)
                    LevelSet_gradient_Y[itN]=-(yNcurr-yc2)/(LevelSet[itN]+radiusC2)
                    LevelSet_gradient_R[itN]=0
                    LevelSet_gradient_T[itN]=-radius_hcircle*np.sin(angleU)*LevelSet_gradient_X[itN]+radius_hcircle*np.cos(angleU)*LevelSet_gradient_Y[itN]
                    IndicZone[itN]=3
                elif zoneB:
                    radiusC4=thicknessU/2.
                    LevelSet[itN]=np.sqrt((xNcurr-xc4)**2+(yNcurr-yc4)**2)-radiusC4
                    LevelSet_gradient_X[itN]=-(xNcurr-xc4)/(LevelSet[itN]+radiusC4)
                    LevelSet_gradient_Y[itN]=-(yNcurr-yc4)/(LevelSet[itN]+radiusC4)
                    LevelSet_gradient_R[itN]=0
                    LevelSet_gradient_T[itN]=radius_hcircle*np.sin(angleU)*LevelSet_gradient_X[itN]-radius_hcircle*np.cos(angleU)*LevelSet_gradient_Y[itN]
                    IndicZone[itN]=4
            #defined level-set tangent
            LevelSetTangent=nodes[:,1]-max(nodes[:,1])

            #store data
            LevelSetGradient=[LevelSet_gradient_X,LevelSet_gradient_Y,LevelSet_gradient_R,LevelSet_gradient_T]
        
        elif self.name == 'u':
            pass

        elif self.name == '2D_thin_x_low_wall':
            pass

        elif self.name == '2D_thin_x_up_wall':
            pass

        elif self.name == '2D_thick_x_up_wall':

            #load values of parameters
            x_pos_struc=paraval[0]
            y_pos_struc=paraval[1]
            radius_hcircle=paraval[2]

            # create coordinates of nodes of the structure (half circle)
            nbNodesHC=50
            nbNodesWall=1
            thetaHC=np.linspace(-np.pi/2,np.pi/2,nbNodesHC)
            
            xNodesHC=x_pos_struc-radius_hcircle*np.sin(thetaHC)
            yNodesHC=y_pos_struc-radius_hcircle*np.cos(thetaHC)

            strucNodes=np.vstack([xNodesHC,yNodesHC]).transpose()
            strucNodes=np.vstack([strucNodes,[x_pos_struc-radius_hcircle,max(nodes[:,1])],[x_pos_struc+radius_hcircle,max(nodes[:,1])]])
            # print(strucNodes)

            lCA=np.linspace(1,nbNodesHC-1,nbNodesHC-1)
            lCB=np.linspace(2,nbNodesHC,nbNodesHC-1)


            strucElem=np.vstack([lCA,lCB]).transpose()
            strucElem=np.vstack([strucElem,[nbNodesHC,nbNodesHC+1],[1,nbNodesHC+2]])
            # print(strucElem)

            NbNodesFluid=nodes.shape[0]
            LevelSet=np.zeros(NbNodesFluid)
            LevelSet_gradient_X=np.zeros(NbNodesFluid)
            LevelSet_gradient_Y=np.zeros(NbNodesFluid)
            LevelSet_gradient_R=np.zeros(NbNodesFluid)
            IndicZone=np.zeros(NbNodesFluid)
            # level set 3 cases
            for itN in range(NbNodesFluid):
                if nodes[itN,1]<y_pos_struc:
                    LevelSet[itN]=np.sqrt((nodes[itN,0]-x_pos_struc)**2+(nodes[itN,1]-y_pos_struc)**2)-radius_hcircle
                    LevelSet_gradient_X[itN]=-(nodes[itN,0]-x_pos_struc)/(LevelSet[itN]+radius_hcircle)
                    LevelSet_gradient_Y[itN]=-(nodes[itN,1]-y_pos_struc)/(LevelSet[itN]+radius_hcircle)
                    LevelSet_gradient_R[itN]=-1    
                    IndicZone[itN]=1        
                elif nodes[itN,0]<x_pos_struc:
                    LevelSet[itN]=-(nodes[itN,0]-x_pos_struc+radius_hcircle)
                    LevelSet_gradient_X[itN]=1
                    LevelSet_gradient_Y[itN]=0
                    LevelSet_gradient_R[itN]=-1
                    IndicZone[itN]=2
                elif nodes[itN,0]>x_pos_struc:
                    LevelSet[itN]=nodes[itN,0]-x_pos_struc-radius_hcircle
                    LevelSet_gradient_X[itN]=-1
                    LevelSet_gradient_Y[itN]=0
                    LevelSet_gradient_R[itN]=-1
                    IndicZone[itN]=3


            LevelSetTangent=nodes[:,1]-max(nodes[:,1])

            #store data
            LevelSetGradient=[LevelSet_gradient_X,LevelSet_gradient_Y,LevelSet_gradient_R]
        else:
            logging.error('Undefined geometry')
        #
        return LevelSet,LevelSetTangent,LevelSetU,LevelSetGrad
