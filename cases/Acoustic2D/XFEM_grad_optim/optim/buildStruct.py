import scipy

### Build levelset of the structure and gradients of it


def buildStruct(paraVal,fluidNodes,caseDefine):

    print('##### build Levelset (case: '+caseDefine+')')

    #depending on the chosen case
    if caseDefine == 'thick_u':

        #load values of parameters
        x_pos_struc=paraVal[0]
        y_pos_struc=paraVal[1]
        radius_hcircle=paraVal[2]
        angleStruct=paraVal[3]

        # create coordinates of nodes of the structure (half circle)
        nbNodesHC=50
        nbNodesSC=25
        nbNodesWall=1
        #
        angleU=angleStruct*scipy.pi/180
        thicknessU=1.
        # points
        xc2=x_pos_struc-radius_hcircle*scipy.cos(angleU)
        yc2=y_pos_struc-radius_hcircle*scipy.sin(angleU)
        xc4=x_pos_struc+radius_hcircle*scipy.cos(angleU)
        yc4=y_pos_struc+radius_hcircle*scipy.sin(angleU)
        #parameter for circles
        thetaHC=scipy.linspace(-angleU-scipy.pi/2,-angleU+scipy.pi/2,nbNodesHC)
        thetaSC=scipy.linspace(-angleU-scipy.pi/2,-angleU+scipy.pi/2,nbNodesSC)
        #inner large circle
        xNodesIHC=x_pos_struc-(radius_hcircle-thicknessU/2.)*scipy.sin(thetaHC)
        yNodesIHC=y_pos_struc-(radius_hcircle-thicknessU/2.)*scipy.cos(thetaHC)
        #outer larger circle
        xNodesOHC=x_pos_struc-(radius_hcircle+thicknessU/2.)*scipy.sin(thetaHC[::-1])
        yNodesOHC=y_pos_struc-(radius_hcircle+thicknessU/2.)*scipy.cos(thetaHC[::-1])
        #small circle 2
        xNodesC2=xc2-thicknessU/2*scipy.sin(thetaSC[::-1]+scipy.pi)
        yNodesC2=yc2-thicknessU/2*scipy.cos(thetaSC[::-1]+scipy.pi)
        #small circle 4
        xNodesC4=xc4+thicknessU/2*scipy.sin(thetaSC[::-1])
        yNodesC4=yc4+thicknessU/2*scipy.cos(thetaSC[::-1])

        strucNodes=scipy.vstack([scipy.hstack([xNodesIHC,xNodesC2[1:],xNodesOHC[1:],xNodesC4[1:]]),scipy.hstack([yNodesIHC,yNodesC2[1:],yNodesOHC[1:],yNodesC4[1:]])]).transpose()
        print(strucNodes)

        nbNodesAllStruct=2*nbNodesHC+2*nbNodesSC-4
        lCA=scipy.linspace(1,nbNodesAllStruct-1,nbNodesAllStruct-1)
        lCB=scipy.linspace(2,nbNodesAllStruct,nbNodesAllStruct-1)

        strucElem=scipy.vstack([lCA,lCB]).transpose()
        strucElem=scipy.vstack([strucElem,[1,nbNodesAllStruct]])
        print(strucElem)

        NbNodesFluid=fluidNodes.shape[0]
        LevelSet=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_X=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_Y=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_R=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_T=scipy.zeros(NbNodesFluid)
        IndicZone=scipy.zeros(NbNodesFluid)
        # level set 3 cases
        for itN in range(NbNodesFluid):
            xNcurr=fluidNodes[itN,0]
            yNcurr=fluidNodes[itN,1]
            #check areas
            if angleU != scipy.pi/2 or angleU != 3*scipy.pi/2:
                zoneA=yNcurr-scipy.tan(angleU)*(xNcurr-x_pos_struc)-y_pos_struc<0
            else:
                zoneA=xNcurr-x_pos_struc<0

            zoneAA=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle>=0
            zoneAB=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radius_hcircle<0
            if angleU != scipy.pi and angleU != 0. and angleU != 2*scipy.pi:
                zoneB=yNcurr-scipy.tan(angleU+scipy.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc>0
                zoneC=yNcurr-scipy.tan(angleU+scipy.pi/2)*(xNcurr-x_pos_struc)-y_pos_struc<0
            else:
                zoneB=xNcurr-x_pos_struc>0
                zoneC=xNcurr-x_pos_struc<0

            if zoneA:
                if zoneAA:
                    radiusOC=radius_hcircle+thicknessU/2.
                    LevelSet[itN]=scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusOC
                    LevelSet_gradient_X[itN]=-(xNcurr-x_pos_struc)/(LevelSet[itN]+radiusOC)
                    LevelSet_gradient_Y[itN]=-(yNcurr-y_pos_struc)/(LevelSet[itN]+radiusOC)
                    LevelSet_gradient_R[itN]=-1 
                    LevelSet_gradient_T[itN]=0 
                    IndicZone[itN]=1
                elif zoneAB: 
                    radiusIC=radius_hcircle-thicknessU/2.
                    LevelSet[itN]=-(scipy.sqrt((xNcurr-x_pos_struc)**2+(yNcurr-y_pos_struc)**2)-radiusIC)
                    LevelSet_gradient_X[itN]=(xNcurr-x_pos_struc)/(-LevelSet[itN]+radiusIC)
                    LevelSet_gradient_Y[itN]=(yNcurr-y_pos_struc)/(-LevelSet[itN]+radiusIC)
                    LevelSet_gradient_R[itN]=1 
                    LevelSet_gradient_T[itN]=0
                    IndicZone[itN]=2
            elif zoneC:
                radiusC2=thicknessU/2.
                LevelSet[itN]=scipy.sqrt((xNcurr-xc2)**2+(yNcurr-yc2)**2)-radiusC2
                LevelSet_gradient_X[itN]=-(xNcurr-xc2)/(LevelSet[itN]+radiusC2)
                LevelSet_gradient_Y[itN]=-(yNcurr-yc2)/(LevelSet[itN]+radiusC2)
                LevelSet_gradient_R[itN]=-scipy.cos(angleU)*LevelSet_gradient_X[itN]-scipy.sin(angleU)*LevelSet_gradient_Y[itN]
                LevelSet_gradient_T[itN]=radius_hcircle*scipy.sin(angleU)*LevelSet_gradient_X[itN]-radius_hcircle*scipy.cos(angleU)*LevelSet_gradient_Y[itN]
                IndicZone[itN]=3
            elif zoneB:
                radiusC4=thicknessU/2.
                LevelSet[itN]=scipy.sqrt((xNcurr-xc4)**2+(yNcurr-yc4)**2)-radiusC4
                LevelSet_gradient_X[itN]=-(xNcurr-xc4)/(LevelSet[itN]+radiusC4)
                LevelSet_gradient_Y[itN]=-(yNcurr-yc4)/(LevelSet[itN]+radiusC4)
                LevelSet_gradient_R[itN]=scipy.cos(angleU)*LevelSet_gradient_X[itN]+scipy.sin(angleU)*LevelSet_gradient_Y[itN]
                LevelSet_gradient_T[itN]=-radius_hcircle*scipy.sin(angleU)*LevelSet_gradient_X[itN]+radius_hcircle*scipy.cos(angleU)*LevelSet_gradient_Y[itN]
                IndicZone[itN]=4
            
        LevelSetTangent=fluidNodes[:,1]-max(fluidNodes[:,1])

        #store date
        LevelSetGradient=[LevelSet_gradient_X,LevelSet_gradient_Y,LevelSet_gradient_R,LevelSet_gradient_T]
        NamePara=['X','Y','R','T']

    elif caseDefine == 'u':
        pass

    elif caseDefine == 'thin_x_low_wall':
        pass

    elif caseDefine == 'thin_x_up_wall':
        pass

    elif caseDefine == 'thick_x_up_wall':

        #load values of parameters
        x_pos_struc=paraVal[0]
        y_pos_struc=paraVal[1]
        radius_hcircle=paraVal[2]

        # create coordinates of nodes of the structure (half circle)
        nbNodesHC=50
        nbNodesWall=1
        thetaHC=scipy.linspace(-scipy.pi/2,scipy.pi/2,nbNodesHC)
        
        xNodesHC=x_pos_struc-radius_hcircle*scipy.sin(thetaHC)
        yNodesHC=y_pos_struc-radius_hcircle*scipy.cos(thetaHC)

        strucNodes=scipy.vstack([xNodesHC,yNodesHC]).transpose()
        strucNodes=scipy.vstack([strucNodes,[x_pos_struc-radius_hcircle,max(fluidNodes[:,1])],[x_pos_struc+radius_hcircle,max(fluidNodes[:,1])]])
        print(strucNodes)

        lCA=scipy.linspace(1,nbNodesHC-1,nbNodesHC-1)
        lCB=scipy.linspace(2,nbNodesHC,nbNodesHC-1)


        strucElem=scipy.vstack([lCA,lCB]).transpose()
        strucElem=scipy.vstack([strucElem,[nbNodesHC,nbNodesHC+1],[1,nbNodesHC+2]])
        print(strucElem)

        NbNodesFluid=fluidNodes.shape[0]
        LevelSet=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_X=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_Y=scipy.zeros(NbNodesFluid)
        LevelSet_gradient_R=scipy.zeros(NbNodesFluid)
        IndicZone=scipy.zeros(NbNodesFluid)
        # level set 3 cases
        for itN in range(NbNodesFluid):
            if fluidNodes[itN,1]<y_pos_struc:
                LevelSet[itN]=scipy.sqrt((fluidNodes[itN,0]-x_pos_struc)**2+(fluidNodes[itN,1]-y_pos_struc)**2)-radius_hcircle
                LevelSet_gradient_X[itN]=-(fluidNodes[itN,0]-x_pos_struc)/(LevelSet[itN]+radius_hcircle)
                LevelSet_gradient_Y[itN]=-(fluidNodes[itN,1]-y_pos_struc)/(LevelSet[itN]+radius_hcircle)
                LevelSet_gradient_R[itN]=-1    
                IndicZone[itN]=1        
            elif fluidNodes[itN,0]<x_pos_struc:
                LevelSet[itN]=-(fluidNodes[itN,0]-x_pos_struc+radius_hcircle)
                LevelSet_gradient_X[itN]=1
                LevelSet_gradient_Y[itN]=0
                LevelSet_gradient_R[itN]=-1
                IndicZone[itN]=2
            elif fluidNodes[itN,0]>x_pos_struc:
                LevelSet[itN]=fluidNodes[itN,0]-x_pos_struc-radius_hcircle
                LevelSet_gradient_X[itN]=-1
                LevelSet_gradient_Y[itN]=0
                LevelSet_gradient_R[itN]=-1
                IndicZone[itN]=3


        LevelSetTangent=fluidNodes[:,1]-max(fluidNodes[:,1])

        #store date
        LevelSetGradient=[LevelSet_gradient_X,LevelSet_gradient_Y,LevelSet_gradient_R]
        NamePara=['X','Y','R']


    return strucNodes,strucElem,LevelSet,LevelSetGradient,NamePara,LevelSetTangent,IndicZone

