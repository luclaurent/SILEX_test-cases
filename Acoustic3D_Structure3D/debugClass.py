#### debug SILEX class


import SILEXclass

# create SILEX object
PB = SILEXclass.SILEX()

#load case properties
caseBD = dict()
caseBD['name'] = 'test'
caseBD['freqMax'] = 100            # maximum frequency of the range
caseBD['freqMin'] = 0.1         # minimum frequency of the range
caseBD['nbSteps'] = 5          # number of frequency steps
caseBD['modal'] = False         # modal version of the computation (building of modal basis and use it for gradients)
caseBD['computeFRF'] = True     # computation of the FRF in control volume
caseBD['typeLS'] = 'manual'           # type of Level-Set (FromMesh or manual)
caseBD['typeGEOstruct'] = '3D_sphere'   # type of geometry of the structure (in the case of manual declaration (see structTools.py))
PB.loadComputeProperties(caseBD)

#parameters values
paraBD = dict()
paraBD['oldval'] = list()       # previous values of parameters
paraBD['val'] = []              # current values of parameters
paraBD['name'] = []             # name of parameters
paraBD['nb'] = []               # number of parameters
paraBD['nameGrad'] = []         # name of parameters for gradients
paraBD['nbGrad'] = []           # number of gradients
paraBD['gradCompute'] = False   # compute gradients or not
PB.loadPara(paraBD)

#load mechanical properties
mechaBD = dict()
mechaBD['celerity'] = 340
mechaBD['rho'] = 1.2
mechaBD['fluid_damping'] = 1+0.01j
PB.loadMechaProperties(mechaBD)

#load data
dataBD = dict()
dataBD['geomFolder'] = 'geom'             # folder of geometry and meshes
dataBD['resultsFolder'] = 'results'        # folder of results
#
dataBD['originalFluidMeshFile'] = 'cavity_acou3D_struc_3D_v3_air.msh'      # provided fluid mesh file
dataBD['originalStructGeoFile'] = ''      # provided structure geometry file (parametric, gmsh format...)
dataBD['originalStructMeshFile'] = ''     # provided structure mesh file
dataBD['currentStructMeshFile'] = ''      # build structure mesh file (generated from geometry)
dataBD['resultsFile'] = ''                # current results file basename
dataBD['prefixResults'] = 'SILEX'         # prefix use on results files
#
dataBD['exportData'] = 'mat'              # default format to export data (as FRF, mat or pickle)
dataBD['exportMesh'] = 'msh'              # default format to export fields and mesh (msh, msh+, vtk)
PB.loadData(dataBD)



#load boundary conditions
bcdef= dict()
bcdef['disp']={'type':'bbx','data':[0,0,5,5,0,0],'values':3.1250E-05}
PB.loadBC(bcdef)

#solve the problem
PB.solvePb([3.,3.,1.,1.])

raise