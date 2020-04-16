#### debug SILEX class


import SILEXclass

# create SILEX object
PB = SILEXclass.SILEX()

#load case properties
caseDict = dict()
caseDict['name'] = 'test'
caseDict['freqMax'] = 100            # maximum frequency of the range
caseDict['freqMin'] = 0.1         # minimum frequency of the range
caseDict['nbSteps'] = 5          # number of frequency steps
caseDict['modal'] = False         # modal version of the computation (building of modal basis and use it for gradients)
caseDict['computeFRF'] = True     # computation of the FRF in control volume
caseDict['typeLS'] = ''           # type of Level-Set (FromMesh or manual)
PB.loadComputeProperties(caseDict)

#load mechanical properties
mechaProp = dict()
mechaProp['celerity'] = 340
mechaProp['rho'] = 1.2
mechaProp['fluid_damping'] = 1+0.01j
PB.loadMechaProperties(mechaProp)

#load data
dataLoad = dict()
dataLoad['geomFolder'] = 'geom'             # folder of geometry and meshes
dataLoad['resultsFolder'] = 'results'        # folder of results
#
dataLoad['originalFluidMeshFile'] = ''      # provided fluid mesh file
dataLoad['originalStructGeoFile'] = ''      # provided structure geometry file (parametric, gmsh format...)
dataLoad['originalStructMeshFile'] = ''     # provided structure mesh file
dataLoad['currentStructMeshFile'] = ''      # build structure mesh file (generated from geometry)
dataLoad['resultsFile'] = ''                # current results file basename
dataLoad['prefixResults'] = 'SILEX'         # prefix use on results files
#
dataLoad['exportData'] = 'mat'              # default format to export data (as FRF, mat or pickle)
dataLoad['exportMesh'] = 'msh'              # default format to export fields and mesh (msh, msh+, vtk)
PB.loadData(dataLoad)

#load boundary conditions
bcdef= dict()
bcdef['disp']={'type':'bbx','data':[0,0,0,0,0,0],'values':3.1250E-05}
PB.loadBC(bcdef)

#solve the problem
PB.solvePb([3.,3.,1.,1.])