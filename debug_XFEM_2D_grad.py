# debug SILEX class
from SILEX import SILEXclass

# create SILEX object
PB = SILEXclass.SILEX()

PB.debug = False


##
# load case properties
caseBD = dict()
caseBD['name'] = 'test'
# dimension of geometry
caseBD['dim'] = 2
# maximum frequency of the range
caseBD['freqMax'] = 150
# minimum frequency of the range
caseBD['freqMin'] = 100
# number of frequency steps
caseBD['nbSteps'] = 20
# modal version of the computation
# (building of modal basis and use it for gradients)
caseBD['modal'] = False
# computation of the FRF in control volume
caseBD['computeFRF'] = True
# type of Level-Set (FromMesh or manual)
caseBD['typeLS'] = 'manual'
# type of geometry of the structure
# (in the case of manual declaration (see structTools.py))
caseBD['typeGEOstruct'] = '2D_thin_x_low_wall'
# IFS coupling
caseBD['coupling'] = False
# consider poromaterials or not
caseBD['poromaterial'] = False
# consider structure (if not coupling is deactivated aswell)
caseBD['structure'] = True
PB.loadComputeProperties(caseBD)

##
# parameters values
paraBD = dict()
paraBD['oldval'] = list()       # previous values of parameters
paraBD['val'] = []              # current values of parameters
paraBD['name'] = []             # name of parameters
paraBD['nameGrad'] = []         # name of parameters for gradients
paraBD['gradCompute'] = False   # compute gradients or not
PB.loadPara(paraBD)

##
# load mechanical properties
mechaBD = dict()
mechaBD['celerity'] = 340
mechaBD['rho'] = 1.2
mechaBD['fluid_damping'] = 0 #1+0.01j
mechaBD['pressRef'] = 20e-6
PB.loadMechaProperties(mechaBD)

##
# load data
dataBD = dict()
# folder of geometry and meshes
dataBD['geomFolder'] = 'cases/Acoustic2D/geom'
dataBD['resultsFolder'] = 'results'       # folder of results
#
# provided fluid mesh file
dataBD['originalFluidMeshFile'] = 'xfem2_fluid.msh'
# provided structure geometry file (parametric, gmsh format...)
dataBD['originalStructGeoFile'] = ''
dataBD['originalStructMeshFile'] = ''     # provided structure mesh file
# build structure mesh file (generated from geometry)
dataBD['currentStructMeshFile'] = ''
dataBD['resultsFile'] = ''                # current results file basename
dataBD['prefixResults'] = 'SILEX'         # prefix use on results files
#
# default format to export data (as FRF, mat or pickle)
dataBD['exportData'] = 'mat'
# default format to export fields and mesh (msh, msh+, vtk)
dataBD['exportMesh'] = 'mshv2'
PB.loadData(dataBD)


##
# load boundary conditions
bcdef = dict()
# bcdef['disp'] = {'type': 'bbx', 'data': [0, 0, 0, 1], 'values': 3.1250E-05}
bcdef['press'] = {'type': 'bbx', 'data': [0, 0, 0, 1], 'values': 1.0}
PB.loadBC(bcdef)

########################################
# solve the problem
PB.solvePb([[0.82, 3]])#,[0.8201, 3],[0.8202, 3]])  # ,[3.,3.0000000001,1.,1.]])

# raise
