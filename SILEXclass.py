###########################################################
# AIR CAVITY
# 3D RIGID STRUCTURE
# XFEM
###########################################################

# python -m cProfile [-o output_file] [-s sort_order] myscript.py

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_CB_fluid_porous_7-v4.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#
###########################################################

###########################################################
# Libraries
###########################################################
import getopt
import importlib
from importlib import util as utilImport
import string
import time
import logging
from datetime import datetime 
import re
#
import numpy as np
import scipy as sp
#
import pickle
#
import sys
import os
from shutil import copyfile
#
import utils
import structTools
import pre
import post
import solver
import tools
#
from SILEX import silex_lib_xfem_acou_tet4
from SILEX import silex_lib_gmsh







###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# computation class
class SILEX(tools.tools,pre.preProcess,post.postProcess,solver.solverTools):
    #case properties
    caseProp = dict()
    caseProp['freqMax'] = []          # maximum frequency of the range
    caseProp['freqMin'] = []          # minimum frequency of the range
    caseProp['nbSteps'] = []          # number of frequency steps
    caseProp['modal'] = False         # modal version of the computation (building of modal basis and use it for gradients)
    caseProp['computeFRF'] = True     # computation of the FRF in control volume
    #
    caseProp['bcdisp'] = list()       # boundary conditions on displacement
    caseProp['bcpress'] = list()      # boundary conditions on pressure
    #
    caseProp['name'] = ''             # name of the case
    #
    caseProp['typeLS'] = ''           # type of Level-Set (FromMesh or manual)
    caseProp['typeGEOstruct'] = ''    # type of geometry of the structure (in the case of manual declaration (see structTools.py))

    #parameters values
    paraData = dict()
    paraData['oldval'] = list()       # previous values of parameters
    paraData['val'] = None            # current values of parameters
    paraData['name'] = None           # name of parameters
    paraData['nb'] = None             # number of parameters
    paraData['nameGrad'] = None       # name of parameters for gradients
    paraData['nbGrad'] = None         # number of gradients
    paraData['gradCompute'] = False   # compute gradients or not

    #material properties
    # air
    mechaProp = dict()
    mechaProp['celerity'] = 340.0
    mechaProp['rho'] = 1.2
    mechaProp['fluid_damping'] = (1+0.01j)

    #data properties
    data = dict()
    data['geomFolder'] = 'geom'             # folder of geometry and meshes
    data['resultsFolder'] = 'results'        # folder of results
    #
    data['originalFluidMeshFile'] = ''      # provided fluid mesh file
    data['originalStructGeoFile'] = ''      # provided structure geometry file (parametric, gmsh format...)
    data['originalStructMeshFile'] = ''     # provided structure mesh file
    data['currentStructMeshFile'] = ''      # build structure mesh file (generated from geometry)
    data['resultsFile'] = ''                # current results file basename
    data['prefixResults'] = 'SILEX'         # prefix use on results files
    #
    data['exportData'] = 'mat'              # default format to export data (as FRF, mat or pickle)
    data['exportMesh'] = 'vtk'              # default format to export fields and mesh (msh, mshv2, vtk)
    #
    fullPathCurrentResultsFolder = ''
    #
    saveResults = True                      # flag used to declare that the results must be saved
    classSave = None                        # object to store class for saving result
    debug = True                            # debug mode to write additionnal data

    #architecture properties
    commMPI = None
    rankMPI = 0
    nbProcMPI = 0
    #
    mumpsOk = False
    classMumps = None               # object containing the Mumps class (if available)

    #flags
    flags = dict()
    flags['saveResults'] = True  # flag to save results
    flags['edgeEnrichement'] = False  # flag to enrich the edge

    ###
    # storage variables
    LevelSet = []               # nodal values of LS (known at fluid nodes) signed
    LevelSetU = []              # nodal values of LS (known at fluid nodes) unsigned
    LevelSetGradient = list()   # nodal values of the parametric gradients of LS (known at fluid nodes)
    #
    fluidNodes = []             # coordinates of the fluid nodes
    fluidNbNodes = 0            # number of fluid nodes
    fluidNbDofs = 0             # number of dofs in fluid
    #
    fluidElems = []             # array of elements for the fluid volume
    idFluidNodes = []           # 
    fluidNbElems = 0            # number of elements in fluid volume
    fluidNbNodes = 0            # number of nodes in fluid volume
    #
    fluidElemsControl = []      # array of elements for the control volume
    idFluidNodesControl = []    # 
    fluidNbElemsControl = 0     # number of elements in control volume
    fluidNbNodesControl = 0     # number of nodes in control volume
    #
    structNodes = []            # coordinates of the structure nodes 
    structNbNodes = []          # number of structures nodes
    #
    structElems = []            # array of elements of structure
    idStructNodes = []          # 
    structNbElems = 0           # number of elements for structure
    structNbNodes = 0           # number of nodes of structure
    #enriched parts
    EnrichedNodes = []          # nodes associated to enriched elements
    EnrichedElems = []          # enriched elements
    EnrichedNbElems = 0        # number of enriched elements
    ##############
    #operators for fluid part
    KFF = []                    # fluid rigidity matrix
    MFF = []                    # fluid mass matrix
    dKFF = list()               # gradient matrices
    dMFF = list()               # gradient matrices
    SolvedDofF = []             # list of solved fluid dofs
    ##############
    #operators for enriched part
    KAA = []                    # rigidity matrix associated to full enriched dofs
    KFA = []                    # rigidity matrix associated to mixed enriched dofs
    MAA = []                    # mass matrix associated to full enriched dofs
    MFA = []                    # mass matrix associated to mixed enriched dofs
    dKAA = list()               # gradient matrices
    dKFA = list()               # gradient matrices
    dMAA = list()               # gradient matrices
    dMFA = list()               # gradient matrices
    SolvedDofA = []             # list of solved enriched dofs
    ##############
    #operators for the whole coupled problem
    K = []                      # rigidity matrix
    M = []                      # mass matrix
    dK = list()                 # gradient matrices
    dM = list()                 # gradients matrices
    UF = []                     # second member
    dUF = list()                # gradients of second member
    SolvedDof = []
    ##############
    #operators for results of the whole coupled problem
    pressure = []               # pressure field (corrected with enrichment field)
    pressureEnrichment = []    # enrichment part of the pressure field
    pressureUncorrect = []      # uncorrected part of the pressure field
    #
    pressureGrad = list()               # gradient of the pressure field (corrected with enrichment field)
    pressureEnrichmentGrad = list()    # gradient of the enrichment part of the pressure field
    pressureUncorrectGrad = list()      # gradient of the uncorrected part of the pressure field

    ##############
    #operators for FRF
    Frequencies = []            # frequencies for computation
    FRF = []                    # results of frequency response function
    FRFgrad = list()            # gradients of the FRF
    #
    allFRF = list()             #save all FRF along iterations on parameters
    allFRFgrad = list()         # save all FRF gradients along iterations on parameters

    ##############
    #others data
    nbRuns = 0                   # number of runs

    ##################################################################
    ##################################################################
    ##################################################################
    ##################################################################
    def __init__(self):
        """
        ##################################################################
        # Constructor of the class
        ##################################################################
        """        
        #initialize logging
        loggingFile = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")+"_SILEX.log"        
        logging.basicConfig(
            handlers=[
                logging.FileHandler(loggingFile),
                logging.StreamHandler(sys.stdout)
                ],
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
            #
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##################################################")
        logging.info("##               Load SILEX object              ##")
        #create symlink to the file
        fileLogLast = 'last.log'
        if os.path.islink(fileLogLast):
            os.unlink(fileLogLast)
        os.symlink(loggingFile,fileLogLast)

    
        


# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# #####################
# # function for dealing with options
# def manageOpt(argv, dV):
#     # load default values
#     freqMin = dV.freqMin
#     freqMax = dV.freqMax
#     nbStep = dV.nbStep
#     paraVal = sp.array(dV.paraVal)
#     gradCompute = sp.array(dV.gradCompute)
#     #caseDefine = dV.caseDef

#     # load info from MPI
#     nbProc, rank, comm = mpiInfo()
#     # load options
#     opts, args = getopt.getopt(argv, "p:s:F:f:hp:c:g:")
#     for opt, arg in opts:
#         if opt == "-s":
#             nbStep = int(arg)
#         elif opt == "-F":
#             freqMax = float(arg)
#         elif opt == "-f":
#             freqMin = float(arg)
#         elif opt == "-p":
#             tmp = sp.array(arg.split(','), dtype=sp.float32)
#             paraVal = tmp
#         elif opt == "-c":
#             caseDefine = str(arg)
#         elif opt == "-g":
#             tmp = sp.array(arg.split(','), dtype=sp.int32)
#             gradCompute = tmp
#         elif opt == "-h":
#             usage()
#             sys.exit()
#     # print chosen parameters
#     print("Number of processors: ", nbProc)
#     print("Parameters: ", paraVal)
#     print("Number of frequency steps: ", nbStep)
#     print("Maximum frequency: ", freqMax)
#     print("Minimum frequency: ", freqMin)
#     print("Components of grad: ", gradCompute)
#     #print ("Case: ",caseDefine)
#     it = 0
#     for itP in paraVal:
#         print('Parameter num '+str(it)+': '+str(itP))
#         it = it+1
#     print("\n\n")

#     # run computation
#     RunPb(freqMin, freqMax, nbStep, nbProc, rank, comm,
#           paraVal, gradCompute, 1)  # ,caseDefine)

# # usage definition


# def usage():
#     dV = defaultV
#     print("Usage: ", sys.argv[0], "-psFfhg [+arg]")
#     print("\t -p : input parameters (default value ", dV.nbProc, ")")
#     print("\t -s : number of steps in the frequency range (default value ", dV.nbStep, ")")
#     print("\t -F : maximum frequency (default value ", dV.freqMax, ")")
#     print("\t -f : minimum frequency (default value ", dV.freqMin, ")")
#     print("\t -g : Components of grad (default value ", dV.gradCompute, ")")

# # default values


# class defaultV:
#     freqMin = 10.0
#     freqMax = 150.0
#     nbStep = 1000
#     paraVal = [2., 2., 1., 1.]
#     gradCompute = [0, 1, 2, 3]
#     nbProc = 1
#     #caseDef= 'thick_u'


# # Run autonomous
# if __name__ == '__main__':
#     # run with options
#     dV = defaultV
#     manageOpt(sys.argv[1:], dV)
