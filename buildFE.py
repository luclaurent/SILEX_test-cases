###########################################################
# Declaration of useful FE operators depending on type of elements and dimension of the considered problem
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import logging
import sys
#
import numpy as np
#
from SILEXlib.lib import silex_lib_xfem


# activate logger
Logger = logging.getLogger(__name__)


class buildFE:
    def __init__(self,dimPB=2,typeElem=None):
        """
        Contruct object depending on dimension and elements
        """
        self.dim = dimPB            # dimension of the problem
        self.typeElem = typeElem    # type of elements
        self.lib = None             # object associated to Fortran library
        #fix dimension when type type of elements is loaded
        if self.typeElem is not None:
            if self.typeElem is 'TRI3':
                self.dim = 2
            elif self.typeElem is 'TET4':
                self.dim = 3
        #
        self.loadLib()
    
    def getDim(self):
        """
        Get the loaded dimension
        """
        return self.dim

    def getTypeElem(self):
        """
        Get the considered type of element
        """
        typeOut = self.typeElem
        if typeOut is None:
            if self.getDim() == 2:
                typeOut = 'TRI3'
            elif self.getDim() == 3:
                typeOut = 'TET4'
        return typeOut
    
    def loadLib(self):
        """
        Load Fortran's libraries depending on the type of elements
        """
        print(self.getTypeElem())
        self.lib = silex_lib_xfem.SILEXlibXFEM(typeElem=self.getTypeElem())

    def buildLevelSet(self,
        cavityNodes=None,
        immersedNodes=None,
        immersedElems=None):
        """
        Build Level-set from mesh
        """
        LevelSet,LevelSetU = self.lib.ComputeLevelSet(
                cavityNodes,        # nodes of the cavity (the level-set will be returned at these nodes)
                immersedNodes,      # nodes of the immersed volume (no compatibility required)
                immersedElems)      # elements of the immersed volume
        return LevelSet,LevelSetU   # level-set and unsigned level-set

    def getEnrichedElements(self,
                            cavityNodes=None,
                            cavityElems=None,
                            structNodes=None,
                            structElems=None,
                            levelset=None):
        """
        Get the enriched elements from the level-set defined at nodes
        """
        EnrichedElems, EnrichedNbElems = self.lib.getEnrichedElements(
            cavityNodes,
            cavityElems,    # elements of the cavity (defined with nodes number)
            structNodes,
            structElems,
            levelset)       # nodal values defining the signed levelset
        return EnrichedElems, EnrichedNbElems   # enriched elements and number of enriched elements

    def getGlobalAcousticsMatrices(self,
        cavityNodes=None,
        cavityElems=None,
        valueCelerity=None,
        valueRho=None):
        """
        Get the global acoustics matrices defined via Compressed Sparse Column forms
        """
        IIf, JJf, Vffk, Vffm = self.lib.getGlobalAcousticMatrices(
            cavityNodes,    # nodes coordinates of the cavity
            cavityElems,    # elements of the cavity            
            valueCelerity,  # celerity value
            valueRho)       # density value
        return IIf,JJf,Vffk,Vffm    # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix)

    def getEnrichedMatrices(self,
        cavityNodes=None,
        cavityElems=None,
        levelset=None,
        levelsetTg=None,
        valueCelerity=None,
        valueRho=None,
        flagEdge=1):
        """
        Get the enriched matrices via XFEM defined via Compressed Sparse Column forms
        """
        IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = self.lib.getGlobalXFEMAcousticMatrices(
            cavityNodes,    # nodes coordinates of the cavity
            cavityElems,    # elements of the cavity            
            levelset,       # level-set defined at cavity nodes
            levelsetTg,     # tangent level-set defined at cavity nodes
            valueCelerity,  # celerity value
            valueRho,       # density value
            flagEdge)       # flag edge enrichment (0=inactive, 1=active)
        return IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix, a = enrichment,f = fluid)

    def getGradEnrichedMatrices(self,
        cavityNodes=None,
        cavityElems=None,
        levelset=None,
        levelsetGrad=None,
        valueCelerity=None,
        valueRho=None):
        """
        Get the sensitivities of enriched matrices via XFEM defined via Compressed Sparse Column forms
        """
        IIf, JJf, Vfak_gradient, Vfam_gradient =\
            self.lib.getGlobalAcousticGradientMatrices(cavityNodes,    # nodes of cavity
                                                    cavityElems,    # elements of cavity
                                                    levelset,       # level-set defined at cavity nodes
                                                    levelsetGrad,   # sensitivity of the level-set
                                                    valueCelerity,  # celerity value
                                                    valueRho)       # density value
                                                    
        return IIf, JJf, Vfak_gradient, Vfam_gradient # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix, a = enrichment,f = fluid)

    def computeQuadraticPressure(self,
        nodesCavity=None,
        elemsControlVolume=None,
        uncorrectedField=None,
        enrichmentField=None,
        levelset=None):
        """
        Compute the quadratic pressure with an enrichment field depending on the type of variable
        """
        if enrichmentField is None:
            quadPress = self.lib.getQuadraticPressure(
                nodesCavity,            # all nodes of the cavity
                elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                uncorrectedField        # field not corrected via XFEM enrichment
            )
        else:
            quadPress = self.lib.getXfemQuadraticPressure(
                nodesCavity,            # all nodes of the cavity
                elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                uncorrectedField,       # field not corrected via XFEM enrichment
                enrichmentField,        # enrichment part of the field from XFEM
                levelset,               # level-set defined via nodal values
                levelset*0.-1.0,        # tangent level-set
                1                       # flag for modified Heaviside on the edge
            )
        
        return quadPress

    def computeGradQuadraticPressure(self,
        nodesCavity=None,
        elemsControlVolume=None,
        Field=None,
        gradField=None,
        enrichmentField=None,
        enrichmentGradField=None,
        levelset=None):
        """
        Compute sensitivity of the quadratic pressure with an enrichment field depending on the type of variable
        """
        if enrichmentField is None:
            quadPress = self.lib.getGradQuadraticPressure(
                nodesCavity,            # all nodes of the cavity
                elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                levelset,               # level-set defined via nodal values
                Field,                  # field
                gradField               # gradient of the field
            )
        else:
            quadPress = self.lib.getGradXfemQuadraticPressure(
                nodesCavity,            # all nodes of the cavity
                elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                levelset,               # level-set defined via nodal values
                levelset*0.-1.0,        # tangent level-set
                Field,                  # field not corrected via XFEM enrichment
                enrichmentField,        # enrichment part of the field from XFEM
                gradField,              # gradient of the field
                enrichmentGradField,    # gradient of the enrichement field
                1                       # flag for modified Heaviside on the edge
            )
        
        return quadPress
