###########################################################
# Declaration of useful FE operators depending on type of elements and dimension of the considered problem
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import logging
import sys
#
import numpy as np
#
from SILEX import silex_lib_xfem_acou_tet4
from SILEX import silex_lib_tri3_acou

class buildFE:
    def __init__(self,dimPB=2,typeElem=None):
        """
        Contruct object depending on dimension and elements
        """
        self.dim = dimPB            # dimension of the problem
        self.typeElem = typeElem    # type of elements
        self.lib = None             # object associated to Fortran library
        # check if logging is already loaded
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(stream=sys.stdout, level=logging.INFO)
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
        Load Fortran's libraries depending on the dimension
        """
        if self.getTypeElem() == 'TRI3':
            self.lib = silex_lib_tri3_acou
            logging.info('SILEX TRI3 library loaded')
        elif self.getTypeElem() == 'TET4':
            self.lib = silex_lib_xfem_acou_tet4
            logging.info('SILEX TRI3 library loaded')

    def buildLevelSet(self,cavityNodes=None,immersedNodes=None,immersedElems=None):
        """
        Build Level-set from mesh
        """
        LevelSet,LevelSetU = self.lib.computelevelset(
                cavityNodes,        # nodes of the cavity (the level-set will be returned at these nodes)
                immersedNodes,      # nodes of the immersed volume (no compatibility required)
                immersedElems)      # elements of the immersed volume
        return LevelSet,LevelSetU   # level-set and unsigned level-set

    def getEnrichedElements(self,cavityElems=None,levelset=None):
        """
        Get the enriched elements from the level-set defined at nodes
        """
        EnrichedElems, EnrichedNbElems = self.lib.getenrichedelementsfromlevelset(
            cavityElems,    # elements of the cavity (defined with nodes number)
            levelset)       # nodal values defining the signed levelset
        return EnrichedElems, EnrichedNbElems   # enriched elements and number of enriched elements

    def getGlobalAcousticsMatrices(self,cavityNodes=None,cavityElems=None,valueCelerity=None,valueRho=None):
        """
        Get the global acoustics matrices defined via Compressed Sparse Column forms
        """
        IIf, JJf, Vffk, Vffm = self.lib.globalacousticmatrices(
            cavityElems,    # elements of the cavity
            cavityNodes,    # nodes coordinates of the cavity
            valueCelerity,  # celerity value
            valueRho)       # density value
        return IIf,JJf,Vffk,Vffm    # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix)

    def getEnrichedMatrices(self,cavityNodes=None,cavityElems=None,levelset=None,valueCelerity=None,valueRho=None):
        """
        Get the enriched matrices via XFEM defined via Compressed Sparse Column forms
        """
        IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm = self.lib.globalxfemacousticmatrices(
            cavityElems,    # elements of the cavity
            cavityNodes,    # nodes coordinates of the cavity
            levelset,       # level-set defined at cavity nodes
            valueCelerity,  # celerity value
            valueRho)       # density value
        return IIaa, JJaa, IIaf, JJaf, Vaak, Vaam, Vafk, Vafm # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix, a = enrichment,f = fluid)

    def getGradEnrichedMatrices(self,cavityNodes=None,cavityElems=None,levelset=None,levelsetGrad=None,valueCelerity=None,valueRho=None):
        """
        Get the sensitivities of enriched matrices via XFEM defined via Compressed Sparse Column forms
        """
        IIf, JJf, Vfak_gradient, Vfam_gradient =\
            silex_lib_xfem_acou_tet4.globalacousticgradientmatrices(cavityNodes,    # nodes of cavity
                                                                    cavityElems,    # elements of cavity
                                                                    levelset,       # level-set defined at cavity nodes
                                                                    valueCelerity,  # celerity value
                                                                    valueRho,       # density value
                                                                    levelsetGrad)   # sensitivity of the level-set
        return IIf, JJf, Vfak_gradient, Vfam_gradient # indices for rows and columns, and matrices values at corresponding indices (k = stiffness matrix, m = mass matrix, a = enrichment,f = fluid)

    def computeQuadraticPressure(self, nodesCavity=None, elemsControlVolume=None, uncorrectedField=None, enrichmentField=None, levelset=None):
        """
        Compute the quadratic pressure with an enrichment field depending on the type of variable
        """
        if np.iscomplexobj(uncorrectedField):   # check if input values are complex
             quadPress = self.lib.computexfemcomplexquadratiquepressure(                 
                    elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                    nodesCavity,            # all nodes of the cavity
                    uncorrectedField,       # field not corrected via XFEM enrichment
                    enrichmentField,        # enrichment part of the field from XFEM
                    levelset,               # level-set defined via nodal values
                    levelset*0.-1.0)        #
        else:
            logging.error('Float version not usable')
        
        return quadPress

    def computeGradQuadraticPressure(self, nodesCavity=None, elemsControlVolume=None, Field=None, gradField=None, levelset=None, levelsetgrad=None):
        """
        Compute sensitivity of the quadratic pressure with an enrichment field depending on the type of variable
        """
        if np.iscomplexobj(Field):   # check if input values are complex
             quadPress = self.lib.computegradientcomplexquadratiquepressure(                 
                    elemsControlVolume,     # elements defining the volume in which the quadratic pressure must be computed
                    nodesCavity,            # all nodes of the cavity
                    Field,                  # nodal values of the field (corrected version)
                    gradField,              # nodal values of sensitivity of the field (corrected version)
                    levelset)               # level-set defined via nodal values
        else:
            logging.error('Float version not usable')
        
        return quadPress