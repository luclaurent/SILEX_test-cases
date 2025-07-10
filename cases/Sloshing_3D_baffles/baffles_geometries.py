import gmsh
from loguru import logger as log
import numpy as np

class create_baffle_geometry:
    
    def __init__(self, parameters=None, active_parameters=None):
        """
        Initialize the baffle geometry creation.
        
        Parameters:
        - parameters: Optional dictionary of parameters for the geometry.
        """
        self.parameters = dict()
        if not active_parameters:
            active_parameters = ['X', 'Y', 'Z', 'hy', 'hz']
        self.active_parameters = active_parameters
        self.default_parameters = { 'X': 0.45, 
                                   'Y': 0.8, 
                                   'Z': 0.6, 
                                   'hy': 0.1, 
                                   'hz': 0.1, 
                                   'density': 1e-1,
                                   'mesh_order': 1}        
        self.filter_parameters(parameters=parameters)
        
    def filter_parameters(self, parameters= None):
        if parameters is not None:
            if not isinstance(parameters, dict):
                if len(parameters) != len(self.active_parameters):
                    raise ValueError("Parameters must be a dictionary with keys matching active parameters.")
                params = dict()
                for i, key in enumerate(parameters):
                    if key in self.active_parameters:
                        params[key] = parameters[i]
                    else:
                        raise ValueError(f"Parameter '{key}' is not an active parameter.")

            else:
                params = parameters
            # update parameters
            self.parameters.update(params)
            log.info(f"Parameters for baffle geometry: {self.parameters}")
        
    def create_geometry(self, parameters=None):

        if parameters is not None:
            self.filter_parameters(parameters)
        # Initialize GMSH
        gmsh.initialize()
        gmsh.model.add("baffle")
        
        # load parameters
        X = self.parameters.get('X', self.default_parameters['X'])
        Y = self.parameters.get('Y', self.default_parameters['Y'])
        Z = self.parameters.get('Z', self.default_parameters['Z'])
        hy = self.parameters.get('hy', self.default_parameters['hy'])
        hz = self.parameters.get('hz', self.default_parameters['hz'])
        #
        density = self.parameters.get('density', self.default_parameters['density'])
        

        # Create a simple baffle geometry
        # Define points
        p = []
        p.append(gmsh.model.geo.addPoint(X, 0, 0, density))
        p.append(gmsh.model.geo.addPoint(X, Y, 0, density))
        p.append(gmsh.model.geo.addPoint(X, Y, Z, density))
        p.append(gmsh.model.geo.addPoint(X, 0, Z, density))
        #
        pi = list()
        pi.append(gmsh.model.geo.addPoint(X, hy, hz, density))
        pi.append(gmsh.model.geo.addPoint(X, Y-hy, hz, density))
        pi.append(gmsh.model.geo.addPoint(X, Y-hy, Z-hz, density))
        pi.append(gmsh.model.geo.addPoint(X, 0+hy, Z-hz, density))

        # Create lines
        le = gmsh.model.geo.addPolyline([*p, p[0]])
        li = gmsh.model.geo.addPolyline([*pi, pi[0]])
        loe = gmsh.model.geo.addCurveLoop([le])
        loi = gmsh.model.geo.addCurveLoop([li])


        # Create a surface
        surface = gmsh.model.geo.addPlaneSurface([loe, loi])

        # Synchronize the geometry
        gmsh.model.geo.synchronize()

        # Define physical groups
        gmsh.model.addPhysicalGroup(2, [surface], name="Baffle")
        
        log.info("Baffle geometry created with parameters: "
                 f"X={X}, Y={Y}, Z={Z}, hy={hy}, hz={hz}, density={density}")
        
        # meshing
        gmsh.option.setNumber('Mesh.ElementOrder', self.parameters.get('mesh_order', self.default_parameters['mesh_order']) )
        gmsh.model.mesh.generate(2)
        log.info("Mesh done")
        
    def getNodes(self):
        """
        Get the nodes of the created geometry.
        
        Returns:
        - List of nodes in the geometry.
        """
        nodes = gmsh.model.mesh.getNodes()
        return nodes[1].reshape(-1, len(nodes[0]))
    
    def getElements(self):
        """
        Get the elements of the created geometry.
        
        Returns:
        - List of elements in the geometry.
        """
        elements = gmsh.model.mesh.getElements()
        # seek for tri3
        ixTri3 = np.where(elements[0] == 2)[0]
        if len(ixTri3) > 0:
            ixTri3 = ixTri3[0]
        else:
            log.warning("No Tri3 elements found in the mesh.")
            return None
        tagsTri3 = elements[1][ixTri3]
        log.info(f"Found {len(tagsTri3)} Tri3 elements")
        # get connectivity
        connectivityTri3 = elements[2][ixTri3].reshape(-1, 3)
        #
        return connectivityTri3
        
    def show(self):
        """
        Display the created geometry.
        """
        gmsh.fltk.run()

baffle = create_baffle_geometry()
baffle.create_geometry()
baffle.getNodes()
baffle.getElements()
baffle.show()