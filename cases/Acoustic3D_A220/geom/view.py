from cabin import CabinInterior
import gmsh

gmsh.initialize()
gmsh.model.add("A220_view")

c = CabinInterior()
c.build()

# Geometry view in Gmsh
c.show_geo_gmsh()

# Mesh view in Gmsh (inside air example)
c.show_mesh_gmsh(mesh_size=0.25, inside_air=True, subtract_internal_solids=False)

# Geometry view in PyVista
c.show_geo_pyvista(surface_mesh_size=0.35)

# Mesh view in PyVista
c.show_mesh_pyvista(mesh_size=0.25, inside_air=True, subtract_internal_solids=False)

gmsh.finalize()