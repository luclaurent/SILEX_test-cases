import gmsh

gmsh.initialize()

c=gmsh.model.occ.addEllipse(0, 2, -1, 1,0.5)
r=gmsh.model.occ.rotate([(1,c)], 0,0,0, 0,1,0, 3.14/2)
b=gmsh.model.occ.extrude([(1,c)], 2, 0, 0)

gmsh.model.occ.addCylinder(0, 0, 0, 1, 0, 0, 1)
gmsh.model.occ.addRectangle(0.5, -0.5, -0.5, 3, 3)
gmsh.model.occ.synchronize()
gmsh.fltk.run()