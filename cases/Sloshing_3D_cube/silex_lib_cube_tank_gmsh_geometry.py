import gmsh
import sys
from pathlib import Path

def xfem_fluid_and_tank(lx1,ly1,lz1,h,ElementOrder,file_name):


    #// Parameters: acoustic cavity
    #lx1 = 1.0;
    #ly1 = 0.8;
    #lz1 = 0.6;

    #// size of elements
    #h =  lx1/20;

    gmsh.initialize()
    gmsh.model.add('toto')

    #Mesh.CharacteristicLengthMax=10*h;
    #gmsh.model.mesh.setOrder(ElementOrder)
    gmsh.option.setNumber("Mesh.ElementOrder", ElementOrder)

    #// Cavity: Corners
    gmsh.model.geo.addPoint(0,     0  , 0, h,1   )
    gmsh.model.geo.addPoint(lx1,    0  , 0, h,2  )
    gmsh.model.geo.addPoint(lx1,    ly1 , 0, h,3 )
    gmsh.model.geo.addPoint(0 ,    ly1 , 0, h,4  )
    gmsh.model.geo.addPoint(0,     0  , lz1, h,5 )
    gmsh.model.geo.addPoint(lx1,    0  , lz1, h,6)
    gmsh.model.geo.addPoint(lx1,   ly1 , lz1, h,7)
    gmsh.model.geo.addPoint(0 ,    ly1 , lz1, h,8)


    #// Cavity: lines
    gmsh.model.geo.addLine(1, 5,1)
    gmsh.model.geo.addLine(5, 6,2)
    gmsh.model.geo.addLine(6, 2,3)
    gmsh.model.geo.addLine(2, 1,4)
    gmsh.model.geo.addLine(1, 4,5)
    gmsh.model.geo.addLine(4, 8,6)
    gmsh.model.geo.addLine(8, 5,7)
    gmsh.model.geo.addLine(8, 7,8)
    gmsh.model.geo.addLine(7, 3,9)
    gmsh.model.geo.addLine(3, 4,10)
    gmsh.model.geo.addLine(2, 3,11)
    gmsh.model.geo.addLine(6, 7,12)

    gmsh.model.geo.addCurveLoop([-7, 1, -5, -6],17)
    gmsh.model.geo.addSurfaceFilling([17],18)
    gmsh.model.geo.addCurveLoop([-2, -3, -4, -1],19)
    gmsh.model.geo.addSurfaceFilling([19],20)
    gmsh.model.geo.addCurveLoop([3, 11, -9, -12],21)
    gmsh.model.geo.addSurfaceFilling([21],22)
    gmsh.model.geo.addCurveLoop([-11, -10, 5, 4],23)
    gmsh.model.geo.addSurfaceFilling([23],24)
    gmsh.model.geo.addCurveLoop([2, 12, -8, 7],25)
    gmsh.model.geo.addSurfaceFilling([25],26)
    gmsh.model.geo.addCurveLoop([8, 9, 10, 6],27)
    gmsh.model.geo.addSurfaceFilling([27],28)


    gmsh.model.geo.addSurfaceLoop([18, 26, 20, 22, 24, 28],1025)

    gmsh.model.geo.addVolume([1025],1027)

    gmsh.model.addPhysicalGroup(3, [1027], 10, "fluid")

    #// Structure surface
    gmsh.model.addPhysicalGroup(2, [20, 22, 24, 28, 18], 20, "tank surface")

    #// Free fluid surface
    gmsh.model.addPhysicalGroup(2, [26],30,"free surface")
    gmsh.model.geo.synchronize() 
    gmsh.model.mesh.generate(3)

    #
    #gmsh.fltk.run()
    file_name_save=file_name

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'.msh')
    #gmsh.write(file_name_save.as_posix()+'.geo_unrolled')

    #gmsh.fltk.run()

    gmsh.finalize()
    return

def Stiffener_DKT_bis(lx_nominal,
                      ly_nominal,
                      lz_nominal,
                      lxashift_up_A,
                      lxashift_down_A,
                      lxashift_up_B,
                      lxashift_down_B,
                      lzashift_up_A,
                      lzashift_up_B,
                      h,ElementOrder,file_name):

    #// h : size of elements

    gmsh.initialize()
    gmsh.model.add('titi')
    gmsh.option.setNumber("Mesh.ElementOrder",ElementOrder)

    #// 
    gmsh.model.geo.addPoint(lx_nominal + lxashift_down_A,   0  , 0             , h,10   )
    gmsh.model.geo.addPoint(lx_nominal + lxashift_up_A,   0  , lz_nominal+lzashift_up_A           , h,12   )  
    gmsh.model.geo.addPoint(lx_nominal + lxashift_down_B,    ly_nominal, 0   , h,20)
    gmsh.model.geo.addPoint(lx_nominal + lxashift_up_B,    ly_nominal, lz_nominal+lzashift_up_B , h,22)

    gmsh.model.geo.addLine(10, 12,13)
    gmsh.model.geo.addLine(12, 22,14)
    gmsh.model.geo.addLine(22, 20,15)
    gmsh.model.geo.addLine(20, 10,16)

    gmsh.model.geo.addCurveLoop([16, 13, 14, 15],1)
    #gmsh.model.geo.addSurfaceFilling([1],1)
    gmsh.model.geo.addSurfaceFilling([1],1)

    gmsh.model.addPhysicalGroup(2, [1], 50, "Stiffener baffle surface")

    gmsh.model.addPhysicalGroup(1, [14], 60, "Stiffener baffle edge")

    gmsh.model.addPhysicalGroup(1, [16], 70, "Stiffener baffle base edge to impose acceleration, velocity or displacement")

    gmsh.model.geo.synchronize() 
    gmsh.model.mesh.generate(2)
    file_name_save=file_name
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'.msh')
    
    gmsh.finalize()

    return

def Stiffener_DKT(lx1,ly1,lz1,lxa,lxashift_up,lxashift_down,lza,h,ElementOrder,file_name):

    #// h : size of elements

    gmsh.initialize()
    gmsh.model.add('titi')
    gmsh.option.setNumber("Mesh.ElementOrder",ElementOrder)

    #// 
    gmsh.model.geo.addPoint(lxa,   0  , 0             , h,10   )
    gmsh.model.geo.addPoint(lxa,   0  , lza           , h,12   )  
    gmsh.model.geo.addPoint(lxa+lxashift_down,    ly1, 0   , h,20)
    gmsh.model.geo.addPoint(lxa+lxashift_up,    ly1, lza , h,22)

    gmsh.model.geo.addLine(10, 12,13)
    gmsh.model.geo.addLine(12, 22,14)
    gmsh.model.geo.addLine(22, 20,15)
    gmsh.model.geo.addLine(20, 10,16)

    gmsh.model.geo.addCurveLoop([16, 13, 14, 15],1)
    #gmsh.model.geo.addSurfaceFilling([1],1)
    gmsh.model.geo.addSurfaceFilling([1],1)

    gmsh.model.addPhysicalGroup(2, [1], 50, "Stiffener baffle surface")

    gmsh.model.addPhysicalGroup(1, [14], 60, "Stiffener baffle edge")

    gmsh.model.addPhysicalGroup(1, [16], 70, "Stiffener baffle base edge to impose acceleration, velocity or displacement")

    gmsh.model.geo.synchronize() 
    gmsh.model.mesh.generate(2)
    file_name_save=file_name
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'.msh')
    
    gmsh.finalize()

    return

def classic_fluid_and_tank(lx1,ly1,lz1,lxa,lxashift_up,lxashift_down,lza,th,h,ElementOrder,file_name):

    gmsh.initialize()
    gmsh.model.add('tata')

    #Mesh.CharacteristicLengthMax=10*h;
    #gmsh.model.mesh.setOrder(ElementOrder)
    gmsh.option.setNumber("Mesh.ElementOrder", ElementOrder)

    #// Cavity: Corners
    gmsh.model.geo.addPoint(0,     0  , 0, h,1   )
    gmsh.model.geo.addPoint(lx1,    0  , 0, h,2  )
    gmsh.model.geo.addPoint(lx1,    ly1 , 0, h,3 )
    gmsh.model.geo.addPoint(0 ,    ly1 , 0, h,4  )
    gmsh.model.geo.addPoint(0,     0  , lz1, h,5 )
    gmsh.model.geo.addPoint(lx1,    0  , lz1, h,6)
    gmsh.model.geo.addPoint(lx1,   ly1 , lz1, h,7)
    gmsh.model.geo.addPoint(0 ,    ly1 , lz1, h,8)

    gmsh.model.geo.addPoint(lxa-th/2,   0  , 0, h   ,10)
    gmsh.model.geo.addPoint(lxa+th/2,   0  , 0, h   ,11)
    gmsh.model.geo.addPoint(lxa-th/2,   0  , lza , h,12)
    gmsh.model.geo.addPoint(lxa+th/2,   0  , lza , h,13)

    gmsh.model.geo.addPoint(lxa+lxashift_down-th/2,   ly1  , 0, h   ,20)
    gmsh.model.geo.addPoint(lxa+lxashift_down+th/2,   ly1  , 0, h   ,21)
    gmsh.model.geo.addPoint(lxa+lxashift_up-th/2,   ly1  , lza , h,22)
    gmsh.model.geo.addPoint(lxa+lxashift_up+th/2,   ly1  , lza , h,23)


    #// Cavity: lines
    gmsh.model.geo.addLine(1, 5,1)
    gmsh.model.geo.addLine(5, 6,2)
    gmsh.model.geo.addLine(6, 2,3)
    #gmsh.model.geo.addLine(2, 1,4)
    gmsh.model.geo.addLine(1, 4,5)
    gmsh.model.geo.addLine(4, 8,6)
    gmsh.model.geo.addLine(8, 5,7)
    gmsh.model.geo.addLine(8, 7,8)
    gmsh.model.geo.addLine(7, 3,9)
    #gmsh.model.geo.addLine(3, 4,10)
    gmsh.model.geo.addLine(2, 3,11)
    gmsh.model.geo.addLine(6, 7,12)

    gmsh.model.geo.addLine(1, 10  , 13)
    gmsh.model.geo.addLine(10,  12 , 14)
    gmsh.model.geo.addLine(12,  13 , 15)
    gmsh.model.geo.addLine(13,  11 , 16)
    gmsh.model.geo.addLine(11,  2 , 17) 
    gmsh.model.geo.addLine(4, 20  , 18) 
    gmsh.model.geo.addLine(20,  22 , 19)
    gmsh.model.geo.addLine(22,  23 , 20)
    gmsh.model.geo.addLine(23,  21 , 21)
    gmsh.model.geo.addLine(21,  3 , 22) 
    gmsh.model.geo.addLine(10,  20 , 23)
    gmsh.model.geo.addLine(12,  22 , 24)
    gmsh.model.geo.addLine(13,  23 , 25)
    gmsh.model.geo.addLine(11,  21 , 26)


    gmsh.model.geo.addCurveLoop([5, 18, -23, -13 ], 1)
    gmsh.model.geo.addSurfaceFilling([1 ],1)
    
    gmsh.model.geo.addCurveLoop([ 23, 19, -24, -14], 2)
    gmsh.model.geo.addSurfaceFilling([2 ],2)

    gmsh.model.geo.addCurveLoop([-15, -25, 20, 24 ], 3)
    gmsh.model.geo.addSurfaceFilling([3 ],3)

    gmsh.model.geo.addCurveLoop([-16, -26, 21, 25 ], 4)
    gmsh.model.geo.addSurfaceFilling([4 ],4)

    gmsh.model.geo.addCurveLoop([ -17, -11, 22, 26], 5)
    gmsh.model.geo.addSurfaceFilling([ 5],5)

    gmsh.model.geo.addCurveLoop([11, -9, -12, 3 ], 6)
    gmsh.model.geo.addSurfaceFilling([6 ],6)

    gmsh.model.geo.addCurveLoop([ -5, -6, -7, 1],7 )
    gmsh.model.geo.addSurfaceFilling([ 7],7)

    gmsh.model.geo.addCurveLoop([ -18, -19, -20, -21, -22, 9, 8, 6],8 )
    gmsh.model.geo.addPlaneSurface([ 8],8)

    gmsh.model.geo.addCurveLoop([ 13, 14, 15, 16, 17, -3, -2, -1],9 )
    gmsh.model.geo.addPlaneSurface([ 9],9)

    gmsh.model.geo.addCurveLoop([ 2, 12, -8, 7], 10)
    gmsh.model.geo.addSurfaceFilling([ 10],10)
    
    gmsh.model.geo.addSurfaceLoop([7, 1, 8, 2, 3, 9, 4, 5, 6, 10],1)

    gmsh.model.geo.addVolume([1],1)

    gmsh.model.addPhysicalGroup(3, [1], 10, "fluid")

    #// Structure surface
    gmsh.model.addPhysicalGroup(2, [7, 1, 2, 3, 4, 5, 6, 9, 8], 20, "tank surface")

    #// Free fluid surface
    gmsh.model.addPhysicalGroup(2, [10],30,"free surface")
    gmsh.model.geo.synchronize() 
    #gmsh.fltk.run()
    gmsh.model.mesh.generate(3)

    #
    #gmsh.fltk.run()
    file_name_save=file_name

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'.msh')
    #gmsh.write(file_name_save.as_posix()+'.geo_unrolled')

    #gmsh.fltk.run()

    gmsh.finalize()
    return


#xfem_fluid_and_tank( 2.5 , 1.5 , 0.5 , 0.1 , 2 , 'titi')

