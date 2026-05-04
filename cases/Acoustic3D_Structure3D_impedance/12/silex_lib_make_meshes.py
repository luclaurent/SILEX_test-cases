import gmsh
import sys
from pathlib import Path
import numpy as np

def xfem_fluid_cavity(dataPb):


    #// Parameters: acoustic cavity
    lx1=dataPb['lx1'] # large length of cavity
    ly1=dataPb['ly1'] # small length of cavity
    lz1=dataPb['lz1'] # height of cavity
    ly2=dataPb['ly2']
    lx2=dataPb['lx2']
    lx5=dataPb['lx5']
    ly5=dataPb['ly5']
    lz5=dataPb['lz5']
    h5= dataPb['h5'] 
    h= dataPb['h']
    ElementOrder= dataPb['ElementOrder']

    gmsh.initialize()
    gmsh.model.add('toto')

    #Mesh.CharacteristicLengthMax=10*h;
    #gmsh.model.mesh.setOrder(ElementOrder)
    gmsh.option.setNumber("Mesh.ElementOrder", ElementOrder)

    # Cavity: Corners
    gmsh.model.geo.addPoint(0,     0  , 0, h,1)
    gmsh.model.geo.addPoint(lx1,    0  , 0, h,2)
    gmsh.model.geo.addPoint(lx1,    ly1 , 0, h,3)
    gmsh.model.geo.addPoint(0 ,    ly1 , 0, h,4)
    gmsh.model.geo.addPoint(0,     0  , lz1, h,5)
    gmsh.model.geo.addPoint(lx1,    0  , lz1, h,6)
    gmsh.model.geo.addPoint(lx1,    ly1 , lz1, h,7)
    gmsh.model.geo.addPoint(0 ,    ly1 , lz1, h,8)
    gmsh.model.geo.addPoint(0,     ly1+ly2  , 0, h,9)
    gmsh.model.geo.addPoint(lx2,     ly1+ly2  , 0, h,10)
    gmsh.model.geo.addPoint(lx2,     ly1  , 0, h,11)
    gmsh.model.geo.addPoint(0,     ly1+ly2  , lz1, h,12)
    gmsh.model.geo.addPoint(lx2,     ly1+ly2  , lz1, h,13)
    gmsh.model.geo.addPoint(lx2,     ly1  , lz1, h,14)

    # Cavity: lines
    gmsh.model.geo.addLine( 1,  4   , 1 )       
    gmsh.model.geo.addLine( 4,  9   , 2 )       
    gmsh.model.geo.addLine( 9,  10  , 3 )      
    gmsh.model.geo.addLine( 10,  11 , 4 )       
    gmsh.model.geo.addLine( 11,  3  , 5 )       
    gmsh.model.geo.addLine( 3,  2   , 6 )       
    gmsh.model.geo.addLine( 2,  1   , 7 )       
    gmsh.model.geo.addLine( 1,  5   , 8 )       
    gmsh.model.geo.addLine( 4,  8   , 9 )         
    gmsh.model.geo.addLine( 9,  12  , 10)      
    gmsh.model.geo.addLine( 10,  13 , 11)       
    gmsh.model.geo.addLine( 11,  14 , 12)        
    gmsh.model.geo.addLine( 3,  7   , 13)   
    gmsh.model.geo.addLine( 2,  6   , 14)   
    gmsh.model.geo.addLine( 5,  8   , 15)   
    gmsh.model.geo.addLine( 8,  12  , 16)    
    gmsh.model.geo.addLine( 12,  13 , 17)     
    gmsh.model.geo.addLine( 13,  14 , 18)     
    gmsh.model.geo.addLine( 14,  7  , 19)    
    gmsh.model.geo.addLine( 7,  6   , 20)   
    gmsh.model.geo.addLine( 6,  5   , 21)   

    # Cavity surfaces
    gmsh.model.geo.addCurveLoop([5, 13, -19, -12],115)
    gmsh.model.geo.addPlaneSurface([115],116)
    gmsh.model.geo.addCurveLoop([4, 12, -18, -11],117)
    gmsh.model.geo.addPlaneSurface([117],118)
    gmsh.model.geo.addCurveLoop([3, 11, -17, -10],119)
    gmsh.model.geo.addPlaneSurface([119],120)
    gmsh.model.geo.addCurveLoop([2, 10, -16, -9],121)
    gmsh.model.geo.addPlaneSurface([121],122)
    gmsh.model.geo.addCurveLoop([1, 9, -15, -8],123)
    gmsh.model.geo.addPlaneSurface([123],124)
    gmsh.model.geo.addCurveLoop([7, 8, -21, -14],125)
    gmsh.model.geo.addPlaneSurface([125],126)
    gmsh.model.geo.addCurveLoop([7,   1,  2,  3,  4,  5,  6],127)
    gmsh.model.geo.addPlaneSurface([127],156)
    gmsh.model.geo.addCurveLoop([21, 15, 16, 17, 18, 19, 20],130)
    gmsh.model.geo.addPlaneSurface([130],131)
    gmsh.model.geo.addCurveLoop([-6,13,20,-14],132)
    gmsh.model.geo.addPlaneSurface([132],133)

    ###gmsh.model.geo.addSurfaceLoop([126, 156, 124, 122, 120, 118, 116, 114, 131],157)

    # control volume : cube center lx5,ly5,lz5 / side h5
    gmsh.model.geo.addPoint(lx5-h5/2,     ly5-h5/2  , lz5-h5/2,   h/2.5 , 31)
    gmsh.model.geo.addPoint(lx5-h5/2,    ly5-h5/2   , lz5+h5/2,   h/2.5 , 32)
    gmsh.model.geo.addPoint(lx5+h5/2,     ly5-h5/2  , lz5-h5/2,   h/2.5 , 33)
    gmsh.model.geo.addPoint(lx5+h5/2,    ly5-h5/2   , lz5+h5/2,   h/2.5 , 34)
    gmsh.model.geo.addPoint(lx5-h5/2,     ly5+h5/2  , lz5-h5/2,   h/2.5 , 35)
    gmsh.model.geo.addPoint(lx5-h5/2,    ly5+h5/2   , lz5+h5/2,   h/2.5 , 36)
    gmsh.model.geo.addPoint(lx5+h5/2,     ly5+h5/2  , lz5-h5/2,   h/2.5 , 37)
    gmsh.model.geo.addPoint(lx5+h5/2,    ly5+h5/2   , lz5+h5/2,   h/2.5 , 38)

    gmsh.model.geo.addLine(33, 37 , 75 )
    gmsh.model.geo.addLine(37, 35 , 76 )
    gmsh.model.geo.addLine(35, 31 , 77 )
    gmsh.model.geo.addLine(31, 33 , 78 )
    gmsh.model.geo.addLine(33, 34 , 79 )
    gmsh.model.geo.addLine(37, 38 , 80 )
    gmsh.model.geo.addLine(38, 36 , 81 )
    gmsh.model.geo.addLine(36, 32 , 82 )
    gmsh.model.geo.addLine(32, 34 , 83 )
    gmsh.model.geo.addLine(34, 38 , 84 )
    gmsh.model.geo.addLine(35, 36 , 85 )
    gmsh.model.geo.addLine(31, 32 , 86 )
    
    gmsh.model.geo.addCurveLoop([75, 80, -84, -79],87)
    gmsh.model.geo.addPlaneSurface([87],88)
    gmsh.model.geo.addCurveLoop([76, 85, -81, -80],89)
    gmsh.model.geo.addPlaneSurface([89],90)
    gmsh.model.geo.addCurveLoop([77, 86, -82, -85],91)
    gmsh.model.geo.addPlaneSurface([91],92)
    gmsh.model.geo.addCurveLoop([78, 79, -83, -86],93)
    gmsh.model.geo.addPlaneSurface([93],94)
    gmsh.model.geo.addCurveLoop([75, 76, 77, 78],95)
    gmsh.model.geo.addPlaneSurface([95],96)
    gmsh.model.geo.addCurveLoop([84, 81, 82, 83],97)
    gmsh.model.geo.addPlaneSurface([97],98)

    gmsh.model.geo.addSurfaceLoop([94, 96, 88, 90, 92, 98],184)
    gmsh.model.geo.addVolume([184],185)

    # Acoustic cavity volume with control volume removed from it
    gmsh.model.geo.addSurfaceLoop([126, 156, 124, 122, 120, 118, 116 , 131, 133],186)
    gmsh.model.geo.addVolume([184,186],187)

    gmsh.model.addPhysicalGroup(3, [187], 10, "acoustic fluid with no control volume")
    gmsh.model.addPhysicalGroup(3, [187,185], 20, "acoustic fluid + control volume")
    gmsh.model.addPhysicalGroup(3, [185], 30, "control volume")
    gmsh.model.addPhysicalGroup(2, [131], 40, "impedance interface, roof")

    gmsh.model.geo.synchronize() 
    gmsh.model.mesh.generate(3)
    #
    #gmsh.fltk.run()
    file_name_save=dataPb['mesh_file']

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'.msh')
    gmsh.write(file_name_save.as_posix()+'.geo_unrolled')

    #gmsh.fltk.run()

    gmsh.finalize()
    return

###################################################
def structure(dataPb):


    ElementOrder= dataPb['ElementOrder_struc']

    lx3 = dataPb['lx3']
    ly3 = dataPb['ly3']
    l4 = dataPb['l4']
    lz4 = dataPb['lz4']
    r4 = dataPb['r4']
    deg = dataPb['deg']
    angle = deg*np.pi/180
    h= dataPb['h_struc']

    gmsh.initialize()
    gmsh.model.add('toto')

    #Mesh.CharacteristicLengthMax=10*h;
    #gmsh.model.mesh.setOrder(ElementOrder)
    gmsh.option.setNumber("Mesh.ElementOrder", ElementOrder)

    # Cavity: Corners
    gmsh.model.geo.addPoint(0,     0  , 0, h,1)
    gmsh.model.geo.addPoint(lx3,     ly3  , 0, h,15)
    gmsh.model.geo.addPoint(lx3+l4,     ly3  , 0, h,16)  
    gmsh.model.geo.addPoint(lx3,     ly3  , lz4-r4 , h,17)   
    gmsh.model.geo.addPoint(lx3+l4,     ly3  , lz4-r4, h,18)    
    gmsh.model.geo.addPoint(lx3+r4,     ly3  , lz4-r4 , h,19)   
    gmsh.model.geo.addPoint(lx3+l4-r4,     ly3  , lz4-r4, h,20) 
    gmsh.model.geo.addPoint(lx3+r4,     ly3  , lz4 , h,21)
    gmsh.model.geo.addPoint(lx3+l4-r4,     ly3  , lz4, h,22)


#// rotate wall
    gmsh.model.geo.rotate([(0,18),(0,20), (0,21), (0,19), (0,17), (0,16), (0,15), (0,22)],
lx3+l4/2 , ly3 , 0 ,     0,0,1,angle)
#Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
#  Point{ 18, 20, 21, 19, 17, 16, 15, 22};
#}
    # lines
    gmsh.model.geo.addLine( 15,  16   , 140 )       
    gmsh.model.geo.addLine( 16,18,141)
    gmsh.model.geo.addLine( 22,21,142)
    gmsh.model.geo.addLine( 17,15,143)

    gmsh.model.geo.addCircleArc(18,20,22)
    gmsh.model.geo.addCircleArc(21,19,17)
 
    gmsh.model.geo.addCurveLoop([141, 144, 142, 145, 143, 140],154)

    gmsh.model.geo.addPlaneSurface([154],155)

    gmsh.model.addPhysicalGroup(2, [155], 20, "structure")
    gmsh.model.addPhysicalGroup(1, [143, 145, 142, 144, 141], 30, "structure edge of the free boundary in air")
    gmsh.model.addPhysicalGroup(1, [140], 40, "fixed structure edge")

    gmsh.model.geo.synchronize() 
    gmsh.model.mesh.generate(3)
    #
    #gmsh.fltk.run()
    file_name_save=dataPb['mesh_file']

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
    gmsh.write(file_name_save.as_posix()+'_struc.msh')
    gmsh.write(file_name_save.as_posix()+'_struc.geo_unrolled')

    #gmsh.fltk.run()

    gmsh.finalize()
    return




