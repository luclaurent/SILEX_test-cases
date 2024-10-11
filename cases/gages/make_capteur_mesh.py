def WriteGmshGeoTri3(filename,parameters):
    a=parameters[0];
    b=parameters[1];

    n1=parameters[2];
    n2=parameters[3];

    f=open(filename+'.geo','w')
    f.write('n1=%9.4g;\n' % n1)
    f.write('n2=%9.4g;\n' % n2)
    f.write('a=%9.4g;\n' % a)
    f.write('b=%9.4g;\n' % b)

    f.write('Point(1) = {0.0, 0.0, 0.0, n1};\n')
    f.write('Point(2) = {192.0, 0.0, 0.0, n1};\n')
    f.write('Point(3) = {192.0, 40.0, 0.0, n1};\n')
    f.write('Point(4) = {0.0, 40.0, 0.0, n1};\n')

    f.write('Point(5) = {58.5-b, 24.0+a, 0.0, n1};\n')
    f.write('Point(6) = {58.5-b, 16.0-a, 0.0, n1};\n')
    f.write('Point(7) = {133.5+b, 24.0+a, 0.0, n1};\n')
    f.write('Point(8) = {133.5+b, 16.0-a, 0.0, n1};\n')
    f.write('Point(9) = {47.5-b, 24.0+a, 0.0, n2};\n')
    f.write('Point(10) = {47.5-b, 16.0-a, 0.0, n2};\n')
    f.write('Point(11) = {144.5+b, 24.0+a, 0.0, n2};\n')
    f.write('Point(12) = {144.5+b, 16.0-a, 0.0, n2};\n')
    f.write('Point(13) = {58.5-b, 0.0, 0.0, n2};\n')
    f.write('Point(14) = {58.5, 40.0, 0.0, n2};\n')
    f.write('Point(15) = {133.5+b, 0.0, 0.0, n2};\n')
    f.write('Point(16) = {133.5+b, 40.0, 0.0, n2};\n')
    f.write('Line(1) = {1, 13};\n')
    f.write('Line(2) = {13, 15};\n')
    f.write('Line(3) = {15, 2};\n')
    f.write('Line(4) = {2, 3};\n')
    f.write('Line(5) = {3, 16};\n')
    f.write('Line(6) = {16, 14};\n')
    f.write('Line(7) = {14, 4};\n')
    f.write('Line(8) = {4, 1};\n')
    f.write('x = Sqrt(11*11-6*6);\n')
    f.write('Point(17) = {58.5+x-b, 24.0+6.0+a, 0.0, n2};\n')
    f.write('Point(18) = {58.5+x-b, 16.0-6.0-a, 0.0, n2};\n')
    f.write('Point(19) = {133.5-x+b, 24.0+6.0+a, 0.0, n2};\n')
    f.write('Point(20) = {133.5-x+b, 16.0-6.0-a, 0.0, n2};\n')
    f.write('Circle(9) = {9, 5, 17};\n')
    f.write('Circle(10) = {10, 6, 18};\n')
    f.write('Circle(11) = {20, 8, 12};\n')
    f.write('Circle(12) = {11, 7, 19};\n')
    f.write('Line(13) = {18, 20};\n')
    f.write('Line(14) = {12, 11};\n')
    f.write('Line(15) = {19, 17};\n')
    f.write('Line(16) = {9, 10};\n')
    f.write('Line Loop(17) = {7, 8, 1, 2, 3, 4, 5, 6};\n')
    f.write('Line Loop(18) = {9, -15, -12, -14, -11, -13, -10, -16};\n')
    f.write('Plane Surface(19) = {17, 18};\n')
    f.write('Physical Surface(1) = {19};\n')
    f.write('Physical Line(2) = {8};\n')
    f.write('Physical Line(3) = {4};\n')

    f.close()
    import os
    os.system('gmsh -2 %s.geo' % filename)

    return

##################################################
def WriteGmshGeoTri6(filename,parameters):
    a=parameters[0];
    b=parameters[1];

    n1=parameters[2];
    n2=parameters[3];

    f=open(filename+'.geo','w')
    f.write('n1=%9.4g;\n' % n1)
    f.write('n2=%9.4g;\n' % n2)
    f.write('a=%9.4g;\n' % a)
    f.write('b=%9.4g;\n' % b)
    f.write('Mesh.ElementOrder = 2;\n')

    f.write('Point(1) = {0.0, 0.0, 0.0, n1};\n')
    f.write('Point(2) = {192.0, 0.0, 0.0, n1};\n')
    f.write('Point(3) = {192.0, 40.0, 0.0, n1};\n')
    f.write('Point(4) = {0.0, 40.0, 0.0, n1};\n')

    f.write('Point(5) = {58.5-b, 24.0+a, 0.0, n1};\n')
    f.write('Point(6) = {58.5-b, 16.0-a, 0.0, n1};\n')
    f.write('Point(7) = {133.5+b, 24.0+a, 0.0, n1};\n')
    f.write('Point(8) = {133.5+b, 16.0-a, 0.0, n1};\n')
    f.write('Point(9) = {47.5-b, 24.0+a, 0.0, n2};\n')
    f.write('Point(10) = {47.5-b, 16.0-a, 0.0, n2};\n')
    f.write('Point(11) = {144.5+b, 24.0+a, 0.0, n2};\n')
    f.write('Point(12) = {144.5+b, 16.0-a, 0.0, n2};\n')
    f.write('Point(13) = {58.5-b, 0.0, 0.0, n2};\n')
    f.write('Point(14) = {58.5, 40.0, 0.0, n2};\n')
    f.write('Point(15) = {133.5+b, 0.0, 0.0, n2};\n')
    f.write('Point(16) = {133.5+b, 40.0, 0.0, n2};\n')
    f.write('Line(1) = {1, 13};\n')
    f.write('Line(2) = {13, 15};\n')
    f.write('Line(3) = {15, 2};\n')
    f.write('Line(4) = {2, 3};\n')
    f.write('Line(5) = {3, 16};\n')
    f.write('Line(6) = {16, 14};\n')
    f.write('Line(7) = {14, 4};\n')
    f.write('Line(8) = {4, 1};\n')
    f.write('x = Sqrt(11*11-6*6);\n')
    f.write('Point(17) = {58.5+x-b, 24.0+6.0+a, 0.0, n2};\n')
    f.write('Point(18) = {58.5+x-b, 16.0-6.0-a, 0.0, n2};\n')
    f.write('Point(19) = {133.5-x+b, 24.0+6.0+a, 0.0, n2};\n')
    f.write('Point(20) = {133.5-x+b, 16.0-6.0-a, 0.0, n2};\n')
    f.write('Circle(9) = {9, 5, 17};\n')
    f.write('Circle(10) = {10, 6, 18};\n')
    f.write('Circle(11) = {20, 8, 12};\n')
    f.write('Circle(12) = {11, 7, 19};\n')
    f.write('Line(13) = {18, 20};\n')
    f.write('Line(14) = {12, 11};\n')
    f.write('Line(15) = {19, 17};\n')
    f.write('Line(16) = {9, 10};\n')
    f.write('Line Loop(17) = {7, 8, 1, 2, 3, 4, 5, 6};\n')
    f.write('Line Loop(18) = {9, -15, -12, -14, -11, -13, -10, -16};\n')
    f.write('Plane Surface(19) = {17, 18};\n')
    f.write('Physical Surface(1) = {19};\n')
    f.write('Physical Line(2) = {8};\n')
    f.write('Physical Line(3) = {4};\n')

    f.close()
    import os
    os.system('gmsh -2 %s.geo' % filename)

    return

############################
#import os
#WriteGmshGeo('capteur',[-2.0,20.0,5.0,1.0])




############################
#import os
#WriteGmshGeo('capteur',[-2.0,20.0,5.0,1.0])



