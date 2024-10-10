import sys, getopt, time, os

#print(sys.argv[0])
#print('Argument List:', str(sys.argv))

angle=float(sys.argv[1])
positionX=float(sys.argv[2])

#parameters=[float(angle)]
#print("parameters=",parameters)
#print(type(parameters[0]))

if not os.path.exists('geom'):
    os.makedirs('geom')



filename='geom/cavity7_with_porous'

f=open(filename+'.geo','w')


f.write('// Parameters: acoustic cavity\n')
f.write('lx1 = 7.0;\n')
f.write('ly1 = 4.0;\n')
f.write('lz1 = 2.5;\n')

f.write('ly2 = 1.0;\n')
f.write('lx2 = 3.5;\n')

f.write('lx3 = %9.4g;\n' % positionX)
f.write('ly3 = 2.0;\n')
f.write('l4 = 2.5;\n')
f.write('lz4 = 2.0;\n')
f.write('r4   = 0.5;\n')
f.write('e4   = 0.2;\n')
f.write('deg=%9.4g;\n' % angle)
f.write('angle = deg*Pi/180;\n')


f.write('// porous\n')
f.write('ee = 0.20;\n')

f.write('// size of elements\n')
f.write('h = lx1/20;\n')
f.write('h2 = ee/1;\n')


f.write('lx5 = 6.0;\n')
f.write('ly5 = 1.0;\n')
f.write('lz5 = 1.5;\n')
f.write('h5  = 1.0;\n')

f.write('Mesh.CharacteristicLengthMax=1.5*h;\n')
f.write('Mesh.ElementOrder = 1;\n')

f.write('// Cavity: Corners\n')
f.write('Point(1) = {0,     0  , 0, h};\n')
f.write('Point(2) = {lx1,    0  , 0, h};\n')
f.write('Point(3) = {lx1,    ly1 , 0, h};\n')
f.write('Point(4) = {0 ,    ly1 , 0, h};\n')

f.write('Point(5) = {0,     0  , lz1, h2};\n')
f.write('Point(6) = {lx1,    0  , lz1, h2};\n')
f.write('Point(7) = {lx1,    ly1 , lz1, h2};\n')
f.write('Point(8) = {0 ,    ly1 , lz1, h2};\n')

f.write('Point(9) = {0,     ly1+ly2  , 0, h};\n')
f.write('Point(10) = {lx2,     ly1+ly2  , 0, h};\n')
f.write('Point(11) = {lx2,     ly1  , 0, h};\n')
f.write('Point(12) = {0,     ly1+ly2  , lz1, h2};\n')
f.write('Point(13) = {lx2,     ly1+ly2  , lz1, h2};\n')
f.write('Point(14) = {lx2,     ly1  , lz1, h2};\n')

f.write('Point(100) = {lx2/2,     ly1+ly2  , 0, h};\n')

f.write('// Cavity: lines\n')
f.write('Line(1) = {1, 4};\n')
f.write('Line(2) = {4, 9};\n')
f.write('Line(3) = {9, 100};\n')
f.write('Line(4) = {10, 11};\n')
f.write('Line(5) = {11, 3};\n')
f.write('Line(6) = {3, 2};\n')
f.write('Line(7) = {2, 1};\n')
f.write('Line(8) = {1, 5};\n')
f.write('Line(9) = {4, 8};\n')
f.write('Line(10) = {9, 12};\n')
f.write('Line(11) = {10, 13};\n')
f.write('Line(12) = {11, 14};\n')
f.write('Line(13) = {3, 7};\n')
f.write('Line(14) = {2, 6};\n')
f.write('Line(15) = {5, 8};\n')
f.write('Line(16) = {8, 12};\n')
f.write('Line(17) = {12, 13};\n')
f.write('Line(18) = {13, 14};\n')
f.write('Line(19) = {14, 7};\n')
f.write('Line(20) = {7, 6};\n')
f.write('Line(21) = {6, 5};\n')

f.write('Line(200)={100,10};\n')


f.write('// wall: corners\n')
f.write('Point(15) = {lx3,     ly3-e4/2  , 0, h/3};\n')
f.write('Point(16) = {lx3+l4,     ly3-e4/2  , 0, h/3};\n')
f.write('Point(17) = {lx3,     ly3-e4/2  , lz4-r4 , h/3};\n')
f.write('Point(18) = {lx3+l4,     ly3-e4/2  , lz4-r4, h/3};\n')
f.write('Point(19) = {lx3+r4,     ly3-e4/2  , lz4-r4 , h/3};\n')
f.write('Point(20) = {lx3+l4-r4,     ly3-e4/2  , lz4-r4, h/3};\n')
f.write('Point(21) = {lx3+r4,     ly3-e4/2  , lz4 , (h2+h)/2};\n')
f.write('Point(22) = {lx3+l4-r4,     ly3-e4/2  , lz4, (h2+h)/2};\n')

f.write('Point(23) = {lx3,     ly3+e4/2  , 0, h/3};\n')
f.write('Point(24) = {lx3+l4,     ly3+e4/2  , 0, h/3};\n')
f.write('Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h/3};\n')
f.write('Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h/3};\n')
f.write('Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h/3};\n')
f.write('Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h/3};\n')
f.write('Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , (h2+h)/2};\n')
f.write('Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, (h2+h)/2};\n')

f.write('// rotate wall\n')
f.write('Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {\n')
f.write('  Point{22, 30, 18, 26, 20, 21, 28, 29, 19, 27, 17, 25, 16, 24, 15, 23};\n')
f.write('}\n')

f.write('// wall : lines\n')
f.write('Line(22) = {15, 16};\n')
f.write('Line(23) = {16, 18};\n')
f.write('Line(24) = {22, 21};\n')
f.write('Line(25) = {17, 15};\n')
f.write('Circle(26) = {18, 20, 22};\n')
f.write('Circle(27) = {21, 19, 17};\n')

f.write('Line(28) = {23, 24};\n')
f.write('Line(29) = {24, 26};\n')
f.write('Line(30) = {30, 29};\n')
f.write('Line(31) = {25, 23};\n')
f.write('Circle(32) = {26, 28, 30};\n')
f.write('Circle(33) = {29, 27, 25};\n')

f.write('Line(34) = {16, 24};\n')
f.write('Line(35) = {18, 26};\n')
f.write('Line(36) = {22, 30};\n')
f.write('Line(37) = {21, 29};\n')
f.write('Line(38) = {17, 25};\n')
f.write('Line(39) = {15, 23};\n')




f.write('// control volume\n')
f.write('Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, h/2.5};\n')
f.write('Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, h/2.5};\n')
f.write('Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, h/2.5};\n')
f.write('Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, h/2.5};\n')
f.write('Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, h/2.5};\n')
f.write('Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, h/2.5};\n')
f.write('Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, h/2.5};\n')
f.write('Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, h/2.5};\n')
f.write('Line(75) = {33, 37};\n')
f.write('Line(76) = {37, 35};\n')
f.write('Line(77) = {35, 31};\n')
f.write('Line(78) = {31, 33};\n')
f.write('Line(79) = {33, 34};\n')
f.write('Line(80) = {37, 38};\n')
f.write('Line(81) = {38, 36};\n')
f.write('Line(82) = {36, 32};\n')
f.write('Line(83) = {32, 34};\n')
f.write('Line(84) = {34, 38};\n')
f.write('Line(85) = {35, 36};\n')
f.write('Line(86) = {31, 32};\n')
f.write('Line Loop(87) = {75, 80, -84, -79};\n')
f.write('Plane Surface(88) = {87};\n')
f.write('Line Loop(89) = {76, 85, -81, -80};\n')
f.write('Plane Surface(90) = {89};\n')
f.write('Line Loop(91) = {77, 86, -82, -85};\n')
f.write('Plane Surface(92) = {91};\n')
f.write('Line Loop(93) = {78, 79, -83, -86};\n')
f.write('Plane Surface(94) = {93};\n')
f.write('Line Loop(95) = {75, 76, 77, 78};\n')
f.write('Plane Surface(96) = {95};\n')
f.write('Line Loop(97) = {84, 81, 82, 83};\n')
f.write('Plane Surface(98) = {97};\n')

f.write('// volume\n')


f.write('Line Loop(99) = {25, 39, -31, -38};\n')
f.write('Plane Surface(100) = {99};\n')
f.write('Line Loop(101) = {25, 22, 23, 26, 24, 27};\n')
f.write('Plane Surface(102) = {101};\n')
f.write('Line Loop(103) = {24, 37, -30, -36};\n')
f.write('Plane Surface(104) = {103};\n')
f.write('Line Loop(105) = {23, 35, -29, -34};\n')
f.write('Plane Surface(106) = {105};\n')
f.write('Line Loop(107) = {31, 28, 29, 32, 30, 33};\n')
f.write('Plane Surface(108) = {107};\n')
f.write('Line Loop(109) = {38, -33, -37, 27};\n')

f.write('Ruled Surface(110) = {109};\n')
f.write('Line Loop(111) = {26, 36, -32, -35};\n')
f.write('Ruled Surface(112) = {111};\n')
f.write('Line Loop(113) = {6, 14, -20, -13};\n')
f.write('Ruled Surface(114) = {113};\n')
f.write('Line Loop(115) = {5, 13, -19, -12};\n')
f.write('Ruled Surface(116) = {115};\n')
f.write('Line Loop(117) = {4, 12, -18, -11};\n')
f.write('Ruled Surface(118) = {117};\n')

f.write('Line Loop(119) = {3, 200, 11, -17, -10};\n')
f.write('Plane Surface(120) = {119};\n')

f.write('Line Loop(121) = {2, 10, -16, -9};\n')
f.write('Ruled Surface(122) = {121};\n')
f.write('Line Loop(123) = {1, 9, -15, -8};\n')
f.write('Ruled Surface(124) = {123};\n')
f.write('Line Loop(125) = {7, 8, -21, -14};\n')
f.write('Ruled Surface(126) = {125};\n')

f.write('Line Loop(127) = {7, 1, 2, 3, 200, 4, 5, 6};\n')
f.write('Line Loop(128) = {28, -34, -22, 39};\n')
f.write('Plane Surface(129) = {127, 128};\n')

f.write('Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};\n')
f.write('Plane Surface(131) = {130};\n')



f.write('Translate {0, 0, ee} {\n')
f.write('  Duplicata { Surface{131}; }\n')
f.write('}\n')




f.write('Line(140) = {6, 101};\n')
f.write('Line(141) = {7, 122};\n')
f.write('Line(142) = {14, 118};\n')
f.write('Line(143) = {13, 114};\n')
f.write('Line(144) = {12, 110};\n')
f.write('Line(145) = {5, 102};\n')
f.write('Line(146) = {8, 106};\n')




f.write('Line Loop(147) = {21, 145, -202, -140};\n')
f.write('Plane Surface(148) = {147};\n')
f.write('Line Loop(149) = {20, 140, -208, -141};\n')
f.write('Plane Surface(150) = {149};\n')
f.write('Line Loop(151) = {19, 141, -207, -142};\n')
f.write('Plane Surface(152) = {151};\n')
f.write('Line Loop(153) = {18, 142, -206, -143};\n')
f.write('Plane Surface(154) = {153};\n')
f.write('Line Loop(155) = {17, 143, -205, -144};\n')
f.write('Plane Surface(156) = {155};\n')
f.write('Line Loop(157) = {16, 144, -204, -146};\n')
f.write('Plane Surface(158) = {157};\n')
f.write('Line Loop(159) = {15, 146, -203, -145};\n')
f.write('Plane Surface(160) = {159};\n')

f.write('Surface Loop(161) = {129, 126, 124, 122, 120, 118, 116, 114, 100, 102, 106, 112, 104, 110, 108, 131};\n')
f.write('Surface Loop(162) = {94, 96, 88, 90, 92, 98};\n')
f.write('Volume(1) = {161, 162};\n')
f.write('Volume(2) = {162};\n')

f.write('Physical Volume(1) = {1,2}; // air cavity + controlled air volume\n')
f.write('Physical Volume(5) = {2}; // controlled air volume (small)\n')

f.write('Surface Loop(163) = {201, 148, 160, 158, 156, 154, 152, 150, 131};\n')
f.write('Volume(3) = {163};\n')

f.write('Physical Volume(2) = {3}; // porous volume\n')

f.write('Physical Surface(3) = {131}; // [porous]-[air cavity] : interface surface\n')

f.write('Physical Surface(4) = {148, 150, 152, 156, 158, 160, 132, 154};// porous external sur\n')

f.close()
import os
os.system('gmsh -3 -optimize %s.geo' % filename)
