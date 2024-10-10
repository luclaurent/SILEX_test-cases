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



filename='geom/cavity8_with_porous_air'

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

f.write('lpy1=0.8;\n')
f.write('lpz1=0.6;\n')
f.write('lpy2=2.4;\n')
f.write('lpz2=1.3;\n')

f.write('// size of elements\n')
f.write('h = lx1/40;\n')
f.write('h2 = ee/2;\n')


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

f.write('Point(5) = {0,     0  , lz1, h};\n')
f.write('Point(6) = {lx1,    0  , lz1, h};\n')
f.write('Point(7) = {lx1,    ly1 , lz1, h};\n')
f.write('Point(8) = {0 ,    ly1 , lz1, h};\n')

f.write('Point(9) = {0,     ly1+ly2  , 0, h};\n')
f.write('Point(10) = {lx2,     ly1+ly2  , 0, h};\n')
f.write('Point(11) = {lx2,     ly1  , 0, h};\n')
f.write('Point(12) = {0,     ly1+ly2  , lz1, h};\n')
f.write('Point(13) = {lx2,     ly1+ly2  , lz1, h};\n')
f.write('Point(14) = {lx2,     ly1  , lz1, h};\n')

f.write('// Cavity: lines\n')
f.write('Line(1) = {1, 4};\n')
f.write('Line(2) = {4, 9};\n')
f.write('Line(3) = {9, 10};\n')
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

f.write('// wall: corners\n')
f.write('Point(15) = {lx3,     ly3  , 0, h/3};\n')
f.write('Point(16) = {lx3+l4,     ly3  , 0, h/3};\n')
f.write('Point(17) = {lx3,     ly3  , lz4-r4 , h/3};\n')
f.write('Point(18) = {lx3+l4,     ly3  , lz4-r4, h/3};\n')
f.write('Point(19) = {lx3+r4,     ly3  , lz4-r4 , h/3};\n')
f.write('Point(20) = {lx3+l4-r4,     ly3  , lz4-r4, h/3};\n')
f.write('Point(21) = {lx3+r4,     ly3  , lz4 , (h+h)/2};\n')
f.write('Point(22) = {lx3+l4-r4,     ly3  , lz4, (h+h)/2};\n')

f.write('// rotate wall\n')
f.write('Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {\n')
f.write('  Point{ 18, 20, 21, 19, 17, 16, 15, 22};\n')
f.write('}\n')

f.write('Line(140) = {15, 16};\n')
f.write('Line(141) = {16, 18};\n')
f.write('Line(142) = {22, 21};\n')
f.write('Line(143) = {17, 15};\n')
f.write('Circle(144) = {18, 20, 22};\n')
f.write('Circle(145) = {21, 19, 17};\n')
f.write('Line Loop(154) = {141, 144, 142, 145, 143, 140};\n')
f.write('Plane Surface(155) = {154};\n')

f.write('// volume\n')




f.write('// Line Loop(113) = {6, 14, -20, -13};\n')
f.write('// Ruled Surface(114) = {113};\n')
f.write('Line Loop(115) = {5, 13, -19, -12};\n')
f.write('Ruled Surface(116) = {115};\n')
f.write('Line Loop(117) = {4, 12, -18, -11};\n')
f.write('Ruled Surface(118) = {117};\n')
f.write('Line Loop(119) = {3, 11, -17, -10};\n')
f.write('Ruled Surface(120) = {119};\n')
f.write('Line Loop(121) = {2, 10, -16, -9};\n')
f.write('Ruled Surface(122) = {121};\n')
f.write('Line Loop(123) = {1, 9, -15, -8};\n')
f.write('Ruled Surface(124) = {123};\n')
f.write('Line Loop(125) = {7, 8, -21, -14};\n')
f.write('Ruled Surface(126) = {125};\n')
f.write('Line Loop(127) = {7, 1, 2, 3, 4, 5, 6};\n')
f.write('Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};\n')
f.write('Plane Surface(131) = {130};\n')
f.write('Plane Surface(156) = {127};\n')
f.write('Surface Loop(157) = {126, 156, 124, 122, 120, 118, 116, 114, 131};\n')

f.write('// Porous media\n')
f.write('Point(200) = {lx1,	lpy1		, lpz1,		h2};\n')
f.write('Point(300) = {lx1,	lpy1+lpy2	, lpz1,		h2};\n')
f.write('Point(600) = {lx1,	lpy1		, lpz1+lpz2,	h2};\n')
f.write('Point(700) = {lx1,	lpy1+lpy2  	, lpz1+lpz2,	h2};\n')
f.write('Point(201) = {lx1+ee,	lpy1		, lpz1,		h2};\n')
f.write('Point(301) = {lx1+ee,	lpy1+lpy2	, lpz1,		h2};\n')
f.write('Point(601) = {lx1+ee,	lpy1		, lpz1+lpz2,	h2};\n')
f.write('Point(701) = {lx1+ee,	lpy1+lpy2	, lpz1+lpz2,	h2};\n')

f.write('Line(158) = {200, 300};\n')
f.write('Line(159) = {300, 700};\n')
f.write('Line(160) = {700, 600};\n')
f.write('Line(161) = {600, 200};\n')
f.write('Line Loop(162) = {6, 14, -20, -13};\n')
f.write('Line Loop(163) = {161, 158, 159, 160};\n')
f.write('Plane Surface(164) = {162, 163};\n')
f.write('Plane Surface(165) = {163};\n')

f.write('Line(166) = {600, 601};\n')
f.write('Line(167) = {700, 701};\n')
f.write('Line(168) = {300, 301};\n')
f.write('Line(169) = {200, 201};\n')
f.write('Line(170) = {201, 301};\n')
f.write('Line(171) = {301, 701};\n')
f.write('Line(172) = {701, 601};\n')
f.write('Line(173) = {601, 201};\n')
f.write('Line Loop(174) = {170, -168, -158, 169};\n')
f.write('Plane Surface(175) = {174};\n')
f.write('Line Loop(176) = {168, 171, -167, -159};\n')
f.write('Plane Surface(177) = {176};\n')
f.write('Line Loop(178) = {160, 166, -172, -167};\n')
f.write('Plane Surface(179) = {178};\n')
f.write('Line Loop(180) = {161, 169, -173, -166};\n')
f.write('Plane Surface(181) = {180};\n')
f.write('Line Loop(182) = {173, 170, 171, 172};\n')
f.write('Plane Surface(183) = {182};\n')

f.write('// controlled air volume\n')
f.write('Surface Loop(184) = {94, 96, 88, 90, 92, 98};\n')
f.write('Volume(185) = {184};\n')

f.write('//Physical Volume(1) = {278,276}; // air cavity with no controlled air volume\n')
f.write('Surface Loop(186) = {126, 156, 124, 122, 120, 118, 116, 164, 131, 165};\n')
f.write('Volume(187) = {184, 186};\n')

f.write('// porous volume\n')
f.write('Surface Loop(188) = {165, 179, 181, 175, 183, 177};\n')
f.write('Volume(189) = {188};\n')


f.write('Physical Volume(1) = {185,187}; // air cavity + controlled air volume\n')

f.write('Physical Volume(5) = {185}; // controlled air volume (small)\n')

f.write('Physical Volume(2) = {189}; // porous volume\n')

f.write('Physical Surface(3) = {165}; // [porous]-[air cavity] : interface surface\n')

f.write('Physical Surface(4) = {181, 179, 177, 175, 183};// porous external surface\n')

f.write('Physical Surface(6) = {155}; // Structure surface\n')

f.write('Physical Line(7) = {143, 145, 142, 144, 141};// Structure edge of the free boundary in air\n')



f.close()
import os
os.system('gmsh -3 -optimize %s.geo' % filename)
########################################
########################################""
########################################""
########################################""
########################################""

filename='geom/cavity8_with_porous_air_struc'

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
f.write('h = lx1/60;\n')
f.write('h2 = ee/3;\n')


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

f.write('// Cavity: lines\n')
f.write('Line(1) = {1, 4};\n')
f.write('Line(2) = {4, 9};\n')
f.write('Line(3) = {9, 10};\n')
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
f.write('//Line Loop(87) = {75, 80, -84, -79};\n')
f.write('//Plane Surface(88) = {87};\n')
f.write('//Line Loop(89) = {76, 85, -81, -80};\n')
f.write('//Plane Surface(90) = {89};\n')
f.write('//Line Loop(91) = {77, 86, -82, -85};\n')
f.write('//Plane Surface(92) = {91};\n')
f.write('//Line Loop(93) = {78, 79, -83, -86};\n')
f.write('//Plane Surface(94) = {93};\n')
f.write('//Line Loop(95) = {75, 76, 77, 78};\n')
f.write('//Plane Surface(96) = {95};\n')
f.write('//Line Loop(97) = {84, 81, 82, 83};\n')
f.write('//Plane Surface(98) = {97};\n')

f.write('// wall: corners\n')
f.write('Point(15) = {lx3,     ly3  , 0, h/3};\n')
f.write('Point(16) = {lx3+l4,     ly3  , 0, h/3};\n')
f.write('Point(17) = {lx3,     ly3  , lz4-r4 , h/3};\n')
f.write('Point(18) = {lx3+l4,     ly3  , lz4-r4, h/3};\n')
f.write('Point(19) = {lx3+r4,     ly3  , lz4-r4 , h/3};\n')
f.write('Point(20) = {lx3+l4-r4,     ly3  , lz4-r4, h/3};\n')
f.write('Point(21) = {lx3+r4,     ly3  , lz4 , (h2+h)/2};\n')
f.write('Point(22) = {lx3+l4-r4,     ly3  , lz4, (h2+h)/2};\n')

f.write('// rotate wall\n')
f.write('Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {\n')
f.write('  Point{ 18, 20, 21, 19, 17, 16, 15, 22};\n')
f.write('}\n')

f.write('Line(140) = {15, 16};\n')
f.write('Line(141) = {16, 18};\n')
f.write('Line(142) = {22, 21};\n')
f.write('Line(143) = {17, 15};\n')
f.write('Circle(144) = {18, 20, 22};\n')
f.write('Circle(145) = {21, 19, 17};\n')
f.write('Line Loop(154) = {141, 144, 142, 145, 143, 140};\n')
f.write('Plane Surface(155) = {154};\n')

f.write('// volume\n')




f.write('//Line Loop(113) = {6, 14, -20, -13};\n')
f.write('//Ruled Surface(114) = {113};\n')
f.write('//Line Loop(115) = {5, 13, -19, -12};\n')
f.write('//Ruled Surface(116) = {115};\n')
f.write('//Line Loop(117) = {4, 12, -18, -11};\n')
f.write('//Ruled Surface(118) = {117};\n')
f.write('//Line Loop(119) = {3, 11, -17, -10};\n')
f.write('//Ruled Surface(120) = {119};\n')
f.write('//Line Loop(121) = {2, 10, -16, -9};\n')
f.write('//Ruled Surface(122) = {121};\n')
f.write('//Line Loop(123) = {1, 9, -15, -8};\n')
f.write('//Ruled Surface(124) = {123};\n')
f.write('//Line Loop(125) = {7, 8, -21, -14};\n')
f.write('//Ruled Surface(126) = {125};\n')
f.write('//Line Loop(127) = {7, 1, 2, 3, 4, 5, 6};\n')
f.write('//Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};\n')
f.write('//Plane Surface(131) = {130};\n')

f.write('// Porous media\n')
f.write('//Translate {0, 0, ee} {\n')
f.write('//  Duplicata { Surface{131}; }\n')
f.write('//}\n')
f.write('//Line(256) = {12, 48};\n')
f.write('//Line(257) = {5, 40};\n')
f.write('//Line(258) = {6, 39};\n')
f.write('//Line(259) = {7, 60};\n')
f.write('//Line(260) = {14, 56};\n')
f.write('//Line(261) = {13, 52};\n')
f.write('//Line Loop(262) = {21, 257, -157, -258};\n')
f.write('//Plane Surface(263) = {262};\n')
f.write('//Line Loop(264) = {258, -163, -259, 20};\n')
f.write('//Plane Surface(265) = {264};\n')
f.write('//Line Loop(266) = {259, -162, -260, 19};\n')
f.write('//Plane Surface(267) = {266};\n')
f.write('//Line Loop(268) = {260, -161, -261, 18};\n')
f.write('//Plane Surface(269) = {268};\n')
f.write('//Line Loop(270) = {17, 261, -160, -256};\n')
f.write('//Plane Surface(271) = {270};\n')
f.write('//Line Loop(272) = {15, 16, 256, -159, -158, -257};\n')
f.write('//Plane Surface(273) = {272};\n')

f.write('//Plane Surface(274) = {127};\n')


f.write('//Surface Loop(275) = {96, 88, 90, 92, 94, 98};\n')
f.write('//Volume(276) = {275};\n')

f.write('//Surface Loop(277) = {124, 274, 126, 114, 116, 118, 120, 122, 131};\n')
f.write('//Volume(278) = {275, 277};\n')

f.write('//Surface Loop(279) = {271, 269, 267, 265, 263, 273, 156, 131};\n')
f.write('//Volume(280) = {279};\n')

f.write('//Physical Volume(1) = {278,276}; // air cavity + controlled air volume\n')

f.write('//Physical Volume(5) = {276}; // controlled air volume (small)\n')

f.write('//Physical Volume(2) = {280}; // porous volume\n')

f.write('//Physical Surface(3) = {131}; // [porous]-[air cavity] : interface surface\n')

f.write('//Physical Surface(4) = {265, 267, 269, 271, 273, 263, 156};// porous external surface\n')

f.write('Physical Surface(6) = {155}; // Structure surface\n')

f.write('Physical Line(7) = {143, 145, 142, 144, 141};// Structure edge of the free boundary in air\n')



f.close()
import os
os.system('gmsh -3 -optimize %s.geo' % filename)
