// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.5;

ly2 = 1.0;
lx2 = 3.5;
lx3 =         2;

ly4 = 1.2; // structure thickness

ly3 = 2.0;
l4 = 2.5;
lz4 = 1.2;
r4   = 0.5;
e4   = 0.2;
deg=       30;
angle = deg*Pi/180;


// porous
ee = 0.2;

lpy1=0.8;
lpz1=0.6;
lpy2=2.4;
lpz2=1.3;

// size of elements
// air
h = lx1/45;
// air control volume
hc = h/1.5;
// structure
hs = h/2;
// porous
h2 = ee; // /2;


lx5 = 6.0;
ly5 = 1.0;
lz5 = 1.5;
h5  = 1.0;

Mesh.CharacteristicLengthMax=1.5*h;
Mesh.ElementOrder = 1;

// Cavity: Corners
Point(1) = {0,     0  , 0, h};
Point(2) = {lx1,    0  , 0, h};
Point(3) = {lx1,    ly1 , 0, h};
Point(4) = {0 ,    ly1 , 0, h};

Point(5) = {0,     0  , lz1, h};
Point(6) = {lx1,    0  , lz1, h};
Point(7) = {lx1,    ly1 , lz1, h};
Point(8) = {0 ,    ly1 , lz1, h};

Point(9) = {0,     ly1+ly2  , 0, h};
Point(10) = {lx2,     ly1+ly2  , 0, h};
Point(11) = {lx2,     ly1  , 0, h};
Point(12) = {0,     ly1+ly2  , lz1, h};
Point(13) = {lx2,     ly1+ly2  , lz1, h};
Point(14) = {lx2,     ly1  , lz1, h};

// Cavity: lines
Line(1) = {1, 4};
Line(2) = {4, 9};
Line(3) = {9, 10};
Line(4) = {10, 11};
Line(5) = {11, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 5};
Line(9) = {4, 8};
Line(10) = {9, 12};
Line(11) = {10, 13};
Line(12) = {11, 14};
Line(13) = {3, 7};
Line(14) = {2, 6};
Line(15) = {5, 8};
Line(16) = {8, 12};
Line(17) = {12, 13};
Line(18) = {13, 14};
Line(19) = {14, 7};
Line(20) = {7, 6};
Line(21) = {6, 5};

// control volume
Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, hc};
Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, hc};
Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, hc};
Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, hc};
Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, hc};
Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, hc};
Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, hc};
Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, hc};
Line(75) = {33, 37};
Line(76) = {37, 35};
Line(77) = {35, 31};
Line(78) = {31, 33};
Line(79) = {33, 34};
Line(80) = {37, 38};
Line(81) = {38, 36};
Line(82) = {36, 32};
Line(83) = {32, 34};
Line(84) = {34, 38};
Line(85) = {35, 36};
Line(86) = {31, 32};
Line Loop(87) = {75, 80, -84, -79};
Plane Surface(88) = {87};
Line Loop(89) = {76, 85, -81, -80};
Plane Surface(90) = {89};
Line Loop(91) = {77, 86, -82, -85};
Plane Surface(92) = {91};
Line Loop(93) = {78, 79, -83, -86};
Plane Surface(94) = {93};
Line Loop(95) = {75, 76, 77, 78};
Plane Surface(96) = {95};
Line Loop(97) = {84, 81, 82, 83};
Plane Surface(98) = {97};

// wall: corners
Point(15) = {lx3,     ly3-ly4/2  , 0, hs};
Point(16) = {lx3+l4,     ly3-ly4/2  , 0, hs};
Point(17) = {lx3,     ly3-ly4/2  , lz4-r4 , hs};
Point(18) = {lx3+l4,     ly3-ly4/2  , lz4-r4, hs};
Point(19) = {lx3+r4,     ly3-ly4/2  , lz4-r4 , hs};
Point(20) = {lx3+l4-r4,     ly3-ly4/2  , lz4-r4, hs};
Point(21) = {lx3+r4,     ly3-ly4/2  , lz4 , (h+h)/2};
Point(22) = {lx3+l4-r4,     ly3-ly4/2  , lz4, (h+h)/2};

Point(515) = {lx3,     ly3+ly4/2  , 0, hs};
Point(516) = {lx3+l4,     ly3+ly4/2  , 0, hs};
Point(517) = {lx3,     ly3+ly4/2  , lz4-r4 , hs};
Point(518) = {lx3+l4,     ly3+ly4/2  , lz4-r4, hs};
Point(519) = {lx3+r4,     ly3+ly4/2  , lz4-r4 , hs};
Point(520) = {lx3+l4-r4,     ly3+ly4/2  , lz4-r4, hs};
Point(521) = {lx3+r4,     ly3+ly4/2  , lz4 , (h+h)/2};
Point(522) = {lx3+l4-r4,     ly3+ly4/2  , lz4, (h+h)/2};


Point(615) = {lx3-ly4/2,     ly3  , 0, hs};
Point(616) = {lx3+l4+ly4/2,     ly3  , 0, hs};
Point(617) = {lx3-ly4/2,     ly3  , lz4-r4 , hs};
Point(618) = {lx3+l4+ly4/2,     ly3  , lz4-r4, hs};
Point(619) = {lx3+r4,     ly3  , lz4-r4 , hs};
Point(620) = {lx3+l4-r4,     ly3  , lz4-r4, hs};
Point(621) = {lx3+r4,     ly3  , lz4+ly4/2 , (h+h)/2};
Point(622) = {lx3+l4-r4,     ly3  , lz4+ly4/2, (h+h)/2};

Point(715) = {lx3,     ly3  , 0, hs};
Point(716) = {lx3+l4,     ly3  , 0, hs};
Point(717) = {lx3,     ly3  , lz4-r4 , hs};
Point(718) = {lx3+l4,     ly3  , lz4-r4, hs};
Point(721) = {lx3+r4,     ly3  , lz4 , (h+h)/2};
Point(722) = {lx3+l4-r4,     ly3  , lz4, (h+h)/2};

// rotate wall
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
  Point{ 18, 20, 21, 19, 17, 16, 15, 22};
}
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
  Point{ 518, 520, 521, 519, 517, 516, 515, 522};
}
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
  Point{ 618, 620, 621, 619, 617, 616, 615, 622,715,716,717,718,721,722};
}

Line(140) = {15, 16};
Line(141) = {16, 18};
Line(142) = {22, 21};
Line(143) = {17, 15};
Circle(144) = {18, 20, 22};
Circle(145) = {21, 19, 17};
Line Loop(154) = {141, 144, 142, 145, 143, 140};
Plane Surface(155) = {154};

Line(5140) = {515, 516};
Line(5141) = {516, 518};
Line(5142) = {522, 521};
Line(5143) = {517, 515};
Circle(5144) = {518, 520, 522};
Circle(5145) = {521, 519, 517};
Line Loop(5154) = {5141, 5144, 5142, 5145, 5143, 5140};
Plane Surface(5155) = {5154};

Circle(6000) = {515, 715, 615};
Circle(6001) = {615, 715, 15};

Circle(6002) = {517, 717, 617};
Circle(6003) = {617, 717, 17};

Circle(6004) = {516, 716, 616};
Circle(6005) = {616, 716, 16};

Circle(6006) = {518, 718, 618};
Circle(6007) = {618, 718, 18};

Circle(6008) = {522, 722, 622};
Circle(6009) = {622, 722, 22};

Circle(6010) = {521, 721, 621};
Circle(6011) = {621, 721, 21};

Line(7000) = {616, 618};
Line(7001) = {622, 621};
Line(7002) = {615, 617};

Circle(6012) = {618, 620, 622};
Circle(6013) = {621, 619, 617};

// structure volume

//+
Curve Loop(5155) = {5141, 6006, -7000, -6004};
//+
Surface(5156) = {5155};
//+
Curve Loop(5156) = {7000, 6007, -141, -6005};
//+
Surface(5157) = {5156};
//+
Curve Loop(5157) = {5143, 6000, 7002, -6002};
//+
Surface(5158) = {5157};
//+
Curve Loop(5158) = {7002, 6003, 143, -6001};
//+
Surface(5159) = {5158};
//+
Curve Loop(5159) = {5142, 6010, -7001, -6008};
//+
Surface(5160) = {5159};
//+
Curve Loop(5160) = {7001, 6011, -142, -6009};
//+
Surface(5161) = {5160};
//+
Curve Loop(5161) = {5144, 6008, -6012, -6006};
//+
Surface(5162) = {5161};
//+
Curve Loop(5162) = {6012, 6009, -144, -6007};
//+
Surface(5163) = {5162};
//+
Curve Loop(5163) = {6002, -6013, -6010, 5145};
//+
Surface(5164) = {5163};
//+
Curve Loop(5164) = {6003, -145, -6011, 6013};
//+
Surface(5165) = {5164};
//+
Curve Loop(5165) = {5140, 6004, 6005, -140, -6001, -6000};
//+
Plane Surface(5166) = {5165};

//+
Surface Loop(189) = {5155, 5156, 5162, 5160, 5164, 5158, 5166, 5157, 5163, 5161, 5165, 5159, 155};
//+ structure volume
Volume(190) = {189};


// Line Loop(113) = {6, 14, -20, -13};
// Ruled Surface(114) = {113};
Line Loop(115) = {5, 13, -19, -12};
Ruled Surface(116) = {115};
Line Loop(117) = {4, 12, -18, -11};
Ruled Surface(118) = {117};
Line Loop(119) = {3, 11, -17, -10};
Ruled Surface(120) = {119};
Line Loop(121) = {2, 10, -16, -9};
Ruled Surface(122) = {121};
Line Loop(123) = {1, 9, -15, -8};
Ruled Surface(124) = {123};
Line Loop(125) = {7, 8, -21, -14};
Ruled Surface(126) = {125};
Line Loop(127) = {7, 1, 2, 3, 4, 5, 6};
Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};
Plane Surface(131) = {130};
Plane Surface(156) = {127};
Surface Loop(157) = {126, 156, 124, 122, 120, 118, 116, 114, 131};

// Porous media
Point(200) = {lx1,	lpy1		, lpz1,		h2};
Point(300) = {lx1,	lpy1+lpy2	, lpz1,		h2};
Point(600) = {lx1,	lpy1		, lpz1+lpz2,	h2};
Point(700) = {lx1,	lpy1+lpy2  	, lpz1+lpz2,	h2};
Point(201) = {lx1+ee,	lpy1		, lpz1,		h2};
Point(301) = {lx1+ee,	lpy1+lpy2	, lpz1,		h2};
Point(601) = {lx1+ee,	lpy1		, lpz1+lpz2,	h2};
Point(701) = {lx1+ee,	lpy1+lpy2	, lpz1+lpz2,	h2};

Line(158) = {200, 300};
Line(159) = {300, 700};
Line(160) = {700, 600};
Line(161) = {600, 200};
Line Loop(162) = {6, 14, -20, -13};
Line Loop(163) = {161, 158, 159, 160};
Plane Surface(164) = {162, 163};
Plane Surface(165) = {163};

Line(166) = {600, 601};
Line(167) = {700, 701};
Line(168) = {300, 301};
Line(169) = {200, 201};
Line(170) = {201, 301};
Line(171) = {301, 701};
Line(172) = {701, 601};
Line(173) = {601, 201};
Line Loop(174) = {170, -168, -158, 169};
Plane Surface(175) = {174};
Line Loop(176) = {168, 171, -167, -159};
Plane Surface(177) = {176};
Line Loop(178) = {160, 166, -172, -167};
Plane Surface(179) = {178};
Line Loop(180) = {161, 169, -173, -166};
Plane Surface(181) = {180};
Line Loop(182) = {173, 170, 171, 172};
Plane Surface(183) = {182};

// controlled air volume
Surface Loop(184) = {94, 96, 88, 90, 92, 98};
Volume(185) = {184};

//Physical Volume(1) = {278,276}; // air cavity with no controlled air volume
Surface Loop(186) = {126, 156, 124, 122, 120, 118, 116, 164, 131, 165};
Volume(187) = {184, 186};

// porous volume
Surface Loop(188) = {165, 179, 181, 175, 183, 177};
Volume(189) = {188};


Physical Volume(1) = {185,187}; // air cavity + controlled air volume

Physical Volume(5) = {185}; // controlled air volume (small)

Physical Volume(2) = {189}; // porous volume

Physical Surface(3) = {165}; // [porous]-[air cavity] : interface surface

Physical Surface(4) = {181, 179, 177, 175, 183};// porous external surface

Physical Surface(6) = {155, 5155, 5159, 5158, 5164, 5165, 5161, 5160, 5163, 5162, 5157, 5156}; // Structure surface

