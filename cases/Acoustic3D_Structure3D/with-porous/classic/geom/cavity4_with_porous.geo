// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.5;
ly2 = 3.0;
lx2 = 3.5;
lx3 = 1.5;
ly3 = 2.0;
l4 = 3.0;
lz4 = 2.4;
r4   = 0.5;
e4   = 0.1;
deg=      0;
angle = deg*Pi/180;
// porous on the moving wall
ee = 0.1;
// porous on the fixed wall
ee6 = 0.2;
dy6 = 0.1;
dz6 = 0.1;
ly6 = 6.8;
lz6 = 2.3;

// size of elements
h = lx1/40;
h2 = ee/4;
h6 = ee6/4;

lx5 = 6.0;
ly5 = 1.0;
lz5 = 1.5;
h5  = 1.0;
Mesh.CharacteristicLengthMax=2*h;
Mesh.ElementOrder = 1;
// Cavity: Corners
Point(1) = {0,     0  , 0, h};
Point(2) = {lx1,    0  , 0, h};
Point(3) = {lx1,    ly1 , 0, h};
Point(5) = {0,     0  , lz1, h};
Point(6) = {lx1,    0  , lz1, h};
Point(7) = {lx1,    ly1 , lz1, h};
Point(9) = {0,     ly1+ly2  , 0, h};
Point(10) = {lx2,     ly1+ly2  , 0, h};
Point(11) = {lx2,     ly1  , 0, h};
Point(12) = {0,     ly1+ly2  , lz1, h};
Point(13) = {lx2,     ly1+ly2  , lz1, h};
Point(14) = {lx2,     ly1  , lz1, h};
// Cavity: lines
Line(2) = {1, 9};
Line(3) = {9, 10};
Line(4) = {10, 11};
Line(5) = {11, 3};
Line(6) = {3, 2};
Line(7) = {2, 1};
Line(8) = {1, 5};
Line(10) = {9, 12};
Line(11) = {10, 13};
Line(12) = {11, 14};
Line(13) = {3, 7};
Line(14) = {2, 6};
Line(16) = {5, 12};
Line(17) = {12, 13};
Line(18) = {13, 14};
Line(19) = {14, 7};
Line(20) = {7, 6};
Line(21) = {6, 5};
// wall: corners
Point(15) = {lx3,     ly3-e4/2  , 0, h};
Point(16) = {lx3+l4,     ly3-e4/2  , 0, h};
Point(17) = {lx3,     ly3-e4/2  , lz4-r4 , h};
Point(18) = {lx3+l4,     ly3-e4/2  , lz4-r4, h};
Point(19) = {lx3+r4,     ly3-e4/2  , lz4-r4 , h};
Point(20) = {lx3+l4-r4,     ly3-e4/2  , lz4-r4, h};
Point(21) = {lx3+r4,     ly3-e4/2  , lz4 , h};
Point(22) = {lx3+l4-r4,     ly3-e4/2  , lz4, h};
Point(23) = {lx3,     ly3+e4/2  , 0, h2};
Point(24) = {lx3+l4,     ly3+e4/2  , 0, h2};
Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h2};
Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h2};
Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h2};
Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h2};
Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , h2};
Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, h2};
Point(123) = {lx3,     ly3+e4/2+ee  , 0, h2};
Point(124) = {lx3+l4,     ly3+e4/2+ee  , 0, h2};
Point(125) = {lx3,     ly3+e4/2+ee  , lz4-r4 , h2};
Point(126) = {lx3+l4,     ly3+e4/2+ee  , lz4-r4, h2};
Point(127) = {lx3+r4,     ly3+e4/2+ee  , lz4-r4 , h2};
Point(128) = {lx3+l4-r4,     ly3+e4/2+ee  , lz4-r4, h2};
Point(129) = {lx3+r4,     ly3+e4/2+ee  , lz4 , h2};
Point(130) = {lx3+l4-r4,     ly3+e4/2+ee  , lz4, h2};
// rotate wall
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
Point{22, 30, 18, 26, 20, 21, 28, 29, 19, 27, 17, 25, 16, 24, 15, 23,123,124,125,126,127,128,129,130};
}
// wall : lines
Line(22) = {15, 16};
Line(23) = {16, 18};
Line(24) = {22, 21};
Line(25) = {17, 15};
Circle(26) = {18, 20, 22};
Circle(27) = {21, 19, 17};
Line(28) = {23, 24};
Line(29) = {24, 26};
Line(30) = {30, 29};
Line(31) = {25, 23};
Circle(32) = {26, 28, 30};
Circle(33) = {29, 27, 25};
Line(128) = {123, 124};
Line(129) = {124, 126};
Line(130) = {130, 129};
Line(131) = {125, 123};
Circle(132) = {126, 128, 130};
Circle(133) = {129, 127, 125};
Line(34) = {16, 24};
Line(35) = {18, 26};
Line(36) = {22, 30};
Line(37) = {21, 29};
Line(38) = {17, 25};
Line(39) = {15, 23};
Line(134) = {123, 23};
Line(135) = {125, 25};
Line(136) = {129, 29};
Line(137) = {130, 30};
Line(138) = {126, 26};
Line(139) = {124, 24};

// Porous, surfaces and volume, on the moving wall
Line Loop(140) = {32, 30, 33, 31, 28, 29};
Plane Surface(141) = {140};
Line Loop(142) = {138, 32, -137, -132};
Ruled Surface(143) = {142};
Line Loop(144) = {130, 136, -30, -137};
Plane Surface(145) = {144};
Line Loop(146) = {133, 135, -33, -136};
Ruled Surface(147) = {146};
Line Loop(148) = {134, -31, -135, 131};
Plane Surface(149) = {148};
Line Loop(150) = {134, 28, -139, -128};
Plane Surface(151) = {150};
Line Loop(152) = {139, 29, -138, -129};
Plane Surface(153) = {152};
Line Loop(154) = {128, 129, 132, 130, 133, 131};
Plane Surface(155) = {154};
Surface Loop(156) = {155, 151, 149, 141, 143, 153, 145, 147};
Volume(157) = {156};

// porous on the fixed wall
Point(300) = {-ee6,  dy6  , dz6, h6};
Point(301) = {-ee6,  dy6  , dz6+lz6, h6};
Point(302) = {-ee6,  dy6+ly6  , dz6+lz6, h6};
Point(303) = {-ee6,  dy6+ly6  , dz6, h6};
Point(304) = {0,  dy6  , dz6, h6};
Point(305) = {0,  dy6  , dz6+lz6, h6};
Point(306) = {0,  dy6+ly6  , dz6+lz6, h6};
Point(307) = {0,  dy6+ly6  , dz6, h6};

Line(188) = {304, 305};
Line(189) = {305, 306};
Line(190) = {306, 307};
Line(191) = {307, 304};
Line(192) = {300, 301};
Line(193) = {301, 302};
Line(194) = {302, 303};
Line(195) = {303, 300};
Line(196) = {300, 304};
Line(197) = {301, 305};
Line(198) = {302, 306};
Line(199) = {303, 307};
Line Loop(200) = {188, 189, 190, 191};
Plane Surface(201) = {200};
Line Loop(202) = {197, 189, -198, -193};
Plane Surface(203) = {202};
Line Loop(204) = {196, 188, -197, -192};
Plane Surface(205) = {204};
Line Loop(206) = {195, 196, -191, -199};
Plane Surface(207) = {206};
Line Loop(208) = {199, -190, -198, 194};
Plane Surface(209) = {208};
Line Loop(210) = {195, 192, 193, 194};
Plane Surface(211) = {210};

Surface Loop(212) = {207, 211, 205, 203, 209, 201};
Volume(213) = {212};


// surfaces, cavity
Line Loop(170) = {6, 7, 2, 3, 4, 5};
Line Loop(171) = {22, 34, -139, -128, 134, -39};
Plane Surface(172) = {170, 171};
Line Loop(173) = {6, 14, -20, -13};
Plane Surface(174) = {173};
Line Loop(175) = {5, 13, -19, -12};
Plane Surface(176) = {175};
Line Loop(177) = {4, 12, -18, -11};
Plane Surface(178) = {177};
Line Loop(179) = {3, 11, -17, -10};
Plane Surface(180) = {179};
//Line Loop(181) = {2, 10, -16, -8};
//Plane Surface(182) = {181};
//Line Loop(181) = {8, 16, -10, -2};
//Plane Surface(182) = {200, 211};

Line Loop(181) = {16, -10, -2, 8};
Plane Surface(182) = {200, 181};

Line Loop(183) = {7, 8, -21, -14};
Plane Surface(184) = {183};
Line Loop(185) = {21, 16, 17, 18, 19, 20};
Plane Surface(186) = {185};
// surface, internal wall
Line Loop(158) = {39, -31, -38, 25};
Plane Surface(159) = {158};
Line Loop(160) = {38, -33, -37, 27};
Ruled Surface(161) = {160};
Line Loop(162) = {37, -30, -36, 24};
Plane Surface(163) = {162};
Line Loop(164) = {36, -32, -35, 26};
Ruled Surface(165) = {164};
Line Loop(166) = {35, -29, -34, 23};
Plane Surface(167) = {166};
Line Loop(168) = {22, 23, 26, 24, 27, 25};
Plane Surface(169) = {168};

// control volume
Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, h};
Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, h};
Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, h};
Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, h};
Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, h};
Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, h};
Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, h};
Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, h};
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

// controlled volume
Surface Loop(99) = {98, 88, 96, 90, 92, 94};
Volume(2) = {99};

// volume: air cavity + controlled volume
Surface Loop(187) = {155, 153, 143, 145, 147, 149, 163, 161, 159, 172, 174, 184, 182, 180, 178, 176, 186, 167, 165, 169,201};
Volume(1) = {99, 187};
//Surface Loop(101) = {58, 54, 52, 50, 48, 46, 44, 42, 56, 64, 72, 60, 66, 62, 68, 70};
//Volume(1) = {99, 101};

Physical Volume(1) = {1,2}; // air cavity + controlled air volume
Physical Volume(5) = {2}; // controlled air volume (small)
// Porous
Physical Volume(2) = {157,213}; // porous volume
Physical Surface(4) = {153, 143, 145, 147, 149, 151, 141,211,203,209,207,205}; // porous external surface
Physical Surface(3) = {155,201}; // [porous]-[air cavity] : interface surface
