// Parameters: acoustic cavity
lx = 5.0;
ly = 2.0;
lz = 2.0;
 
bx1  = 2.0;
lys = 0.7;
h1  = 0.6;
a1 = -0.1;
r  = 0.1;
a2 = 0.6;
h2 = 0.8;
h3 = 0.2;
a3 = 0.1;

bx2 = 4.0;
by2 = ly-lys;

//nf = lx/40;
nf = lx/12;
//ns = lys/25;
ns = lys/12;

Mesh.CharacteristicLengthMax=nf;

// Corners
Point(1) = {0,     0  , 0, nf};
Point(2) = {lx,    0  , 0, nf};
Point(3) = {lx,    ly , 0, nf};
Point(4) = {0 ,    ly , 0, nf};

Point(5) = {0,     0  , lz, nf};
Point(6) = {lx,    0  , lz, nf};
Point(7) = {lx,    ly , lz, nf};
Point(8) = {0 ,    ly , lz, nf};


// Structure 1
Point(9) = {bx1,     0    , 0, ns};
Point(10) = {bx1,     lys  , 0, ns};
Point(11) = {bx1-a1,     0    , h1, ns};
Point(12) = {bx1-a1,     lys  , h1, ns};
Point(13) = {bx1-a1-r,     0    , h1, ns};
Point(14) = {bx1-a1-r,     lys  , h1, ns};
Point(15) = {bx1-a1-r,     0    , h1+r, ns};
Point(16) = {bx1-a1-r,     lys  , h1+r, ns};
Point(17) = {bx1-a1-r-a2,     0    , h1+r, ns};
Point(18) = {bx1-a1-r-a2,     lys  , h1+r, ns};
Point(19) = {bx1-a1-r-a2,     0    , h1+r+r, ns};
Point(20) = {bx1-a1-r-a2,     lys  , h1+r+r, ns};
Point(21) = {bx1-a1-r-a2-r,     0    , h1+r+r, ns};
Point(22) = {bx1-a1-r-a2-r,     lys  , h1+r+r, ns};
Point(23) = {bx1-a1-r-a2-r-a3,     0    , h1+r+r+h2, ns};
Point(24) = {bx1-a1-r-a2-r-a3,     lys  , h1+r+r+h2, ns};
Point(25) = {bx1-a1-r-a2-r-a3,     lys/2  , h1+r+r+h2-h3, ns};
Point(26) = {bx1-a1-r-a2/2,     lys/2    , h1+r, ns};


Circle(1) = {23, 25, 24};

Line(2) = {9, 10};
Line(3) = {10, 12};
Line(4) = {12, 11};
Line(5) = {11, 9};

Line(6) = {15, 17};
Line(7) = {17, 18};
Line(8) = {18, 16};
Line(9) = {16, 15};
Line(10) = {22, 21};
Line(11) = {24, 22};
Line(12) = {23, 24};
Line(13) = {21, 23};
Circle(14) = {11, 13, 15};
Circle(15) = {16, 14, 12};
Circle(16) = {21, 19, 17};
Circle(17) = {18, 20, 22};

Line Loop(38) = {2, 3, 4, 5};
Plane Surface(39) = {38};
//Line Loop(40) = {-7, -8, -9, -6};
//Plane Surface(41) = {40};
Line Loop(42) = {-10, -13, -12, -11};
Plane Surface(43) = {42};
Line Loop(44) = {-1, 12};
Plane Surface(45) = {44};
Line Loop(46) = {-4, -14, 9, -15};
Ruled Surface(47) = {46};
Line Loop(48) = {7, 17, 10, 16};
Ruled Surface(49) = {48};

Line(202) = {17, 26};
Line(203) = {18, 26};
Line(204) = {16, 26};
Line(205) = {15, 26};


Line Loop(206) = {-203, 202, -7};
Plane Surface(207) = {206};
Line Loop(208) = {-6, -202, 205};
Plane Surface(209) = {208};
Line Loop(210) = {-9, -205, 204};
Plane Surface(211) = {210};
Line Loop(212) = {-204, 203, -8};
Plane Surface(213) = {212};


n=20;
// Structure 2
Point(9+n) = {bx2,     by2    , 0, ns};
Point(10+n) = {bx2,     lys+by2  , 0, ns};
Point(11+n) = {bx2-a1,     by2    , h1, ns};
Point(12+n) = {bx2-a1,     lys+by2  , h1, ns};
Point(13+n) = {bx2-a1-r,     by2   , h1, ns};
Point(14+n) = {bx2-a1-r,     lys+by2  , h1, ns};
Point(15+n) = {bx2-a1-r,     by2    , h1+r, ns};
Point(16+n) = {bx2-a1-r,     lys+by2  , h1+r, ns};
Point(17+n) = {bx2-a1-r-a2,     by2    , h1+r, ns};
Point(18+n) = {bx2-a1-r-a2,     lys+by2  , h1+r, ns};
Point(19+n) = {bx2-a1-r-a2,     by2    , h1+r+r, ns};
Point(20+n) = {bx2-a1-r-a2,     lys+by2  , h1+r+r, ns};
Point(21+n) = {bx2-a1-r-a2-r,     by2    , h1+r+r, ns};
Point(22+n) = {bx2-a1-r-a2-r,     lys+by2  , h1+r+r, ns};
Point(23+n) = {bx2-a1-r-a2-r-a3,     by2    , h1+r+r+h2, ns};
Point(24+n) = {bx2-a1-r-a2-r-a3,     lys+by2  , h1+r+r+h2, ns};
Point(25+n) = {bx2-a1-r-a2-r-a3,     lys/2+by2  , h1+r+r+h2-h3, ns};

Circle(1+n) = {23+n, 25+n, 24+n};

Line(2+n) = {9+n, 10+n};
Line(3+n) = {10+n, 12+n};
Line(4+n) = {12+n, 11+n};
Line(5+n) = {11+n, 9+n};

Line(6+n) = {15+n, 17+n};
Line(7+n) = {17+n, 18+n};
Line(8+n) = {18+n, 16+n};
Line(9+n) = {16+n, 15+n};
Line(10+n) = {22+n, 21+n};
Line(11+n) = {24+n, 22+n};
Line(12+n) = {23+n, 24+n};
Line(13+n) = {21+n, 23+n};
Circle(14+n) = {11+n, 13+n, 15+n};
Circle(15+n) = {16+n, 14+n, 12+n};
Circle(16+n) = {21+n, 19+n, 17+n};
Circle(17+n) = {18+n, 20+n, 22+n};

Line Loop(38+n) = {2+n, 3+n, 4+n, 5+n};
Plane Surface(39+n) = {38+n};
Line Loop(40+n) = {-7-n, -8-n, -9-n, -6-n};
Plane Surface(41+n) = {40+n};
Line Loop(42+n) = {-10-n, -13-n, -12-n, -11-n};
Plane Surface(43+n) = {42+n};
Line Loop(44+n) = {-1-n, 12+n};
Plane Surface(45+n) = {44+n};
Line Loop(46+n) = {-4-n, -14-n, 9+n, -15-n};
Ruled Surface(47+n) = {46+n};
Line Loop(48+n) = {7+n, 17+n, 10+n, 16+n};
Ruled Surface(49+n) = {48+n};





// fluid
Line(70) = {1, 9};
//Line(71) = {9, 2};
Line(72) = {2, 3};
Line(73) = {3, 30};
Line(75) = {4, 1};
Line(76) = {1, 5};
Line(77) = {4, 8};
Line(78) = {2, 6};
Line(79) = {3, 7};
Line(80) = {6, 7};

Point(100) = {bx1-a1-r-a2-r-a3,     0    , lz, nf};
Point(101) = {bx2-a1-r-a2-r-a3,     ly    , lz, nf};
Point(102) = {bx1-a1-r-a2-r-a3,     ly    , lz, nf};
Point(103) = {bx2-a1-r-a2-r-a3,     0    , lz, nf};
Point(104) = {bx1,     ly  , 0, nf};

Point(105) = {bx1-a1,     ly  , h1, nf};
Point(106) = {bx1-a1-r,     ly  , h1, nf};
Point(107) = {bx1-a1-r,     ly  , h1+r, nf};

Point(108) = {bx1-a1-r-a2,     ly  , h1+r, nf};
Point(109) = {bx1-a1-r-a2,     ly  , h1+r+r, nf};
Point(110) = {bx1-a1-r-a2-r,     ly  , h1+r+r, nf};

Point(111) = {bx1-a1-r-a2-r-a3,     ly  , h1+r+r+h2, nf};

Circle(81) = {107, 106, 105};
Circle(82) = {108, 109, 110};

Line(83) = {104, 105};
Line(84) = {107, 108};
Line(85) = {110, 111};
Line(86) = {111, 102};
Line(87) = {102, 100};
Line(88) = {100, 23};
Line(89) = {5, 100};
Line(90) = {8, 102};
Line(91) = {24, 111};
Line(92) = {22, 110};
Line(93) = {18, 108};
Line(94) = {16, 107};
Line(95) = {12, 105};
Line(96) = {10, 104};
Line(97) = {5, 8};
Line Loop(98) = {96, 83, -95, -3};
Plane Surface(99) = {98};
Line Loop(100) = {94, 84, -93, 8};
Plane Surface(101) = {100};
Line Loop(102) = {92, 85, -91, 11};
Plane Surface(103) = {102};
Line Loop(104) = {95, -81, -94, 15};
Ruled Surface(105) = {104};
Line Loop(106) = {93, 82, -92, -17};
Ruled Surface(107) = {106};
Line Loop(108) = {91, 86, 87, 88, 1};
Plane Surface(109) = {108};

//Point(200) = {bx1-a1-r-a2-r-a3,     0    , lz, nf};
//Point(201) = {bx2-a1-r-a2-r-a3,     ly    , lz, nf};
//Point(202) = {bx1-a1-r-a2-r-a3,     ly    , lz, nf};
//Point(203) = {bx2-a1-r-a2-r-a3,     0    , lz, nf};
Point(204) = {bx2,     0  , 0, nf};

Point(205) = {bx2-a1,     0  , h1, nf};
Point(206) = {bx2-a1-r,     0  , h1, nf};
Point(207) = {bx2-a1-r,     0  , h1+r, nf};

Point(208) = {bx2-a1-r-a2,     0  , h1+r, nf};
Point(209) = {bx2-a1-r-a2,     0  , h1+r+r, nf};
Point(210) = {bx2-a1-r-a2-r,     0  , h1+r+r, nf};

Point(211) = {bx2-a1-r-a2-r-a3,     0  , h1+r+r+h2, nf};

// curved surface on top
Point(300)= {0 ,    ly/2 , 0.8*lz, nf};
Point(301)= {lx ,    ly/2 , 0.8*lz, nf};
Point(302) = {bx1-a1-r-a2-r-a3,     ly/2    , 0.8*lz, nf};
Point(303) = {bx2-a1-r-a2-r-a3,     ly/2    , 0.8*lz, nf};

Circle(178) = {8, 300, 5};
Circle(179) = {102, 302, 100};
Circle(180) = {101, 303, 103};
Circle(181) = {7, 301, 6};



Line(110) = {2, 204};
Line(111) = {204, 29};
Line(112) = {204, 205};
Line(113) = {205, 31};
Line(114) = {207, 35};
Line(115) = {208, 207};
Line(116) = {208, 37};
Line(117) = {41, 210};
Line(118) = {210, 211};
Line(119) = {211, 103};
Line(120) = {103, 101};
Line(121) = {101, 44};
Line(130) = {211, 43};

Circle(122) = {205, 206, 207};
Circle(123) = {208, 209, 210};
Line Loop(124) = {72, 73, -22, -111, -110};
Plane Surface(125) = {124};
Line Loop(126) = {111, -25, -113, -112};
Plane Surface(127) = {126};
Line Loop(128) = {114, 26, -116, 115};
Plane Surface(129) = {128};

Line Loop(131) = {117, 118, 130, -33};
Plane Surface(132) = {131};
Line Loop(133) = {130, 21, -121, -120, -119};
Plane Surface(134) = {133};
Line Loop(135) = {117, -123, 116, -36};
Ruled Surface(136) = {135};
Line Loop(137) = {122, 114, -34, -113};
Ruled Surface(138) = {137};
Line(139) = {4, 104};
Line(140) = {104, 30};
Line(141) = {204, 9};
Line(142) = {100, 103};
Line(143) = {102, 101};
Line(144) = {101, 7};
Line(145) = {103, 6};
Line Loop(146) = {70, -5, 14, 6, -16, 13, -88, -89, -76};
Plane Surface(147) = {146};
Line Loop(148) = {70, 2, 96, -139, 75};
Plane Surface(149) = {148};
Line Loop(150) = {75, 76, 97, -77};
Plane Surface(151) = {150};
Line Loop(152) = {89, -87, -90, -97};
Plane Surface(153) = {152};
Line Loop(154) = {142, 120, -143, 87};
Plane Surface(155) = {154};
Line Loop(156) = {145, 80, -144, -120};
Plane Surface(157) = {156};
Line Loop(158) = {78, 80, -79, -72};
Plane Surface(159) = {158};
Line Loop(160) = {139, 83, -81, 84, 82, 85, 86, -90, -77};
Plane Surface(161) = {160};
Line Loop(162) = {140, 23, -35, -28, 37, -31, -121, -143, -86, -85, -82, -84, 81, -83};
Plane Surface(163) = {162};
Line Loop(164) = {73, 23, -35, -28, 37, -31, -121, 144, -79};
Plane Surface(165) = {164};
Line Loop(166) = {110, 112, 122, -115, 123, 118, 119, 145, -78};
Plane Surface(167) = {166};
Line Loop(168) = {141, -5, 14, 6, -16, 13, -88, 142, -119, -118, -123, 115, -122, -112};
Plane Surface(169) = {168};
Line Loop(170) = {141, 2, 96, 140, -22, -111};
Plane Surface(171) = {170};


Surface Loop(172) = {149, 147, 153, 161, 151, 109, 103, 107, 101, 105, 99, 39, 47, 207,209,211,213, 49, 43, 45};
Volume(173) = {172};
Surface Loop(174) = {155, 169, 171, 163, 134, 132, 136, 129, 138, 67, 59, 127, 61, 69, 63, 65, 109, 103, 107, 101, 105, 99, 39, 47, 207,209,211,213, 49, 43, 45};
Volume(175) = {174};
Surface Loop(176) = {134, 132, 136, 129, 138, 67, 59, 127, 61, 69, 63, 65, 157, 167, 125, 159, 165};
Volume(177) = {176};

Line Loop(182) = {90, 179, -89, -178};
Ruled Surface(183) = {182};
Line Loop(184) = {143, 180, -142, -179};
Ruled Surface(185) = {184};
Line Loop(186) = {144, 181, -145, -180};
Ruled Surface(187) = {186};
Line Loop(188) = {178, 97};
Plane Surface(189) = {188};
Line Loop(190) = {179, -87};
Plane Surface(191) = {190};
Line Loop(192) = {180, 120};
Plane Surface(193) = {192};
Line Loop(194) = {181, 80};
Plane Surface(195) = {194};
Surface Loop(196) = {189, 183, 153, 191};
Volume(197) = {196};
Surface Loop(198) = {191, 185, 193, 155};
Volume(199) = {198};
Surface Loop(200) = {193, 187, 195, 157};
Volume(201) = {200};


Physical Volume(3) = {173, 175, 177, 197,199,201};
Physical Surface(2) = {45, 43, 49, 207,209,211,213, 47, 39, 65, 63, 69, 61, 67, 59};
Physical Line(1) = {1, 11, 17, 8, 15, 3, 21, 33, 36, 26, 34, 25};




