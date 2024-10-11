h=0.002;

L=0.1;
r=0.005;
e=0.010;
e2=0.015;
r2=0.004;
h2=0.011;

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

Point(5) = {0, 0, L, h};
Point(6) = {L, 0, L, h};
Point(7) = {L, L, L, h};
Point(8) = {0, L, L, h};

// face XZ - Y=0 
Point(10) = {e+r, 0, e+r, h};
Point(11) = {e+r, 0, e, h};
Point(12) = {e, 0, e+r, h};

Point(13) = {L-(e+r), 0, e+r, h};
Point(14) = {L-(e+r), 0, e, h};
Point(15) = {L-e, 0, e+r, h};

Point(16) = {L-(e+r), 0, L-(e+r), h};
Point(17) = {L-(e+r), 0, L-e, h};
Point(18) = {L-e, 0, L-(e+r), h};

Point(19) = {e+r, 0, L-(e+r), h};
Point(20) = {e+r, 0, L-e, h};
Point(21) = {e, 0, L-(e+r), h};

// face YZ - X=0
Point(22) = {0,e+r,  e+r, h};
Point(23) = {0,e+r,  e, h};
Point(24) = {0,e,  e+r, h};

Point(25) = {0,L-(e+r), e+r, h};
Point(26) = {0,L-(e+r), e, h};
Point(27) = {0,L-e, e+r, h};

Point(28) = {0,L-(e+r),  L-(e+r), h};
Point(29) = {0,L-(e+r),  L-e, h};
Point(30) = {0,L-e, L-(e+r), h};

Point(31) = {0,e+r,  L-(e+r), h};
Point(32) = {0,e+r,  L-e, h};
Point(33) = {0,e,  L-(e+r), h};

// face XZ - Y=L
Point(34) = {e+r, L, e+r, h};
Point(35) = {e+r, L, e, h};
Point(36) = {e, L, e+r, h};

Point(37) = {L-(e+r), L, e+r, h};
Point(38) = {L-(e+r), L, e, h};
Point(39) = {L-e, L, e+r, h};

Point(40) = {L-(e+r), L, L-(e+r), h};
Point(41) = {L-(e+r), L, L-e, h};
Point(42) = {L-e, L, L-(e+r), h};

Point(43) = {e+r, L, L-(e+r), h};
Point(44) = {e+r, L, L-e, h};
Point(45) = {e, L, L-(e+r), h};

// face YZ - X=L
Point(46) = {L,e+r,  e+r, h};
Point(47) = {L,e+r,  e, h};
Point(48) = {L,e,  e+r, h};

Point(49) = {L,L-(e+r), e+r, h};
Point(50) = {L,L-(e+r), e, h};
Point(51) = {L,L-e, e+r, h};

Point(52) = {L,L-(e+r),  L-(e+r), h};
Point(53) = {L,L-(e+r),  L-e, h};
Point(54) = {L,L-e, L-(e+r), h};

Point(55) = {L,e+r,  L-(e+r), h};
Point(56) = {L,e+r,  L-e, h};
Point(57) = {L,e,  L-(e+r), h};


// Surfaces des plots elastomere
Point(58) = {e2,0,  0, h};
Point(59) = {0,e2,  0, h};
Point(60) = {e2,e2,  0, h};

Point(61) = {L-e2,0,  0, h};
Point(62) = {L,e2,  0, h};
Point(63) = {L-e2,e2,  0, h};

Point(64) = {e2,L,  0, h};
Point(65) = {0,L-e2,  0, h};
Point(66) = {e2,L-e2,  0, h};

Point(67) = {L-e2,L,  0, h};
Point(68) = {L,L-e2,  0, h};
Point(69) = {L-e2,L-e2,  0, h};

// plots
Point(70) = {e2/2,e2/2,  0, h};
Point(71) = {e2/2+r2,e2/2,  0, h};
Point(72) = {e2/2,e2/2+r2,  0, h};
Point(73) = {e2/2-r2,e2/2,  0, h};
Point(74) = {e2/2,e2/2-r2,  0, h};
Point(75) = {e2/2,e2/2,  -h2, h};
Point(76) = {e2/2+r2,e2/2,  -h2, h};
Point(77) = {e2/2,e2/2+r2,  -h2, h};
Point(78) = {e2/2-r2,e2/2,  -h2, h};
Point(79) = {e2/2,e2/2-r2,  -h2, h};

Point(80) = {L-e2/2,e2/2,  0, h};
Point(81) = {L-(e2/2+r2),e2/2,  0, h};
Point(82) = {L-e2/2,e2/2+r2,  0, h};
Point(83) = {L-(e2/2-r2),e2/2,  0, h};
Point(84) = {L-e2/2,e2/2-r2,  0, h};
Point(85) = {L-e2/2,e2/2,  -h2, h};
Point(86) = {L-(e2/2+r2),e2/2,  -h2, h};
Point(87) = {L-e2/2,e2/2+r2,  -h2, h};
Point(88) = {L-(e2/2-r2),e2/2,  -h2, h};
Point(89) = {L-e2/2,e2/2-r2,  -h2, h};

Point(90) = {e2/2,L-e2/2,  0, h};
Point(91) = {e2/2+r2,L-e2/2,  0, h};
Point(92) = {e2/2,L-(e2/2+r2),  0, h};
Point(93) = {e2/2-r2,L-e2/2,  0, h};
Point(94) = {e2/2,L-(e2/2-r2),  0, h};
Point(95) = {e2/2,L-e2/2,  -h2, h};
Point(96) = {e2/2+r2,L-e2/2,  -h2, h};
Point(97) = {e2/2,L-(e2/2+r2),  -h2, h};
Point(98) = {e2/2-r2,L-e2/2,  -h2, h};
Point(99) = {e2/2,L-(e2/2-r2),  -h2, h};

Point(100) = {L-e2/2,L-e2/2,  0, h};
Point(101) = {L-(e2/2+r2),L-e2/2,  0, h};
Point(102) = {L-e2/2,L-(e2/2+r2),  0, h};
Point(103) = {L-(e2/2-r2),L-e2/2,  0, h};
Point(104) = {L-e2/2,L-(e2/2-r2),  0, h};
Point(105) = {L-e2/2,L-e2/2,  -h2, h};
Point(106) = {L-(e2/2+r2),L-e2/2,  -h2, h};
Point(107) = {L-e2/2,L-(e2/2+r2),  -h2, h};
Point(108) = {L-(e2/2-r2),L-e2/2,  -h2, h};
Point(109) = {L-e2/2,L-(e2/2-r2),  -h2, h};

Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {8, 5};
Line(10) = {5, 6};
Line(11) = {6, 7};
Line(12) = {7, 8};
Line(13) = {11, 14};
Line(14) = {15, 18};
Line(15) = {17, 20};
Line(16) = {21, 12};
Line(17) = {24, 33};
Line(18) = {32, 29};
Line(19) = {30, 27};
Line(20) = {26, 23};
Line(21) = {48, 57};
Line(22) = {56, 53};
Line(23) = {54, 51};
Line(24) = {50, 47};
Line(25) = {42, 39};
Line(26) = {38, 35};
Line(27) = {36, 45};
Line(28) = {44, 41};
Circle(29) = {38, 37, 39};
Circle(30) = {42, 40, 41};
Circle(31) = {44, 43, 45};
Circle(32) = {36, 34, 35};
Circle(33) = {29, 28, 30};
Circle(34) = {27, 25, 26};
Circle(35) = {23, 22, 24};
Circle(36) = {33, 31, 32};
Circle(37) = {20, 19, 21};
Circle(38) = {12, 10, 11};
Circle(39) = {14, 13, 15};
Circle(40) = {18, 16, 17};
Circle(41) = {56, 55, 57};
Circle(42) = {48, 46, 47};
Circle(43) = {50, 49, 51};
Circle(44) = {54, 52, 53};


Line(45) = {1, 59};
Line(46) = {59, 60};
Line(47) = {60, 58};
Line(48) = {58, 1};
Line(49) = {59, 65};
Line(50) = {65, 66};
Line(51) = {66, 64};
Line(52) = {64, 4};
Line(53) = {4, 65};
Line(54) = {64, 67};
Line(55) = {67, 69};
Line(56) = {69, 68};
Line(57) = {68, 3};
Line(58) = {3, 67};
Line(59) = {68, 62};
Line(60) = {62, 2};
Line(61) = {62, 63};
Line(62) = {63, 61};
Line(63) = {61, 2};
Line(64) = {58, 61};
Line(65) = {83, 88};
Line(66) = {82, 87};
Line(67) = {81, 86};
Line(68) = {84, 89};
Line(69) = {103, 108};
Line(70) = {104, 109};
Line(71) = {101, 106};
Line(72) = {102, 107};
Line(73) = {91, 96};
Line(74) = {94, 99};
Line(75) = {92, 97};
Line(76) = {93, 98};
Line(77) = {72, 77};
Line(78) = {73, 78};
Line(79) = {71, 76};
Line(80) = {74, 79};
Circle(81) = {74, 70, 71};
Circle(82) = {71, 70, 72};
Circle(83) = {72, 70, 73};
Circle(84) = {73, 70, 74};
Circle(85) = {79, 75, 76};
Circle(86) = {76, 75, 77};
Circle(87) = {77, 75, 78};
Circle(88) = {78, 75, 79};
Circle(89) = {92, 90, 91};
Circle(90) = {91, 90, 94};
Circle(91) = {94, 90, 93};
Circle(92) = {93, 90, 92};
Circle(93) = {97, 95, 96};
Circle(94) = {96, 95, 99};
Circle(95) = {99, 95, 98};
Circle(96) = {98, 95, 97};
Circle(97) = {102, 100, 103};
Circle(98) = {103, 100, 104};
Circle(99) = {104, 100, 101};
Circle(100) = {101, 100, 102};
Circle(101) = {107, 105, 108};
Circle(102) = {108, 105, 109};
Circle(103) = {109, 105, 106};
Circle(104) = {106, 105, 107};
Circle(105) = {89, 85, 88};
Circle(106) = {88, 85, 87};
Circle(107) = {87, 85, 86};
Circle(108) = {86, 85, 89};
Circle(109) = {84, 80, 83};
Circle(110) = {83, 80, 82};
Circle(111) = {82, 80, 81};
Circle(112) = {81, 80, 84};
Line Loop(113) = {12, -8, -52, 54, -58, 7};
Line Loop(114) = {25, -29, 26, -32, 27, -31, 28, -30};
Plane Surface(115) = {113, 114};
Line Loop(116) = {11, -7, -57, 59, 60, 6};
Line Loop(117) = {23, -43, 24, -42, 21, -41, 22, -44};
Plane Surface(118) = {116, 117};
Line Loop(119) = {10, -6, -63, -64, 48, 5};
Line Loop(120) = {14, 40, 15, 37, 16, 38, 13, 39};
Plane Surface(121) = {119, 120};
Line Loop(122) = {9, -5, 45, 49, -53, 8};
Line Loop(123) = {18, 33, 19, 34, 20, 35, 17, 36};
Plane Surface(124) = {122, 123};
Line Loop(125) = {56, 57, 58, 55};
Line Loop(126) = {98, 99, 100, 97};
Plane Surface(127) = {125, 126};
Plane Surface(128) = {126};
Line Loop(129) = {60, -63, -62, -61};
Line Loop(130) = {112, 109, 110, 111};
Plane Surface(131) = {129, 130};
Line Loop(132) = {48, 45, 46, 47};
Line Loop(133) = {82, 83, 84, 81};
Plane Surface(134) = {132, 133};
Plane Surface(135) = {133};
Line Loop(136) = {52, 53, 50, 51};
Line Loop(137) = {89, 90, 91, 92};
Plane Surface(138) = {136, 137};
Plane Surface(139) = {137};
Plane Surface(140) = {130};
Line Loop(141) = {103, 104, 101, 102};
Plane Surface(142) = {141};
Line Loop(143) = {94, 95, 96, 93};
Plane Surface(144) = {143};
Line Loop(145) = {87, 88, 85, 86};
Plane Surface(146) = {145};
Line Loop(147) = {107, 108, 105, 106};
Plane Surface(148) = {147};
Line Loop(149) = {107, -67, -111, 66};
Ruled Surface(150) = {149};
Line Loop(151) = {67, 108, -68, -112};
Ruled Surface(152) = {151};
Line Loop(153) = {68, 105, -65, -109};
Ruled Surface(154) = {153};
Line Loop(155) = {106, -66, -110, 65};
Ruled Surface(156) = {155};
Line Loop(157) = {103, -71, -99, 70};
Ruled Surface(158) = {157};
Line Loop(159) = {104, -72, -100, 71};
Ruled Surface(160) = {159};
Line Loop(161) = {97, 69, -101, -72};
Ruled Surface(162) = {161};
Line Loop(163) = {102, -70, -98, 69};
Ruled Surface(164) = {163};
Line Loop(165) = {95, -76, -91, 74};
Ruled Surface(166) = {165};
Line Loop(167) = {76, 96, -75, -92};
Ruled Surface(168) = {167};
Line Loop(169) = {75, 93, -73, -89};
Ruled Surface(170) = {169};
Line Loop(171) = {94, -74, -90, 73};
Ruled Surface(172) = {171};
Line Loop(173) = {87, -78, -83, 77};
Ruled Surface(174) = {173};
Line Loop(175) = {88, -80, -84, 78};
Ruled Surface(176) = {175};
Line Loop(177) = {85, -79, -81, 80};
Ruled Surface(178) = {177};
Line Loop(179) = {86, -77, -82, 79};
Ruled Surface(180) = {179};


Surface Loop(181) = {146, 174, 176, 178, 180, 135};
Volume(182) = {181};
Surface Loop(183) = {148, 150, 152, 154, 156, 140};
Volume(184) = {183};
Surface Loop(185) = {164, 142, 158, 160, 162, 128};
Volume(186) = {185};
Surface Loop(187) = {172, 144, 166, 168, 170, 139};
Volume(188) = {187};

// surfaces superieures des plots
Physical Surface(1) = {140};
Physical Surface(2) = {128};
Physical Surface(3) = {139};
Physical Surface(4) = {135};

// surfaces inferieures des plots (meme ordre que les superieures)
Physical Surface(5) = {148};
Physical Surface(6) = {142};
Physical Surface(7) = {144};
Physical Surface(8) = {146};

// surface du cadre de la structure du cubesat
Physical Surface(10) = {131, 127, 138, 134, 118, 115, 124, 121};

// volumes des plots (meme ordre que les surfaces)

Physical Volume(1) = {184};
Physical Volume(2) = {186};
Physical Volume(3) = {188};
Physical Volume(4) = {182};

// toutes les surfaces pour la visu
Physical Surface(20) = {164, 162, 128, 158, 127, 142, 160, 118, 156, 154, 140, 150, 131, 148, 152, 115, 172, 170, 139, 166, 138, 144, 168, 121, 124, 180, 178, 135, 174, 134, 146, 176};
