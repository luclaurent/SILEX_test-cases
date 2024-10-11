h=1.5;

r1=5;
e=10;
lx=30;
lz=25;
r2=15;
ly=30;
a=-15;
r3=7;
r4=5;
rayon=Sqrt(ly*ly+(r1+e+lz-a)*(r1+e+lz-a));
s1=5;

Point(1) = {r1+e,	0, 	r1+e,		h};
Point(2) = {r1+e,	0, 	e,		h};
Point(3) = {r1+e,	0,	0,		h};
Point(4) = {r1+e+lx-r2,	0,	e,		h};
Point(5) = {r1+e+lx-r2,	0,	0,		h};
Point(6) = {e,		0,	r1+e,		h};
Point(7) = {0,		0,	r1+e,		h};
Point(8) = {e,		0,	r1+e+lz,	h};
Point(9) = {0,		0,	r1+e+lz,	h};

Point(11) = {r1+e,	ly, 	r1+e,		h};
Point(12) = {r1+e,	ly, 	e,		h};
Point(13) = {r1+e,	ly,	0,		h};
Point(14) = {r1+e+lx-r2,ly,	e,		h};
Point(15) = {r1+e+lx-r2,ly,	0,		h};
Point(16) = {e,		ly,	r1+e,		h};
Point(17) = {0,		ly,	r1+e,		h};
Point(18) = {e,		ly,	r1+e+lz-s1,	h};
Point(19) = {0,		ly,	r1+e+lz-s1,	h};

Point(20) = {r1+e+lx,	r2,	0,		h};
Point(21) = {r1+e+lx-r2,r2,	0,		h};
Point(22) = {r1+e+lx,	r2,	e,		h};
Point(23) = {r1+e+lx-r2,r2,	e,		h};

Point(24) = {r1+e+lx,	ly,	0,		h};
Point(25) = {r1+e+lx,	ly,	e,		h};

Point(30) = {e,		ly,	a+rayon,	h};
Point(33) = {e,		ly,	a,		h};
Point(40) = {0,		ly,	a+rayon,	h};
Point(43) = {0,		ly,	a,		h};

Point(51) = {r1+e,	2*ly, 	r1+e,		h};
Point(52) = {r1+e,	2*ly, 	e,		h};
Point(53) = {r1+e,	2*ly,	0,		h};
Point(54) = {r1+e+lx-r2,2*ly,	e,		h};
Point(55) = {r1+e+lx-r2,2*ly,	0,		h};
Point(56) = {e,		2*ly,	r1+e,		h};
Point(57) = {0,		2*ly,	r1+e,		h};
Point(58) = {e,		2*ly,	r1+e+lz,	h};
Point(59) = {0,		2*ly,	r1+e+lz,	h};

Point(60) = {r1+e+lx,	2*ly-r2,	0,	h};
Point(61) = {r1+e+lx-r2,2*ly-r2,	0,	h};
Point(62) = {r1+e+lx,	2*ly-r2,	e,	h};
Point(63) = {r1+e+lx-r2,2*ly-r2,	e,	h};

Point(70) = {e,		ly+r3,	r1+e+lz-s1,	h};
Point(71) = {e,		ly,	r1+e+lz-s1+r3,	h};
Point(72) = {e,		ly,	r1+e+lz-s1-r3,	h};
Point(73) = {e,		ly-r3,	r1+e+lz-s1,	h};

Point(80) = {0,		ly+r3,	r1+e+lz-s1,	h};
Point(81) = {0,		ly,	r1+e+lz-s1+r3,	h};
Point(82) = {0,		ly,	r1+e+lz-s1-r3,	h};
Point(83) = {0,		ly-r3,	r1+e+lz-s1,	h};

Point(90) = {r1+e+lx-r2+r4,r2,	e,		h};
Point(91) = {r1+e+lx-r2-r4,r2,	e,		h};
Point(92) = {r1+e+lx-r2,r2+r4,	e,		h};
Point(93) = {r1+e+lx-r2,r2-r4,	e,		h};

Point(100) = {r1+e+lx-r2+r4,r2,	0,		h};
Point(101) = {r1+e+lx-r2-r4,r2,	0,		h};
Point(102) = {r1+e+lx-r2,r2+r4,	0,		h};
Point(103) = {r1+e+lx-r2,r2-r4,	0,		h};

Point(110) = {r1+e+lx-r2+r4,2*ly-r2,	0,	h};
Point(111) = {r1+e+lx-r2-r4,2*ly-r2,	0,	h};
Point(112) = {r1+e+lx-r2,2*ly-r2+r4,	0,	h};
Point(113) = {r1+e+lx-r2,2*ly-r2-r4,	0,	h};

Point(120) = {r1+e+lx-r2+r4,2*ly-r2,	e,	h};
Point(121) = {r1+e+lx-r2-r4,2*ly-r2,	e,	h};
Point(122) = {r1+e+lx-r2,2*ly-r2+r4,	e,	h};
Point(123) = {r1+e+lx-r2,2*ly-r2-r4,	e,	h};

Point(130) = {0,	0,	0,	h};
Point(131) = {0,	ly,	0,	h};
Point(132) = {0,	2*ly,	0,	h};


Line(1) = {4, 2};
Line(2) = {5, 3};
Line(3) = {6, 8};
Line(4) = {9, 7};
Circle(5) = {2, 1, 6};
Line(7) = {8, 9};
Line(8) = {5, 4};

Line(9) = {14, 12};
Line(10) = {15, 13};
Line(12) = {16, 72};
Line(13) = {19, 19};
Line(14) = {17, 82};
Line(15) = {19, 18};
Circle(16) = {12, 11, 16};
Circle(18) = {4, 23, 22};
Circle(19) = {5, 21, 20};
Line(20) = {14, 25};
Line(21) = {15, 24};
Line(22) = {24, 25};
Line(23) = {22, 20};
Line(24) = {20, 24};
Line(25) = {25, 22};

Circle(26) = {8, 33, 30};
Circle(27) = {9, 43, 40};

Line(28) = {40, 30};
Line(29) = {40, 81};
Line(30) = {30, 71};
Line(31) = {52, 54};
Line(32) = {55, 53};
Line(33) = {57, 59};
Line(34) = {58, 56};
Circle(35) = {52, 51, 56};
Circle(37) = {58, 33, 30};
Circle(38) = {59, 43, 40};
Circle(39) = {55, 61, 60};
Circle(40) = {54, 63, 62};
Line(41) = {60, 24};
Line(42) = {62, 25};
Line(43) = {62, 60};
Line(44) = {54, 55};
Line(45) = {58, 59};
Circle(46) = {71, 18, 73};
Circle(47) = {73, 18, 72};
Circle(48) = {72, 18, 70};
Circle(49) = {70, 18, 71};
Circle(50) = {81, 19, 83};
Circle(51) = {83, 19, 82};
Circle(52) = {82, 19, 80};
Circle(53) = {80, 19, 81};
Line(54) = {81, 71};
Line(55) = {83, 73};
Line(56) = {82, 72};
Line(57) = {80, 70};
Line(58) = {3, 13};
Line(59) = {13, 53};
Line(60) = {57, 17};
Line(61) = {17, 7};
Line(62) = {2, 12};
Line(63) = {12, 52};
Line(64) = {6, 16};
Line(65) = {16, 56};
Circle(66) = {93, 23, 90};
Circle(67) = {90, 23, 92};
Circle(68) = {92, 23, 91};
Circle(69) = {91, 23, 93};
Circle(70) = {103, 21, 100};
Circle(71) = {100, 21, 102};
Circle(72) = {102, 21, 101};
Circle(73) = {101, 21, 103};

Circle(74) = {113, 61, 110};
Circle(75) = {110, 61, 112};
Circle(76) = {112, 61, 111};
Circle(77) = {111, 61, 113};

Circle(78) = {123, 63, 120};
Circle(79) = {120, 63, 122};
Circle(80) = {122, 63, 121};
Circle(81) = {121, 63, 123};

Line(129) = {122, 112};
Line(130) = {120, 110};
Line(131) = {123, 113};
Line(132) = {121, 111};
Line(133) = {92, 102};
Line(134) = {90, 100};
Line(135) = {93, 103};
Line(136) = {91, 101};

Line(161) = {7, 6};
Line(162) = {3, 2};
Line(163) = {13, 12};
Line(164) = {16, 17};
Line(165) = {53, 52};
Line(166) = {57, 56};

Line Loop(84) = {12, 48, 49, -30, -37, 34, -65};
Plane Surface(85) = {84};
Line Loop(92) = {14, 52, 53, -29, -38, -33, 60};
Plane Surface(93) = {92};
Line Loop(94) = {1, 62, -9, 20, 25, -18};
Line Loop(95) = {68, 69, 66, 67};
Plane Surface(96) = {94, 95};
Line Loop(97) = {63, 31, 40, 42, -20, 9};
Line Loop(98) = {80, 81, 78, 79};
Plane Surface(99) = {97, 98};
Line Loop(100) = {59, -32, 39, 41, -21, 10};
Line Loop(101) = {76, 77, 74, 75};
Plane Surface(102) = {100, 101};
Line Loop(105) = {18, 23, -19, 8};
Ruled Surface(106) = {105};
Line Loop(107) = {43, -39, -44, 40};
Ruled Surface(108) = {107};
Line Loop(109) = {35, -65, -16, 63};
Ruled Surface(110) = {109};
Line Loop(111) = {16, -64, -5, 62};
Ruled Surface(112) = {111};
Line Loop(113) = {27, 28, -26, 7};
Ruled Surface(114) = {113};
Line Loop(115) = {28, -37, 45, 38};
Ruled Surface(116) = {115};
Line Loop(117) = {53, 54, -49, -57};
Ruled Surface(118) = {117};
Line Loop(119) = {50, 55, -46, -54};
Ruled Surface(120) = {119};
Line Loop(121) = {55, 47, -56, -51};
Ruled Surface(122) = {121};
Line Loop(123) = {56, 48, -57, -52};
Ruled Surface(124) = {123};

Line Loop(137) = {76, -132, -80, 129};
Ruled Surface(138) = {137};
Line Loop(139) = {129, -75, -130, 79};
Ruled Surface(140) = {139};
Line Loop(141) = {130, -74, -131, 78};
Ruled Surface(142) = {141};
Line Loop(143) = {81, 131, -77, -132};
Ruled Surface(144) = {143};
Line Loop(145) = {71, -133, -67, 134};
Ruled Surface(146) = {145};
Line Loop(147) = {68, 136, -72, -133};
Ruled Surface(148) = {147};
Line Loop(149) = {136, 73, -135, -69};
Ruled Surface(150) = {149};
Line Loop(151) = {135, 70, -134, -66};
Ruled Surface(152) = {151};
Line Loop(153) = {29, 54, -30, -28};
Plane Surface(154) = {153};
Line Loop(157) = {23, 24, 22, 25};
Plane Surface(158) = {157};
Line Loop(159) = {22, -42, 43, 41};
Plane Surface(160) = {159};

Line Loop(167) = {7, 4, 161, 3};
Plane Surface(168) = {167};

Line Loop(171) = {2, 162, -1, -8};
Plane Surface(172) = {171};
Line Loop(173) = {14, 56, -12, 164};
Plane Surface(174) = {173};
Line Loop(177) = {163, -9, 20, -22, -21, 10};
Plane Surface(178) = {177};
Line Loop(179) = {32, 165, 31, 44};
Plane Surface(180) = {179};
Line Loop(183) = {166, -34, 45, -33};
Plane Surface(184) = {183};

Line Loop(185) = {2, 58, -10, 21, -24, -19};
Line Loop(186) = {71, 72, 73, 70};
Plane Surface(187) = {185, 186};
Line Loop(188) = {64, 12, -47, -46, -30, -26, -3};
Plane Surface(189) = {188};
Line Loop(190) = {61, -4, 27, 29, 50, 51, -14};
Plane Surface(191) = {190};


Line(192) = {3, 130};
Line(193) = {130, 7};
Line(194) = {130, 131};
Line(195) = {131, 13};
Line(196) = {131, 17};
Line(197) = {131, 132};
Line(198) = {132, 53};
Line(199) = {132, 57};
Line Loop(200) = {5, -161, -193, -192, 162};
Plane Surface(201) = {200};
Line Loop(202) = {58, -195, -194, -192};
Plane Surface(203) = {202};
Line Loop(204) = {195, 163, 16, 164, -196};
Plane Surface(205) = {204};
Line Loop(206) = {195, 59, -198, -197};
Plane Surface(207) = {206};
Line Loop(208) = {198, 165, 35, -166, -199};
Plane Surface(209) = {208};
Line Loop(210) = {60, -196, 197, 199};
Plane Surface(211) = {210};
Line Loop(212) = {196, 61, -193, 194};
Plane Surface(213) = {212};

Physical Surface(2) = {187, 102, 203, 207};
Physical Surface(4) = {138, 140, 142, 144, 148, 146, 150, 152};
Physical Surface(3) = {124, 118, 120, 122};


Surface Loop(214) = {189, 112, 201, 168, 114, 191, 213, 203, 187, 172, 96, 158, 106, 146, 148, 150, 152, 120, 122, 154, 174, 205, 178};
Volume(215) = {214};
Surface Loop(216) = {85, 124, 118, 93, 116, 184, 209, 207, 102, 180, 99, 110, 108, 160, 140, 138, 144, 142, 211, 174, 205, 178, 154};
Volume(217) = {216};
Physical Volume(1) = {215};
Physical Volume(2) = {217};
Physical Surface(5) = {174, 205, 178, 154};

