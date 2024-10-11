h=3;

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
Point(54) = {r1+e+lx,	2*ly,	e,		h};
Point(55) = {r1+e+lx,	2*ly,	0,		h};
Point(56) = {e,		2*ly,	r1+e,		h};
Point(57) = {0,		2*ly,	r1+e,		h};
Point(58) = {e,		2*ly,	r1+e+lz-s1-r3,	h};
Point(59) = {0,		2*ly,	r1+e+lz-s1-r3,	h};

Point(61) = {r1+e+lx-r2,2*ly-r2,	0,	h};
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

Point(138) = {e,	2*ly,	r1+e+lz-s1,	h};
Point(139) = {0,	2*ly,	r1+e+lz-s1,	h};
Point(140) = {e,	2*ly,	r1+e+lz-s1+r3,	h};
Point(141) = {0,	2*ly,	r1+e+lz-s1+r3,	h};
Point(142) = {e,	2*ly-r3,	r1+e+lz-s1,	h};
Point(143) = {0,	2*ly-r3,	r1+e+lz-s1,	h};
Point(144) = {e,	2*ly,	a+rayon,	h};
Point(145) = {0,	2*ly,	a+rayon,	h};


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
Line(41) = {55, 24};
Line(42) = {54, 25};
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
Line(192) = {3, 130};
Line(193) = {130, 7};
Line(194) = {130, 131};
Line(195) = {131, 13};
Line(196) = {131, 17};
Line(197) = {131, 132};
Line(198) = {132, 53};
Line(199) = {132, 57};

Circle(256) = {58, 138, 142};
Circle(257) = {142, 138, 140};
Circle(258) = {59, 139, 143};
Circle(259) = {143, 139, 141};
Line(260) = {140, 144};
Line(261) = {141, 145};
Line(262) = {140, 141};
Line(263) = {145, 144};
Line(264) = {144, 30};
Line(265) = {40, 145};
Line(276) = {142, 143};

// 1ere partie

Line Loop(727) = {3, 26, 30, 46, 47, -12, -64};
Plane Surface(728) = {727};
Line Loop(729) = {4, -61, 14, -51, -50, -29, -27};
Plane Surface(730) = {729};
Line Loop(731) = {4, 161, 3, 7};
Plane Surface(732) = {731};
Line Loop(735) = {2, 58, -10, 21, -24, -19};
Line Loop(736) = {71, 72, 73, 70};
Plane Surface(737) = {735, 736};
Line Loop(738) = {62, -9, 20, 25, -18, 1};
Line Loop(739) = {68, 69, 66, 67};
Plane Surface(740) = {738, 739};
Line Loop(741) = {25, 23, 24, 22};
Plane Surface(742) = {741};
Line Loop(743) = {27, 28, -26, 7};
Ruled Surface(744) = {743};
Line Loop(745) = {50, 55, -46, -54};
Ruled Surface(746) = {745};
Line Loop(747) = {55, 47, -56, -51};
Ruled Surface(748) = {747};
Line Loop(749) = {68, 136, -72, -133};
Ruled Surface(750) = {749};
Line Loop(751) = {69, 135, -73, -136};
Ruled Surface(752) = {751};
Line Loop(753) = {66, 134, -70, -135};
Ruled Surface(754) = {753};
Line Loop(755) = {67, 133, -71, -134};
Ruled Surface(756) = {755};
Line Loop(757) = {62, 16, -64, -5};
Ruled Surface(758) = {757};
Line Loop(759) = {162, 5, -161, -193, -192};
Plane Surface(760) = {759};
Line Loop(761) = {58, -195, -194, -192};
Plane Surface(762) = {761};
Line Loop(911) = {18, 23, -19, 8};
Ruled Surface(912) = {911};
Line Loop(1061) = {1, -162, -2, 8};
Plane Surface(1062) = {1061};
Line Loop(1063) = {193, -61, -196, -194};
Plane Surface(1064) = {1063};

// 2eme partie


Line Loop(215) = {63, 35, -65, -16};
Ruled Surface(216) = {215};
Line Loop(217) = {129, 76, -132, -80};
Ruled Surface(218) = {217};
Line Loop(219) = {132, 77, -131, -81};
Ruled Surface(220) = {219};
Line Loop(221) = {78, 130, -74, -131};
Ruled Surface(222) = {221};
Line Loop(223) = {130, 75, -129, -79};
Ruled Surface(224) = {223};
Line Loop(202) = {31, 42, -20, 9, 63};
Line Loop(203) = {78, 79, 80, 81};
Plane Surface(204) = {202, 203};
Line Loop(227) = {21, -41, 32, -59, -10};
Line Loop(228) = {75, 76, 77, 74};
Plane Surface(229) = {227, 228};
Line Loop(232) = {198, -59, -195, 197};
Plane Surface(233) = {232};
Line Loop(234) = {197, 199, 60, -196};
Plane Surface(235) = {234};
Line Loop(250) = {53, 54, -49, -57};
Ruled Surface(251) = {250};
Line Loop(252) = {57, -48, -56, 52};
Ruled Surface(253) = {252};
Line Loop(254) = {22, -42, 44, 41};
Plane Surface(255) = {254};


Line Loop(266) = {12, 48, 49, -30, -264, -260, -257, -256, 34, -65};
Plane Surface(267) = {266};
Line Loop(268) = {60, 14, 52, 53, -29, 265, -261, -259, -258, -33};
Plane Surface(269) = {268};
Line Loop(270) = {265, 263, 264, -28};
Plane Surface(271) = {270};

Line Loop(277) = {256, 276, -258, -45};
Ruled Surface(278) = {277};
Line Loop(279) = {276, 259, -262, -257};
Ruled Surface(280) = {279};


Symmetry {0, 1, 0, -45} {
  Duplicata { Surface{732, 760, 1062, 912, 752, 750}; }
}
