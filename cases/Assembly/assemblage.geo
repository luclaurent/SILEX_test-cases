
//================================================
//                                               
//    MAILLAGE D'UNE LIAISON DIABOLO PAULSTRA       
//                                               
//================================================


// PARAMETRES


// NB NOEUDS LIAISON

// maillage 1

//nl = 5;           // moitié du nombre de noeuds dans la liaison suivant l'axe Y
//np = 2;           // nombre de noeuds dans les plaques suivant l'axe Y
//nc = 5;           // nombre de noeuds pour les carrés et les arcs de cercle
//nt = 6;           // nombre de noeuds pour les lignes entre les carrés et les arcs de cercle

// maillage 2

nl = 11;           // moitié du nombre de noeuds dans la liaison suivant l'axe Y
np = 2;           // nombre de noeuds dans les plaques suivant l'axe Y
nc = 6;           // nombre de noeuds pour les carrés et les arcs de cercle
nt = 6;           // nombre de noeuds pour les lignes entre les carrés et les arcs de cercle

// maillage 3

//nl = 14;           // moitié du nombre de noeuds dans la liaison suivant l'axe Y
//np = 2;           // nombre de noeuds dans les plaques suivant l'axe Y
//nc = 8;           // nombre de noeuds pour les carrés et les arcs de cercle
//nt = 8;           // nombre de noeuds pour les lignes entre les carrés et les arcs de cercle

// maillage 4

//nl = 16;           // moitié du nombre de noeuds dans la liaison suivant l'axe Y
//np = 2;           // nombre de noeuds dans les plaques suivant l'axe Y
//nc = 10;           // nombre de noeuds pour les carrés et les arcs de cercle
//nt = 10;           // nombre de noeuds pour les lignes entre les carrés et les arcs de cercle

// NB NOEUDS STRUCTURE SUPPORT

// maillage 1

//nh = 3; // nb noeuds sur plaque superieure direction radiale et transversale 
//ng = 3; // nb noeuds sur plaque superieure direction radiale 2
//nep = 2; // nb noeuds dans l'epaisseur de la plaque
//nep2 = 2; // nb noeuds dans l'epaisseur des montants
//n3 = 7; // nb noeuds dans la hauteur des montants
//n4 = 2; // nb noeuds dans l'epaisseur de la plaque superieure
//n5 = 4; // nb noeuds dans la longueur des coins superieurs
//n6 = 5; // nb noeuds dans la longueur parties droites des plaques superieure et inferieure

// maillage 2

nh = 6; // nb noeuds sur plaque superieure direction radiale et transversale 
ng = 5; // nb noeuds sur plaque superieure direction radiale 2
nep = 3; // nb noeuds dans l'epaisseur de la plaque
nep2 = 3; // nb noeuds dans l'epaisseur des montants
n3 = 13; // nb noeuds dans la hauteur des montants
n4 = 3; // nb noeuds dans l'epaisseur de la plaque superieure
n5 = 9; // nb noeuds dans la longueur des coins superieurs
n6 = 10; // nb noeuds dans la longueur parties droites des plaques superieure et inferieure

// maillage 3

//nh = 9; // nb noeuds sur plaque superieure direction radiale et transversale 
//ng = 9; // nb noeuds sur plaque superieure direction radiale 2
//nep = 4; // nb noeuds dans l'epaisseur de la plaque
//nep2 = 4; // nb noeuds dans l'epaisseur des montants
//n3 = 18; // nb noeuds dans la hauteur des montants
//n4 = 4; // nb noeuds dans l'epaisseur de la plaque superieure
//n5 = 11; // nb noeuds dans la longueur des coins superieurs
//n6 = 14; // nb noeuds dans la longueur parties droites des plaques superieure et inferieure

// maillage 4

//nh = 10; // nb noeuds sur plaque superieure direction radiale et transversale 
//ng = 10; // nb noeuds sur plaque superieure direction radiale 2
//nep = 4; // nb noeuds dans l'epaisseur de la plaque
//nep2 = 4; // nb noeuds dans l'epaisseur des montants
//n3 = 21; // nb noeuds dans la hauteur des montants
//n4 = 4; // nb noeuds dans l'epaisseur de la plaque superieure
//n5 = 15; // nb noeuds dans la longueur des coins superieurs
//n6 = 21; // nb noeuds dans la longueur parties droites des plaques superieure et inferieure

DA = 100.0e-3;    // diamètre des faces inf/sup
A = DA/2.0;       // rayon des faces  inf/sup
Ap = A/2.0;       // A prime
As = A/3.0;       // A second
B = 100.0e-3;     // hauteur de la liaison, sans les plaques en metal
H = B/2.0;        // demi-hauteur de la liaison
DC = 70.0e-3;     // diamètre de la face au milieu de la liaison
C = DC/2.0;       // rayon de la face au milieu de la liaison
E = 2.0e-3;       // epaisseur des plaques dessus et dessous

L1 = 600.0e-3;    // entraxe entre les diabolos
L2 = 1000.0e-3;    // entraxe entre les diabolos
R = 100.0e-3;    // rayon des 4 coins
ep = 30.0e-3;     // epaisseur de la plaque
ep2 = 30.0e-3;     // epaisseur partie verticale de la structure
HH = 900.0e-3; // hauteur des montants
ep3 = 30.0e-3; // epaisseur plaque superieure de la structure

// cf fichier rayon de courbure dans répertoire parent:

X = A-C;
R3 = (X^2+(B/2.0)^2)/(2.0*X);
D = R3+C;          // coordonnée des noeuds centre de rotation pour
                  // les rayons de courbures

Geometry.CopyMeshingMethod = 1;

// DEFINITION DE LA GEOMETRIE
//===========================

// liste des points
Point(1) = {0, 0, 0, 1};
Point(2) = {A, 0, 0, 1};
Point(3) = {0, 0, -A, 1};
Point(4) = {-A, 0, 0, 1};
Point(5) = {0, 0, A, 1};

Point(6) = {0, E, 0, 1};
Point(7) = {A, E, 0, 1};
Point(8) = {0, E, -A, 1};
Point(9) = {-A, E, 0, 1};
Point(10) = {0, E, A, 1};

Point(11) = {0, E+B, 0, 1};
Point(12) = {A, E+B, 0, 1};
Point(13) = {0, E+B, -A, 1};
Point(14) = {-A, E+B, 0, 1};
Point(15) = {0, E+B, A, 1};

Point(16) = {0, 2*E+B, 0, 1};
Point(17) = {A, 2*E+B, 0, 1};
Point(18) = {0, 2*E+B, -A, 1};
Point(19) = {-A, 2*E+B, 0, 1};
Point(20) = {0, 2*E+B, A, 1};

Point(21) = {0, E+H, 0, 1};
Point(22) = {C, E+H, 0, 1};
Point(23) = {0, E+H, -C, 1};
Point(24) = {-C, E+H, 0, 1};
Point(25) = {0, E+H, C, 1};

Point(31) = {D, E+H, 0, 1};
Point(32) = {0, E+H, -D, 1};
Point(33) = {-D, E+H, 0, 1};
Point(34) = {0, E+H, D, 1};

Point(35) = {Ap, 0, 0, 1};
Point(36) = {0, 0, -Ap, 1};
Point(37) = {-Ap, 0, 0, 1};
Point(38) = {0, 0, Ap, 1};

Point(39) = {Ap, E, 0, 1};
Point(40) = {0, E, -Ap, 1};
Point(41) = {-Ap, E, 0, 1};
Point(42) = {0, E, Ap, 1};

Point(43) = {Ap, E+B, 0, 1};
Point(44) = {0, E+B, -Ap, 1};
Point(45) = {-Ap, E+B, 0, 1};
Point(46) = {0, E+B, Ap, 1};

Point(47) = {Ap, 2*E+B, 0, 1};
Point(48) = {0, 2*E+B, -Ap, 1};
Point(49) = {-Ap, 2*E+B, 0, 1};
Point(50) = {0, 2*E+B, Ap, 1};

Point(51) = {As, E+H, 0, 1};
Point(52) = {0, E+H, -As, 1};
Point(53) = {-As, E+H, 0, 1};
Point(54) = {0, E+H, As, 1};

// liste des lignes
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};
Circle(9) = {12, 11, 13};
Circle(10) = {13, 11, 14};
Circle(11) = {14, 11, 15};
Circle(12) = {15, 11, 12};
Circle(13) = {17, 16, 18};
Circle(14) = {18, 16, 19};
Circle(15) = {19, 16, 20};
Circle(16) = {20, 16, 17};
Circle(17) = {22, 21, 23};
Circle(18) = {23, 21, 24};
Circle(19) = {24, 21, 25};
Circle(20) = {25, 21, 22};
Circle(25) = {12, 31, 22};
Circle(26) = {22, 31, 7};
Circle(27) = {13, 32, 23};
Circle(28) = {23, 32, 8};
Circle(29) = {14, 33, 24};
Circle(30) = {24, 33, 9};
Circle(31) = {15, 34, 25};
Circle(32) = {25, 34, 10};
Line(33) = {17, 12};
Line(34) = {18, 13};
Line(35) = {19, 14};
Line(36) = {20, 15};
Line(37) = {7, 2};
Line(38) = {8, 3};
Line(39) = {9, 4};
Line(40) = {10, 5};
Line(105) = {36, 37};
Line(106) = {37, 38};
Line(107) = {38, 35};
Line(108) = {35, 36};
Line(109) = {40, 41};
Line(110) = {41, 42};
Line(111) = {42, 39};
Line(112) = {39, 40};
Line(113) = {44, 45};
Line(114) = {45, 46};
Line(115) = {46, 43};
Line(116) = {43, 44};
Line(117) = {48, 49};
Line(118) = {49, 50};
Line(119) = {50, 47};
Line(120) = {47, 48};
Line(121) = {52, 53};
Line(122) = {53, 54};
Line(123) = {54, 51};
Line(124) = {51, 52};
Line(125) = {36, 40};
Line(126) = {37, 41};
Line(127) = {38, 42};
Line(128) = {35, 39};
Line(129) = {44, 48};
Line(130) = {45, 49};
Line(131) = {46, 50};
Line(132) = {43, 47};
Line(133) = {40, 52};
Line(134) = {41, 53};
Line(135) = {42, 54};
Line(136) = {39, 51};
Line(137) = {52, 44};
Line(138) = {53, 45};
Line(139) = {54, 46};
Line(140) = {51, 43};
Line(141) = {36, 3};
Line(142) = {37, 4};
Line(143) = {38, 5};
Line(144) = {35, 2};
Line(145) = {40, 8};
Line(146) = {41, 9};
Line(147) = {42, 10};
Line(148) = {39, 7};
Line(149) = {52, 23};
Line(150) = {53, 24};
Line(151) = {54, 25};
Line(152) = {51, 22};
Line(153) = {44, 13};
Line(154) = {45, 14};
Line(155) = {46, 15};
Line(156) = {43, 12};
Line(157) = {48, 18};
Line(158) = {49, 19};
Line(159) = {50, 20};
Line(160) = {47, 17};

// liste des surfaces
Line Loop(161) = {141, 2, -142, -105};
Plane Surface(162) = {161};
Line Loop(163) = {142, 3, -143, -106};
Plane Surface(164) = {163};
Line Loop(165) = {143, 4, -144, -107};
Plane Surface(166) = {165};
Line Loop(167) = {144, 1, -141, -108};
Plane Surface(168) = {167};
Line Loop(169) = {145, 6, -146, -109};
Plane Surface(170) = {169};
Line Loop(171) = {146, 7, -147, -110};
Plane Surface(172) = {171};
Line Loop(173) = {147, 8, -148, -111};
Plane Surface(174) = {173};
Line Loop(175) = {148, 5, -145, -112};
Plane Surface(176) = {175};
Line Loop(177) = {149, 18, -150, -121};
Plane Surface(178) = {177};
Line Loop(179) = {150, 19, -151, -122};
Plane Surface(180) = {179};
Line Loop(181) = {151, 20, -152, -123};
Plane Surface(182) = {181};
Line Loop(183) = {152, 17, -149, -124};
Plane Surface(184) = {183};
Line Loop(185) = {153, 10, -154, -113};
Plane Surface(186) = {185};
Line Loop(187) = {154, 11, -155, -114};
Plane Surface(188) = {187};
Line Loop(189) = {155, 12, -156, -115};
Plane Surface(190) = {189};
Line Loop(191) = {156, 9, -153, -116};
Plane Surface(192) = {191};
Line Loop(193) = {157, 14, -158, -117};
Plane Surface(194) = {193};
Line Loop(195) = {158, 15, -159, -118};
Plane Surface(196) = {195};
Line Loop(197) = {159, 16, -160, -119};
Plane Surface(198) = {197};
Line Loop(199) = {160, 13, -157, -120};
Plane Surface(200) = {199};
Line Loop(201) = {105, 106, 107, 108};
Plane Surface(202) = {201};
Line Loop(203) = {109, 110, 111, 112};
Plane Surface(204) = {203};
Line Loop(205) = {121, 122, 123, 124};
Plane Surface(206) = {205};
Line Loop(207) = {113, 114, 115, 116};
Plane Surface(208) = {207};
Line Loop(209) = {117, 118, 119, 120};
Plane Surface(210) = {209};
Line Loop(211) = {125, 145, 38, -141};
Plane Surface(212) = {211};
Line Loop(213) = {126, 146, 39, -142};
Plane Surface(214) = {213};
Line Loop(215) = {127, 147, 40, -143};
Plane Surface(216) = {215};
Line Loop(217) = {128, 148, 37, -144};
Plane Surface(218) = {217};
Line Loop(219) = {125, 109, -126, -105};
Plane Surface(220) = {219};
Line Loop(221) = {126, 110, -127, -106};
Plane Surface(222) = {221};
Line Loop(223) = {127, 111, -128, -107};
Plane Surface(224) = {223};
Line Loop(225) = {128, 112, -125, -108};
Plane Surface(226) = {225};
Line Loop(227) = {129, 157, 34, -153};
Plane Surface(228) = {227};
Line Loop(229) = {130, 158, 35, -154};
Plane Surface(230) = {229};
Line Loop(231) = {131, 159, 36, -155};
Plane Surface(232) = {231};
Line Loop(233) = {132, 160, 33, -156};
Plane Surface(234) = {233};
Line Loop(235) = {129, 117, -130, -113};
Plane Surface(236) = {235};
Line Loop(237) = {130, 118, -131, -114};
Plane Surface(238) = {237};
Line Loop(239) = {131, 119, -132, -115};
Plane Surface(240) = {239};
Line Loop(241) = {132, 120, -129, -116};
Plane Surface(242) = {241};
Line Loop(243) = {38, 2, -39, -6};
Ruled Surface(244) = {243};
Line Loop(245) = {39, 3, -40, -7};
Ruled Surface(246) = {245};
Line Loop(247) = {40, 4, -37, -8};
Ruled Surface(248) = {247};
Line Loop(249) = {37, 1, -38, -5};
Ruled Surface(250) = {249};
Line Loop(251) = {34, 10, -35, -14};
Ruled Surface(252) = {251};
Line Loop(253) = {35, 11, -36, -15};
Ruled Surface(254) = {253};
Line Loop(255) = {36, 12, -33, -16};
Ruled Surface(256) = {255};
Line Loop(257) = {33, 9, -34, -13};
Ruled Surface(258) = {257};
Line Loop(259) = {28, 6, -30, -18};
Ruled Surface(260) = {259};
Line Loop(261) = {30, 7, -32, -19};
Ruled Surface(262) = {261};
Line Loop(263) = {32, 8, -26, -20};
Ruled Surface(264) = {263};
Line Loop(265) = {26, 5, -28, -17};
Ruled Surface(266) = {265};
Line Loop(267) = {10, 29, -18, -27};
Ruled Surface(268) = {267};
Line Loop(269) = {11, 31, -19, -29};
Ruled Surface(270) = {269};
Line Loop(271) = {31, 20, -25, -12};
Ruled Surface(272) = {271};
Line Loop(273) = {25, 17, -27, -9};
Ruled Surface(274) = {273};
Line Loop(285) = {121, 138, -113, -137};
Plane Surface(286) = {285};
Line Loop(287) = {122, 139, -114, -138};
Plane Surface(288) = {287};
Line Loop(289) = {139, 115, -140, -123};
Plane Surface(290) = {289};
Line Loop(291) = {124, 137, -116, -140};
Plane Surface(292) = {291};
Line Loop(293) = {121, -134, -109, 133};
Plane Surface(294) = {293};
Line Loop(295) = {134, 122, -135, -110};
Plane Surface(296) = {295};
Line Loop(297) = {135, 123, -136, -111};
Plane Surface(298) = {297};
Line Loop(299) = {136, 124, -133, -112};
Plane Surface(300) = {299};
Line Loop(301) = {133, 149, 28, -145};
Plane Surface(302) = {301};
Line Loop(303) = {134, 150, 30, -146};
Plane Surface(304) = {303};
Line Loop(305) = {147, -32, -151, -135};
Plane Surface(306) = {305};
Line Loop(307) = {148, -26, -152, -136};
Plane Surface(308) = {307};
Line Loop(309) = {149, -27, -153, -137};
Plane Surface(310) = {309};
Line Loop(311) = {150, -29, -154, -138};
Plane Surface(312) = {311};
Line Loop(313) = {139, 155, 31, -151};
Plane Surface(314) = {313};
Line Loop(315) = {152, -25, -156, -140};
Plane Surface(316) = {315};

// listes des volumes
Surface Loop(317) = {210, 208, 236, 238, 240, 242};
Volume(318) = {317};
Surface Loop(319) = {194, 252, 186, 230, 236, 228};
Volume(320) = {319};
Surface Loop(321) = {196, 254, 188, 232, 238, 230};
Volume(322) = {321};
Surface Loop(323) = {232, 240, 234, 190, 198, 256};
Volume(324) = {323};
Surface Loop(325) = {200, 258, 192, 242, 228, 234};
Volume(326) = {325};
Surface Loop(327) = {292, 274, 310, 192, 184, 316};
Volume(328) = {327};
Surface Loop(329) = {286, 268, 312, 310, 186, 178};
Volume(330) = {329};
Surface Loop(331) = {270, 188, 180, 288, 314, 312};
Volume(332) = {331};
Surface Loop(333) = {290, 272, 314, 316, 182, 190};
Volume(334) = {333};
Surface Loop(335) = {286, 292, 290, 288, 206, 208};
Volume(336) = {335};
Surface Loop(337) = {262, 296, 304, 306, 180, 172};
Volume(338) = {337};
Surface Loop(339) = {174, 264, 306, 308, 298, 182};
Volume(340) = {339};
Surface Loop(341) = {176, 266, 300, 308, 302, 184};
Volume(342) = {341};
Surface Loop(343) = {170, 304, 302, 294, 178, 260};
Volume(344) = {343};
Surface Loop(345) = {298, 300, 296, 294, 206, 204};
Volume(346) = {345};
Surface Loop(347) = {202, 204, 222, 224, 220, 226};
Volume(348) = {347};
Surface Loop(349) = {244, 162, 170, 220, 214, 212};
Volume(350) = {349};
Surface Loop(351) = {168, 250, 176, 218, 226, 212};
Volume(352) = {351};
Surface Loop(353) = {172, 164, 246, 216, 222, 214};
Volume(354) = {353};
Surface Loop(355) = {174, 166, 248, 218, 224, 216};
Volume(356) = {355};


// DEFINITION DU MAILLAGE
//=======================

// nombre de noeuds par ligne
Transfinite Line {116, 120, 117, 113, 114, 118, 115, 119} = nc Using Progression 1;
Transfinite Line {12, 16, 13, 9, 14, 10, 15, 11} = nc Using Progression 1;
Transfinite Line {123, 20, 17, 124, 18, 121, 122, 19} = nc Using Progression 1;
Transfinite Line {111, 107, 108, 112, 109, 105, 106, 110} = nc Using Progression 1;
Transfinite Line {4, 8, 1, 5, 6, 2, 7, 3} = nc Using Progression 1;

Transfinite Line {160, 156, 157, 153, 158, 154, 159, 155} = nt Using Progression 1;
Transfinite Line {151, 152, 149, 150, 148, 144, 145, 141} = nt Using Progression 1;
Transfinite Line {146, 142, 143, 147} = nt Using Progression 1;

Transfinite Line {137, 27, 138, 29, 140, 25, 139, 31, 135} = nl Using Progression 1;
Transfinite Line {32, 136, 26, 30, 134, 133, 28} = nl Using Progression 1;

Transfinite Line {39, 126, 127, 40, 37, 128, 125, 38, 34} = np Using Progression 1;
Transfinite Line {35, 130, 129, 132, 33, 131, 36} = np Using Progression 1;

// maillage des surfaces
Transfinite Surface {162:274:2};
Transfinite Surface {286:316:2};
Recombine Surface {162:274:2};
Recombine Surface {286:316:2};

// maillage des volumes
Transfinite Volume{324} = {12, 43, 46, 15, 50, 20, 17, 47};
Transfinite Volume{320} = {48, 18, 19, 49, 44, 13, 14, 45};
Transfinite Volume{318} = {47, 48, 49, 50, 43, 44, 45, 46};
Transfinite Volume{336} = {45, 46, 43, 44, 53, 54, 51, 52};
Transfinite Volume{332} = {53, 24, 25, 54, 45, 14, 15, 46};
Transfinite Volume{334} = {54, 25, 22, 51, 46, 15, 12, 43};
Transfinite Volume{328} = {51, 22, 23, 52, 43, 12, 13, 44};
Transfinite Volume{330} = {23, 24, 53, 52, 13, 14, 45, 44};
Transfinite Volume{346} = {53, 54, 51, 52, 41, 42, 39, 40};
Transfinite Volume{340} = {42, 10, 7, 39, 54, 25, 22, 51};
Transfinite Volume{342} = {39, 7, 22, 51, 40, 8, 23, 52};
Transfinite Volume{344} = {40, 8, 23, 52, 41, 53, 24, 9};
Transfinite Volume{338} = {10, 42, 54, 25, 53, 24, 9, 41};
Transfinite Volume{354} = {9, 10, 42, 41, 4, 5, 38, 37};
Transfinite Volume{356} = {38, 5, 2, 35, 42, 10, 7, 39};
Transfinite Volume{352} = {35, 2, 3, 36, 39, 7, 8, 40};
Transfinite Volume{350} = {36, 3, 4, 37, 40, 8, 9, 41};
Transfinite Volume{348} = {41, 42, 39, 40, 37, 38, 35, 36};
Transfinite Volume{322} = {45, 14, 15, 46, 49, 19, 20, 50};
Transfinite Volume{326} = {47, 17, 18, 48, 43, 12, 13, 44};
Transfinite Volume{324} = {50, 20, 17, 47, 46, 15, 12, 43};
Transfinite Volume{338} = {9, 10, 42, 41, 24, 25, 54, 53};
Transfinite Volume{344} = {9, 41, 40, 8, 24, 53, 52, 23};



// DESSIN DES ARCS DE CERCLES DES COINS

Point(1000) = {R, 2*E+B, 0, 1};
Point(1001) = {0, 2*E+B, -R, 1};
Point(1002) = {-R, 2*E+B, 0, 1};
Point(1003) = {0, 2*E+B, R, 1};

Point(1004) = {1.5*R, 2*E+B, -R, 1};
Point(1005) = {-R, 2*E+B, 1.5*R, 1};

Point(1006) = {2*R, 2*E+B, R, 1};
Point(1007) = {R, 2*E+B, 2*R, 1};
Point(1008) = {2*R, 2*E+B, 2*R, 1};
Point(1009) = {-R, 2*E+B, 2*R, 1};
Point(1010) = {2*R, 2*E+B, -R, 1};


Circle(2145) = {1001, 16, 1002};
Circle(2161) = {1006, 1008, 1007};
Line(2146) = {18, 1001};
Line(2147) = {19, 1002};
Line(2149) = {1004, 1001};
Line(2150) = {1004, 17};
Line(2156) = {1005, 1002};
Line(2157) = {1005, 20};
Line(2162) = {1005, 1009};
Line(2163) = {1009, 1007};
Line(2164) = {1006, 1010};
Line(2165) = {1010, 1004};
Line(2167) = {1006, 17};
Line(2168) = {1007, 20};

Line Loop(2154) = {2146, 2145, -2147, -14};
Plane Surface(2155) = {2154};
Line Loop(2169) = {2163, 2168, -2157, 2162};
Plane Surface(2170) = {2169};
Line Loop(2171) = {2157, -15, 2147, -2156};
Plane Surface(2172) = {2171};
Line Loop(2173) = {2167, -2150, -2165, -2164};
Plane Surface(2174) = {2173};
Line Loop(2175) = {2161, 2168, 16, -2167};
Plane Surface(2176) = {2175};
Line Loop(2177) = {2150, 13, 2146, -2149};
Plane Surface(2178) = {2177};

Transfinite Line {2145,2156,2149,2161} = nc Using Progression 1;
Transfinite Line {2146,2147,2157,2163,2150,2164} = nh Using Progression 1;
Transfinite Line {2162,2168,2167,2165} = ng Using Progression 1;

Extrude {0, ep, 0} {
  Surface{2170, 2172, 2155, 2178, 2174, 2176, 196, 194, 200, 198, 210};
}

Transfinite Line {2225,2205,2249,2290,2203,2315,2227,2337,2247,2359,2292,2381} = nc Using Progression 1;
Transfinite Line {2224,2204,2182,2180,2246,2271} = nh Using Progression 1;
Transfinite Line {2183,2181,2268,2270} = ng Using Progression 1;
Transfinite Line {2185,2194,2216,2230,2251,2282,2273,2186,2190,2212,2229,2252,2326,2317,2339,2361} = nep Using Progression 1;

Transfinite Line {2312,2334,2356,2314} = nt Using Progression 1;

Transfinite Surface {2275,2253,2155,2170,2172,2174,2176,2178,2297,2287,2283,2265,2235,2221,2187,2191,2195,2217,2243,2213,2327,2331,2353,2420,2319,2305,2222,2266,2310,2288,2200,2199,2231,2397,2244,2257,2375,2398,2363,2327,2341,2332,2376,2354};
Recombine Surface {2275,2253,2155,2170,2172,2174,2176,2178,2297,2287,2283,2265,2235,2221,2187,2191,2195,2217,2243,2213,2327,2331,2353,2420,2319,2305,2222,2266,2310,2288,2200,2199,2231,2397,2244,2257,2375,2398,2363,2327,2341,2332,2376,2354};

// suite pour faire des parties verticales

Point(2007) = {R, 2*E+B, 2*R+ep2, 1};
Point(2009) = {-R, 2*E+B, 2*R+ep2, 1};
Point(2006) = {2*R+ep2, 2*E+B, R, 1};
Point(2010) = {2*R+ep2, 2*E+B, -R, 1};
Point(2107) = {R, 2*E+B+ep, 2*R+ep2, 1};
Point(2109) = {-R, 2*E+B+ep, 2*R+ep2, 1};
Point(2106) = {2*R+ep2, 2*E+B+ep, R, 1};
Point(2110) = {2*R+ep2, 2*E+B+ep, -R, 1};


Line(2421) = {2009, 2109};
Line(2422) = {2109, 1011};
Line(2423) = {2009, 1009};
Line(2424) = {2009, 2007};
Line(2425) = {2007, 1007};
Line(2426) = {2107, 1012};
Line(2427) = {2107, 2109};
Line(2428) = {2107, 2007};
Line(2429) = {1036, 2106};
Line(2430) = {1006, 2006};
Line(2431) = {2006, 2106};
Line(2432) = {2106, 2110};
Line(2433) = {2110, 1045};
Line(2434) = {2010, 2110};
Line(2435) = {2010, 1010};
Line(2436) = {2010, 2006};

Transfinite Line {2424,2432,2436,2427} = nh Using Progression 1;
Transfinite Line {2421,2428,2431,2434} = nep Using Progression 1;
Transfinite Line {2425,2426,2423,2422,2430,2429,2435,2433} = nep2 Using Progression 1;

Line Loop(2437) = {2421, 2422, -2185, -2423};
Plane Surface(2438) = {2437};
Line Loop(2439) = {2423, 2163, -2425, -2424};
Plane Surface(2440) = {2439};
Line Loop(2441) = {2425, 2186, -2426, 2428};
Plane Surface(2442) = {2441};
Line Loop(2443) = {2426, -2180, -2422, -2427};
Plane Surface(2444) = {2443};
Line Loop(2445) = {2424, -2428, 2427, -2421};
Plane Surface(2446) = {2445};
Line Loop(2447) = {2273, 2429, -2431, -2430};
Plane Surface(2448) = {2447};
Line Loop(2449) = {2430, -2436, 2435, -2164};
Plane Surface(2450) = {2449};
Line Loop(2451) = {2434, 2433, -2282, -2435};
Plane Surface(2452) = {2451};
Line Loop(2453) = {2432, 2433, 2271, 2429};
Plane Surface(2454) = {2453};
Line Loop(2455) = {2431, 2432, -2434, 2436};
Plane Surface(2456) = {2455};

Transfinite Surface {2438,2440,2442,2444,2446,2448,2450,2452,2454,2456};
Recombine Surface {2438,2440,2442,2444,2446,2448,2450,2452,2454,2456};

Surface Loop(2457) = {2438, 2446, 2440, 2442, 2444, 2187};
Volume(2458) = {2457};
Surface Loop(2459) = {2448, 2454, 2456, 2452, 2450, 2287};
Volume(2460) = {2459};

Transfinite Volume{2458,2460,357,358,363,359,365,360,361,362,366,364,367};

// MONTANTS VERTICAUX

Extrude {0, HH, 0} {
  Surface{2444, 2454};
}
Transfinite Line {2464,2462,2487,2485} = nep2 Using Progression 1;
Transfinite Line {2490,2494,2489,2498,2468,2467,2472,2476} = n3 Using Progression 1;
Transfinite Line {2484,2486,2463,2465} = nh Using Progression 1;

Transfinite Surface {2491,2499,2494,2498,2486,2467,2473,2481,2477,2465,2503,2495,2504,2482,2469};
Recombine Surface {2491,2499,2494,2498,2486,2467,2473,2481,2477,2465,2503,2495,2504,2482,2469};

Transfinite Volume{2462,2461};

Extrude {0, ep3, 0} {
  Surface{2482, 2504};
}

Transfinite Line {2509,2507,2530,2528} = nh Using Progression 1;
Transfinite Line {2511,2512,2520,2516,2542,2533,2538,2534} = n4 Using Progression 1;
Transfinite Line {2529,2531,2506,2508} = nep2 Using Progression 1;

Transfinite Surface {2521,2517,2513,2525,2526,2547,2535,2548,2543,2539,2513};
Recombine Surface {2521,2517,2513,2525,2526,2547,2535,2548,2543,2539,2513};

Transfinite Volume{2463,2464};

// COIN DU HAUT

Line(2549) = {2130, 2112};
Line(2550) = {2126, 2116};
Line Loop(2551) = {2486, 2549, 2463, -2550};
Plane Surface(2552) = {2551};
Line(2553) = {2146, 2136};
Line(2554) = {2150, 2132};

Transfinite Line {2549,2550,2553,2554} = n5 Using Progression 1;

Line Loop(2555) = {2530, 2554, 2507, -2553};
Plane Surface(2556) = {2555};
Line Loop(2557) = {2542, 2554, -2512, -2549};
Plane Surface(2558) = {2557};
Line Loop(2559) = {2538, 2553, -2516, -2550};
Plane Surface(2560) = {2559};

Surface Loop(2561) = {2560, 2556, 2558, 2552, 2517, 2543};
Volume(2562) = {2561};

Transfinite Surface {2556,2558,2560,2552};
Recombine Surface {2556,2558,2560,2552};
Transfinite Volume{2562};

// SYMMETRIES POUR LES 3 AUTRES COINS
Symmetry {-1, 0, 0, L1/2} {
  Duplicata { Volume{2463, 2464, 2562, 2461, 2462, 2458, 362, 357, 2460, 361, 358, 366, 324, 363, 367, 322, 360, 334, 365, 318, 326, 332, 364, 336, 320, 328, 359, 340, 330, 338, 356, 346, 342, 354, 348, 344, 352, 350}; }
}
Symmetry {0, 0, -1, L2/2} {
  Duplicata { Volume{2563, 2625, 2463, 2594, 2656, 2464, 2562, 2461, 2687, 2462, 2718, 2780, 2873, 2749, 2966, 2904, 3028, 2997, 2935, 3245, 3152, 3214, 3369, 3121, 3307, 3090, 2842, 3183, 3276, 2458, 3431, 3059, 2811, 3462, 3338, 3400, 362, 3524, 3586, 357, 3648, 3493, 3555, 3617, 2460, 3710, 3679, 361, 358, 366, 363, 324, 367, 322, 365, 360, 318, 334, 364, 326, 332, 320, 336, 359, 328, 340, 330, 338, 346, 356, 342, 354, 344, 348, 352, 350}; }
}


Line(6049) = {8935, 10143};
Line(6050) = {8903, 10111};
Line(6051) = {6511, 6607};
Line(6052) = {6510, 6606};
Line(6053) = {6031, 6223};
Line(6054) = {6047, 6239};
Line(6055) = {6030, 6222};
Line(6056) = {6046, 6238};
Line(6057) = {2936, 2006};
Line(6058) = {2535, 2106};
Line(6059) = {2536, 2110};
Line(6060) = {2968, 2010};
Line(6061) = {2247, 2121};
Line(6062) = {2263, 2141};
Line(6063) = {2264, 2142};
Line(6064) = {2248, 2122};
Line(6065) = {2631, 6702};
Line(6066) = {2448, 6135};
Line(6067) = {2648, 6719};
Line(6068) = {2439, 6126};
Line(6069) = {2160, 5751};
Line(6070) = {2176, 5767};
Line(6071) = {2151, 5742};
Line(6072) = {2167, 5758};
Line(6073) = {2111, 5934};
Line(6074) = {2131, 5950};
Line(6075) = {2120, 5943};
Line(6076) = {2140, 5959};
Line(6077) = {8594, 2009};
Line(6078) = {8611, 2007};
Line(6079) = {6423, 2109};
Line(6080) = {6414, 2107};

Transfinite Line {6049,6050,6051,6052,6053,6054,6055,6056,6057,6058,6059,6060,6061} = n6 Using Progression 1;
Transfinite Line {6062,6063,6064,6065,6066,6067,6068,6069,6070,6071,6072,6073,6074} = n6 Using Progression 1;
Transfinite Line {6075,6076,6077,6078,6079,6080} = n6 Using Progression 1;


Line Loop(6081) = {6077, 2421, -6079, -4618};
Plane Surface(6082) = {6081};
Line Loop(6083) = {2424, -6078, -4623, 6077};
Plane Surface(6084) = {6083};
Line Loop(6085) = {6078, -2428, -6080, -4624};
Plane Surface(6086) = {6085};
Line Loop(6087) = {6080, 2427, -6079, 3939};
Plane Surface(6088) = {6087};
Surface Loop(6089) = {6088, 6086, 6084, 6082, 2446, 4622};
Volume(6090) = {6089};

Line Loop(6091) = {6075, 3802, -6076, -2520};
Plane Surface(6092) = {6091};
Line Loop(6093) = {6075, 3784, -6073, -2465};
Plane Surface(6094) = {6093};
Line Loop(6095) = {6073, -3794, -6074, -2511};
Plane Surface(6096) = {6095};
Line Loop(6097) = {6074, -3789, -6076, 2509};
Plane Surface(6098) = {6097};
Surface Loop(6099) = {6092, 6094, 6096, 6098, 3805, 2525};
Volume(6100) = {6099};

Line Loop(6101) = {4724, 6050, -5096, -6049};
Plane Surface(6102) = {6101};
Line Loop(6103) = {6050, -5085, -6052, 4713};
Plane Surface(6104) = {6103};
Line Loop(6105) = {6052, 3998, -6051, -3967};
Plane Surface(6106) = {6105};
Line Loop(6107) = {6051, 5095, -6049, -4723};
Plane Surface(6108) = {6107};
Surface Loop(6109) = {6108, 6106, 6104, 6102, 5092, 4720};
Volume(6110) = {6109};

Line Loop(6111) = {6053, 3885, -6054, -3823};
Plane Surface(6112) = {6111};
Line Loop(6113) = {6053, -3874, -6055, 3812};
Plane Surface(6114) = {6113};
Line Loop(6115) = {6055, -3887, -6056, 3825};
Plane Surface(6116) = {6115};
Line Loop(6117) = {6056, 3879, -6054, -3817};
Plane Surface(6118) = {6117};
Surface Loop(6119) = {6112, 6114, 6116, 6118, 3821, 3883};
Volume(6120) = {6119};

Line Loop(6121) = {2725, 6067, -4034, -6065};
Plane Surface(6122) = {6121};
Line Loop(6123) = {6067, 4035, -6068, -2726};
Plane Surface(6124) = {6123};
Line Loop(6125) = {6068, -3846, -6066, 2661};
Plane Surface(6126) = {6125};
Line Loop(6127) = {6066, -4029, -6065, 2720};
Plane Surface(6128) = {6127};
Surface Loop(6129) = {6124, 6122, 6128, 6126, 2724, 4033};
Volume(6130) = {6129};

Line Loop(6131) = {6071, -3722, -6069, 2568};
Plane Surface(6132) = {6131};
Line Loop(6133) = {2578, 6071, -3732, -6072};
Plane Surface(6134) = {6133};
Line Loop(6135) = {6072, -3727, -6070, 2573};
Plane Surface(6136) = {6135};
Line Loop(6137) = {6070, -3740, -6069, 2586};
Plane Surface(6138) = {6137};
Surface Loop(6139) = {6138, 6136, 6134, 6132, 3743, 2589};
Volume(6140) = {6139};

Line Loop(6141) = {6057, -2436, -6060, 2826};
Plane Surface(6142) = {6141};
Line Loop(6143) = {6058, -2431, -6057, -2815};
Plane Surface(6144) = {6143};
Line Loop(6145) = {6058, 2432, -6059, -2689};
Plane Surface(6146) = {6145};
Line Loop(6147) = {6059, -2434, -6060, -2825};
Plane Surface(6148) = {6147};
Surface Loop(6149) = {6148, 6146, 6144, 6142, 2822, 2456};
Volume(6150) = {6149};

Line Loop(6151) = {6061, 2484, -6064, -2596};
Plane Surface(6152) = {6151};
Line Loop(6153) = {6062, -2533, -6061, -2609};
Plane Surface(6154) = {6153};
Line Loop(6155) = {6062, 2528, -6063, -2601};
Plane Surface(6156) = {6155};
Line Loop(6157) = {6063, -2534, -6064, 2607};
Plane Surface(6158) = {6157};
Surface Loop(6159) = {6152, 6154, 6156, 6158, 2535, 2605};
Volume(6160) = {6159};

Transfinite Surface {6082,6084,6086,6088,6092,6094,6096,6098,6102,6104,6106,6108,6112,6114,6116,6118,6122,6124,6126,6128,6132,6134,6136,6138,6142,6144,6146,6148,6152,6154,6156,6158};
Recombine Surface {6082,6084,6086,6088,6092,6094,6096,6098,6102,6104,6106,6108,6112,6114,6116,6118,6122,6124,6126,6128,6132,6134,6136,6138,6142,6144,6146,6148,6152,6154,6156,6158};

Transfinite Volume{6090,6100,6110,6120,6130,6140,6150,6160};


// DEFINITION DES GROUPES PHYSIQUES
//=================================

// face du haut
Physical Surface(2) = {194, 200, 198, 210, 196, 2905, 2998, 2967, 3246, 3122, 4245, 4152, 4183, 4431, 4307, 5330, 5516, 5268, 5237, 5392};

// face du bas
//Physical Surface(1) = {164, 166, 168, 162, 202, 3680, 3499, 3618, 3716, 3592, 4901, 4994, 5144, 5118, 5051, 5862, 5924, 6048, 6012, 5981};
// pied 2: Physical Surface(6161) = {3680, 3499, 3618, 3716, 3592};
// pied 3: Physical Surface(6162) = {4901, 4994, 5144, 5118, 5051};
// pied 4: Physical Surface(6163) = {5862, 5924, 6048, 6012, 5981};

Physical Surface(1) = {164, 166, 168, 162, 202};
Physical Surface(3) = {3680, 3499, 3618, 3716, 3592, 4901, 4994, 5144, 5118, 5051, 5862, 5924, 6048, 6012, 5981};

// volume des plaque en acier dessus et dessous
Physical Volume(10) = {326, 320, 318, 324, 322, 352, 350, 348, 356, 354, 5918, 6042, 5980, 5856, 6011, 5608, 5546, 5298, 5360, 5453,4988, 5143, 5112, 4895, 5050, 4275, 4213, 4461, 4554, 4337, 3679, 3493, 3617, 3710, 3586, 2935, 3183, 3152, 3307, 3028};

// volume en caoutchouc
Physical Volume(11) = {328, 330, 336, 334, 332, 342, 344, 346, 340, 338, 5794, 5949, 5887, 5732, 5825, 5577, 5763, 5639, 5484, 5701, 4802, 5019, 4957, 4740, 4864, 4492, 4771, 4585, 4647, 4368, 3555, 3400, 3524, 3648, 3462, 3214, 3090, 3276, 3431, 3338};

// volume structure au dessus
Physical Volume(12) = {3903, 3779, 3934, 6100, 3872, 3996, 6120, 2463, 2461, 2562, 3810, 3965, 2464, 3748, 2462, 3717, 6160, 3841, 2594, 6140, 2687, 2563, 2656, 2625, 6090, 6130, 4616, 4926, 4833, 5205, 5422, 5174, 5081, 6110, 4120, 4058, 4027, 5267, 5329, 5236, 5515, 5670, 5391, 4709, 4523, 4182, 4151, 4244, 4678, 4430, 4306, 2718, 2780, 2873, 2966, 2749, 2904, 2997, 3245, 3369, 3121, 2842, 3059, 2811, 2458, 6150, 362, 357, 361, 358, 363, 367, 365, 364, 359, 360, 366, 2460, 4089, 4399};


