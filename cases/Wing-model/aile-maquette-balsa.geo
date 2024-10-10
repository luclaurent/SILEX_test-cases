// Taille elements (m)
h = 3.0e-3  ;
corde = 0.170 ;
d = 0.05;

// Par defaut, triangles a 3 noeuds
// Pour des triangles a 6 noeuds, activer:
// Mesh.ElementOrder = 2;

Point(2 ) = {corde*0.0000,0.0,corde*0.0005 , h };
Point(3 ) = {corde*0.0011,0.0,corde*0.0062 , h };
Point(4 ) = {corde*0.0043,0.0,corde*0.0125 , h };
Point(5 ) = {corde*0.0096,0.0,corde*0.0194 , h };
Point(6 ) = {corde*0.0170,0.0,corde*0.0265 , h };
Point(7 ) = {corde*0.0265,0.0,corde*0.0335 , h };
Point(8 ) = {corde*0.0381,0.0,corde*0.0403 , h };
Point(9 ) = {corde*0.0516,0.0,corde*0.0468 , h };
Point(10 ) = {corde*0.0670,0.0,corde*0.0531 , h };
Point(11 ) = {corde*0.0843,0.0,corde*0.0594 , h };
Point(12 ) = {corde*0.1033,0.0,corde*0.0655 , h };
Point(13 ) = {corde*0.1241,0.0,corde*0.0713 , h };
Point(14 ) = {corde*0.1465,0.0,corde*0.0766 , h };
Point(15 ) = {corde*0.1703,0.0,corde*0.0811 , h };
Point(16 ) = {corde*0.1956,0.0,corde*0.0848 , h };
Point(17 ) = {corde*0.2222,0.0,corde*0.0877 , h };
Point(18 ) = {corde*0.2500,0.0,corde*0.0900 , h };
Point(19 ) = {corde*0.2789,0.0,corde*0.0916 , h };
Point(20 ) = {corde*0.3087,0.0,corde*0.0927 , h };
Point(21 ) = {corde*0.3393,0.0,corde*0.0932 , h };
Point(22 ) = {corde*0.3706,0.0,corde*0.0931 , h };
Point(23 ) = {corde*0.4347,0.0,corde*0.0913 , h };
Point(24 ) = {corde*0.5000,0.0,corde*0.0872 , h };
Point(25 ) = {corde*0.5653,0.0,corde*0.0811 , h };
Point(26 ) = {corde*0.6294,0.0,corde*0.0732 , h };
Point(27 ) = {corde*0.6913,0.0,corde*0.0641 , h };
Point(28 ) = {corde*0.7500,0.0,corde*0.0541 , h };
Point(29 ) = {corde*0.8044,0.0,corde*0.0439 , h };
Point(30 ) = {corde*0.8536,0.0,corde*0.0340 , h };
Point(31 ) = {corde*0.8967,0.0,corde*0.0248 , h };
Point(32 ) = {corde*0.9330,0.0,corde*0.0166 , h };
Point(33 ) = {corde*0.9619,0.0,corde*0.0097 , h };
Point(34 ) = {corde*0.9830,0.0,corde*0.0045 , h };
Point(35 ) = {corde*0.9957,0.0,corde*0.0012 , h };
Point(36 ) = {corde*1.0000,0.0,corde*0.0000 , h };
Point(37 ) = {corde*0.0000,0.0,corde*0.0005 , h };
Point(38 ) = {corde*0.0011,0.0,corde*-0.0045 , h };
Point(39 ) = {corde*0.0043,0.0,corde*-0.0090 , h };
Point(40 ) = {corde*0.0096,0.0,corde*-0.0130 , h };
Point(41 ) = {corde*0.0170,0.0,corde*-0.0165 , h };
Point(42 ) = {corde*0.0265,0.0,corde*-0.0196 , h };
Point(43 ) = {corde*0.0381,0.0,corde*-0.0221 , h };
Point(44 ) = {corde*0.0516,0.0,corde*-0.0241 , h };
Point(45 ) = {corde*0.0670,0.0,corde*-0.0257 , h };
Point(46 ) = {corde*0.0843,0.0,corde*-0.0268 , h };
Point(47 ) = {corde*0.1033,0.0,corde*-0.0276 , h };
Point(48 ) = {corde*0.1241,0.0,corde*-0.0282 , h };
Point(49 ) = {corde*0.1465,0.0,corde*-0.0284 , h };
Point(50 ) = {corde*0.1703,0.0,corde*-0.0283 , h };
Point(51 ) = {corde*0.1956,0.0,corde*-0.0280 , h };
Point(52 ) = {corde*0.2222,0.0,corde*-0.0273 , h };
Point(53 ) = {corde*0.2500,0.0,corde*-0.0265 , h };
Point(54 ) = {corde*0.2789,0.0,corde*-0.0256 , h };
Point(55 ) = {corde*0.3087,0.0,corde*-0.0246 , h };
Point(56 ) = {corde*0.3393,0.0,corde*-0.0235 , h };
Point(57 ) = {corde*0.3706,0.0,corde*-0.0224 , h };
Point(58 ) = {corde*0.4347,0.0,corde*-0.0202 , h };
Point(59 ) = {corde*0.5000,0.0,corde*-0.0179 , h };
Point(60 ) = {corde*0.5653,0.0,corde*-0.0157 , h };
Point(61 ) = {corde*0.6294,0.0,corde*-0.0135 , h };
Point(62 ) = {corde*0.6913,0.0,corde*-0.0113 , h };
Point(63 ) = {corde*0.7500,0.0,corde*-0.0093 , h };
Point(64 ) = {corde*0.8044,0.0,corde*-0.0074 , h };
Point(65 ) = {corde*0.8536,0.0,corde*-0.0058 , h };
Point(66 ) = {corde*0.8967,0.0,corde*-0.0043 , h };
Point(67 ) = {corde*0.9330,0.0,corde*-0.0030 , h };
Point(68 ) = {corde*0.9619,0.0,corde*-0.0019 , h };
Point(69 ) = {corde*0.9830,0.0,corde*-0.0009 , h };
Point(70 ) = {corde*0.9957,0.0,corde*-0.0003 , h };
Point(71 ) = {corde*1.0000,0.0,corde*0.0000 , h };

// bord de fuite
Point(72 ) = {corde*(0.5684*0.8967+0.4315*0.9330),0.0,corde*(-0.5684*0.0043-0.4315*0.0030) , h };
Point(73 ) = {corde*(0.5684*0.8967+0.4315*0.9330),0.0,corde*(0.5684*0.0248+0.4315*0.0166) , h };

// bord d'attaque
Point(74 ) = {corde*(0.307*0.0170+0.693*0.0265),0.0,corde*(-0.307*0.0165-0.693*0.0196) , h };
Point(75 ) = {corde*(0.307*0.0170+0.693*0.0265),0.0,corde*(0.307*0.0265+0.693*0.0335) , h };

// longeron 0.059=distance au bord d'attaque
Point(76 ) = {corde*(0.4598*0.2789+0.5402*0.3087),0.0,corde*(-0.4598*0.0256-0.5402*0.0246) , h };
Point(77 ) = {corde*(0.4598*0.2789+0.5402*0.3087),0.0,corde*(0.4598*0.0916+0.5402*0.0927) , h };

Line(1) = {72, 73};
Line(2) = {74, 75};
Line(3) = {76, 77};
Line(4) = {72, 66};
Line(5) = {66, 65};
Line(6) = {65, 64};
Line(7) = {64, 63};
Line(8) = {63, 62};
Line(9) = {62, 61};
Line(10) = {61, 60};
Line(11) = {60, 59};
Line(12) = {59, 58};
Line(13) = {58, 57};
Line(14) = {57, 56};
Line(15) = {56, 55};
Line(16) = {55, 76};
Line(17) = {76, 54};
Line(18) = {54, 53};
Line(19) = {53, 52};
Line(20) = {52, 51};
Line(21) = {51, 50};
Line(22) = {50, 49};
Line(23) = {49, 48};
Line(24) = {48, 47};
Line(25) = {47, 46};
Line(26) = {46, 45};
Line(27) = {45, 44};
Line(28) = {44, 43};
Line(29) = {43, 42};
Line(30) = {42, 74};
Line(31) = {75, 7};
Line(32) = {7, 8};
Line(33) = {8, 9};
Line(34) = {9, 10};
Line(35) = {10, 11};
Line(36) = {11, 12};
Line(37) = {12, 13};
Line(38) = {13, 14};
Line(39) = {14, 15};
Line(40) = {15, 16};
Line(41) = {16, 17};
Line(42) = {17, 18};
Line(43) = {18, 19};
Line(44) = {19, 77};
Line(45) = {77, 20};
Line(46) = {20, 21};
Line(47) = {21, 22};
Line(48) = {22, 23};
Line(49) = {23, 24};
Line(50) = {24, 25};
Line(51) = {25, 26};
Line(52) = {26, 27};
Line(53) = {27, 28};
Line(54) = {28, 29};
Line(55) = {29, 30};
Line(56) = {30, 31};
Line(57) = {31, 73};



Line Loop(58) = {53, 54, 55, 56, 57, -1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3, 45, 46, 47, 48, 49, 50, 51, 52};
Plane Surface(59) = {58};
Line Loop(60) = {3, -44, -43, -42, -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -2, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17};
Plane Surface(61) = {60};


Translate {0, d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 2*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 3*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 4*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 5*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 6*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 7*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 8*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 9*d, 0} {
  Duplicata { Surface{59, 61}; }
}

Translate {0, 10*d, 0} {
  Duplicata { Surface{59, 61}; }
}

// longeron bord de fuite
Line(662) = {73, 95};
Line(663) = {72, 99};
Line(664) = {95, 321};
Line(665) = {99, 325};
Line(666) = {321, 547};
Line(667) = {325, 551};
Line(668) = {547, 773};
Line(669) = {551, 777};
Line(670) = {773, 999};
Line(671) = {777, 1003};
Line(672) = {999, 1225};
Line(673) = {1003, 1229};
Line(674) = {1225, 1451};
Line(675) = {1229, 1455};
Line(676) = {1451, 1677};
Line(677) = {1455, 1681};
Line(678) = {1677, 1903};
Line(679) = {1681, 1907};
Line(680) = {1903, 2129};
Line(681) = {1907, 2133};
Line Loop(722) = {1, 662, 68, -663};
Plane Surface(723) = {722};
Line Loop(724) = {68, 665, -128, -664};
Plane Surface(725) = {724};
Line Loop(726) = {128, 667, -188, -666};
Plane Surface(727) = {726};
Line Loop(728) = {188, 669, -248, -668};
Plane Surface(729) = {728};
Line Loop(730) = {248, 671, -308, -670};
Plane Surface(731) = {730};
Line Loop(732) = {308, 673, -368, -672};
Plane Surface(733) = {732};
Line Loop(734) = {368, 675, -428, -674};
Plane Surface(735) = {734};
Line Loop(736) = {428, 677, -488, -676};
Plane Surface(737) = {736};
Line Loop(738) = {488, 679, -548, -678};
Plane Surface(739) = {738};
Line Loop(740) = {548, 681, -608, -680};
Plane Surface(741) = {740};

// longeron central
Line(682) = {77, 155};
Line(683) = {155, 381};
Line(684) = {381, 607};
Line(685) = {607, 833};
Line(686) = {833, 1059};
Line(687) = {1059, 1285};
Line(688) = {1285, 1511};
Line(689) = {1511, 1737};
Line(690) = {1737, 1963};
Line(691) = {1963, 2189};
Line(692) = {2185, 1959};
Line(693) = {1959, 1733};
Line(694) = {1507, 1733};
Line(695) = {1507, 1281};
Line(696) = {1281, 1055};
Line(697) = {1055, 829};
Line(698) = {829, 603};
Line(699) = {603, 377};
Line(700) = {377, 151};
Line(701) = {151, 76};
Line Loop(742) = {3, 682, -82, 701};
Plane Surface(743) = {742};
Line Loop(744) = {82, 683, -142, 700};
Plane Surface(745) = {744};
Line Loop(746) = {142, 684, -202, 699};
Plane Surface(747) = {746};
Line Loop(748) = {202, 685, -262, 698};
Plane Surface(749) = {748};
Line Loop(750) = {262, 686, -322, 697};
Plane Surface(751) = {750};
Line Loop(752) = {322, 687, -382, 696};
Plane Surface(753) = {752};
Line Loop(754) = {382, 688, -442, 695};
Plane Surface(755) = {754};
Line Loop(756) = {442, 689, -502, -694};
Plane Surface(757) = {756};
Line Loop(758) = {502, 690, -562, 693};
Plane Surface(759) = {758};
Line Loop(760) = {562, 691, -622, 692};
Plane Surface(761) = {760};


// longeron bord d'attaque 
Line(702) = {75, 247};
Line(703) = {74, 251};
Line(704) = {251, 477};
Line(705) = {473, 247};
Line(706) = {473, 699};
Line(707) = {477, 703};
Line(708) = {699, 925};
Line(709) = {703, 929};
Line(710) = {925, 1151};
Line(711) = {929, 1155};
Line(712) = {1151, 1377};
Line(713) = {1155, 1381};
Line(714) = {1377, 1603};
Line(715) = {1381, 1607};
Line(716) = {1603, 1829};
Line(717) = {1607, 1833};
Line(718) = {1829, 2055};
Line(719) = {1833, 2059};
Line(720) = {2055, 2281};
Line(721) = {2059, 2285};
Line Loop(762) = {2, 702, 107, -703};
Plane Surface(763) = {762};
Line Loop(764) = {107, 704, -167, 705};
Plane Surface(765) = {764};
Line Loop(766) = {167, 707, -227, -706};
Plane Surface(767) = {766};
Line Loop(768) = {227, 709, -287, -708};
Plane Surface(769) = {768};
Line Loop(770) = {287, 711, -347, -710};
Plane Surface(771) = {770};
Line Loop(772) = {347, 713, -407, -712};
Plane Surface(773) = {772};
Line Loop(774) = {407, 715, -467, -714};
Plane Surface(775) = {774};
Line Loop(776) = {467, 717, -527, -716};
Plane Surface(777) = {776};
Line Loop(778) = {527, 719, -587, -718};
Plane Surface(779) = {778};
Line Loop(780) = {587, 721, -647, -720};
Plane Surface(781) = {780};



// longeron bord de fuite
Physical Surface(1) = {723, 725, 727, 729, 731, 733, 735, 737, 739, 741};

// longeron central
Physical Surface(2) = {743, 745, 747, 749, 751, 753, 755, 757, 759, 761};

// longeron bord d'attaque
Physical Surface(3) = {763, 765, 767, 769, 771, 773, 775, 777, 779, 781};

// nervures sauf emplanture
Physical Surface(4) = {62, 91, 122, 151, 182, 211, 242, 271, 302, 331, 362, 391, 422, 451, 482, 511, 542, 571, 602, 631};

// nervure emplanture
Physical Surface(5) = {59, 61};

