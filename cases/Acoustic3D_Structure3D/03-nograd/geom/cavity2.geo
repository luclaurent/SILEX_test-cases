// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.8;

ly2 = 3.0;
lx2 = 3.5;

lx3 = 1.5;
ly3 = 2.0;
l4 = 3.0;
lz4 = 2.0;
r4   = 0.5;
e4   = 0.2;
angle = Pi/2;
h = lx1/30;

lx5 = 6.0;
ly5 = 1.0;
lz5 = 1.5;
h5  = 0.5;

Mesh.CharacteristicLengthMax=h;
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



// wall: corners
Point(15) = {lx3,     ly3-e4/2  , 0, h};
Point(16) = {lx3+l4,     ly3-e4/2  , 0, h};
Point(17) = {lx3,     ly3-e4/2  , lz4-r4 , h};
Point(18) = {lx3+l4,     ly3-e4/2  , lz4-r4, h};
Point(19) = {lx3+r4,     ly3-e4/2  , lz4-r4 , h};
Point(20) = {lx3+l4-r4,     ly3-e4/2  , lz4-r4, h};
Point(21) = {lx3+r4,     ly3-e4/2  , lz4 , h};
Point(22) = {lx3+l4-r4,     ly3-e4/2  , lz4, h};

Point(23) = {lx3,     ly3+e4/2  , 0, h};
Point(24) = {lx3+l4,     ly3+e4/2  , 0, h};
Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h};
Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h};
Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h};
Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h};
Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , h};
Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, h};

// rotate wall
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
  Point{22, 30, 18, 26, 20, 21, 28, 29, 19, 27, 17, 25, 16, 24, 15, 23};
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

Line(34) = {16, 24};
Line(35) = {18, 26};
Line(36) = {22, 30};
Line(37) = {21, 29};
Line(38) = {17, 25};
Line(39) = {15, 23};


// surfaces, cavity
Line Loop(40) = {1, 2, 3, 4, 5, 6, 7};
Line Loop(41) = {22, 34, -28, -39};
Plane Surface(42) = {40, 41};
Line Loop(43) = {4, 12, -18, -11};
Plane Surface(44) = {43};
Line Loop(45) = {11, -17, -10, 3};
Plane Surface(46) = {45};
Line Loop(47) = {10, -16, -9, 2};
Plane Surface(48) = {47};
Line Loop(49) = {9, -15, -8, 1};
Plane Surface(50) = {49};
Line Loop(51) = {8, -21, -14, 7};
Plane Surface(52) = {51};
Line Loop(53) = {14, -20, -13, 6};
Plane Surface(54) = {53};
Line Loop(55) = {13, -19, -12, 5};
Plane Surface(56) = {55};
Line Loop(57) = {20, 21, 15, 16, 17, 18, 19};
Plane Surface(58) = {57};

// surface, internal wall
Line Loop(59) = {34, 29, -35, -23};
Plane Surface(60) = {59};
Line Loop(61) = {30, -37, -24, 36};
Plane Surface(62) = {61};
Line Loop(63) = {31, -39, -25, 38};
Plane Surface(64) = {63};
Line Loop(65) = {35, 32, -36, -26};
Ruled Surface(66) = {65};
Line Loop(67) = {37, 33, -38, -27};
Ruled Surface(68) = {67};
Line Loop(69) = {22, 23, 26, 24, 27, 25};
Plane Surface(70) = {69};
Line Loop(71) = {29, 32, 30, 33, 31, 28};
Plane Surface(72) = {71};

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

// volume
Surface Loop(99) = {98, 88, 96, 90, 92, 94};
Volume(2) = {99};
Surface Loop(101) = {58, 54, 52, 50, 48, 46, 44, 42, 56, 64, 72, 60, 66, 62, 68, 70};
Volume(1) = {99, 101};


Physical Volume(1) = {1};
Physical Volume(2) = {2};
