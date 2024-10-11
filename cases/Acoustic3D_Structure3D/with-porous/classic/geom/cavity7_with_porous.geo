// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.5;

ly2 = 1.0;
lx2 = 3.5;

lx3 =         2;
ly3 = 2.0;
l4 = 2.5;
lz4 = 2.0;
r4   = 0.5;
e4   = 0.05;
deg=       30;
angle = deg*Pi/180;


// porous
ee = 0.20;

// size of elements
h = lx1/25;
h2 = ee/1;


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

Point(5) = {0,     0  , lz1, h2};
Point(6) = {lx1,    0  , lz1, h2};
Point(7) = {lx1,    ly1 , lz1, h2};
Point(8) = {0 ,    ly1 , lz1, h2};

Point(9) = {0,     ly1+ly2  , 0, h};
Point(10) = {lx2,     ly1+ly2  , 0, h};
Point(11) = {lx2,     ly1  , 0, h};
Point(12) = {0,     ly1+ly2  , lz1, h2};
Point(13) = {lx2,     ly1+ly2  , lz1, h2};
Point(14) = {lx2,     ly1  , lz1, h2};

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
Point(15) = {lx3,     ly3-e4/2  , 0, h/3};
Point(16) = {lx3+l4,     ly3-e4/2  , 0, h/3};
Point(17) = {lx3,     ly3-e4/2  , lz4-r4 , h/3};
Point(18) = {lx3+l4,     ly3-e4/2  , lz4-r4, h/3};
Point(19) = {lx3+r4,     ly3-e4/2  , lz4-r4 , h/3};
Point(20) = {lx3+l4-r4,     ly3-e4/2  , lz4-r4, h/3};
Point(21) = {lx3+r4,     ly3-e4/2  , lz4 , (h2+h)/2};
Point(22) = {lx3+l4-r4,     ly3-e4/2  , lz4, (h2+h)/2};

Point(23) = {lx3,     ly3+e4/2  , 0, h/3};
Point(24) = {lx3+l4,     ly3+e4/2  , 0, h/3};
Point(25) = {lx3,     ly3+e4/2  , lz4-r4 , h/3};
Point(26) = {lx3+l4,     ly3+e4/2  , lz4-r4, h/3};
Point(27) = {lx3+r4,     ly3+e4/2  , lz4-r4 , h/3};
Point(28) = {lx3+l4-r4,     ly3+e4/2  , lz4-r4, h/3};
Point(29) = {lx3+r4,     ly3+e4/2  , lz4 , (h2+h)/2};
Point(30) = {lx3+l4-r4,     ly3+e4/2  , lz4, (h2+h)/2};

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



// control volume
Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, h/2.5};
Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, h/2.5};
Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, h/2.5};
Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, h/2.5};
Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, h/2.5};
Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, h/2.5};
Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, h/2.5};
Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, h/2.5};
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


Line Loop(99) = {25, 39, -31, -38};
Plane Surface(100) = {99};
Line Loop(101) = {25, 22, 23, 26, 24, 27};
Plane Surface(102) = {101};
Line Loop(103) = {24, 37, -30, -36};
Plane Surface(104) = {103};
Line Loop(105) = {23, 35, -29, -34};
Plane Surface(106) = {105};
Line Loop(107) = {31, 28, 29, 32, 30, 33};
Plane Surface(108) = {107};
Line Loop(109) = {38, -33, -37, 27};
Ruled Surface(110) = {109};
Line Loop(111) = {26, 36, -32, -35};
Ruled Surface(112) = {111};
Line Loop(113) = {6, 14, -20, -13};
Ruled Surface(114) = {113};
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
Line Loop(128) = {28, -34, -22, 39};
Plane Surface(129) = {127, 128};
Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};
Plane Surface(131) = {130};



Translate {0, 0, ee} {
  Duplicata { Surface{131}; }
}


Line(140) = {6, 39};
Line(141) = {7, 60};
Line(142) = {14, 56};
Line(143) = {13, 52};
Line(144) = {12, 48};
Line(145) = {5, 40};
Line(146) = {8, 44};
Line Loop(147) = {21, 145, -133, -140};
Plane Surface(148) = {147};
Line Loop(149) = {20, 140, -139, -141};
Plane Surface(150) = {149};
Line Loop(151) = {19, 141, -138, -142};
Plane Surface(152) = {151};
Line Loop(153) = {18, 142, -137, -143};
Plane Surface(154) = {153};
Line Loop(155) = {17, 143, -136, -144};
Plane Surface(156) = {155};
Line Loop(157) = {16, 144, -135, -146};
Plane Surface(158) = {157};
Line Loop(159) = {15, 146, -134, -145};
Plane Surface(160) = {159};

Surface Loop(161) = {129, 126, 124, 122, 120, 118, 116, 114, 100, 102, 106, 112, 104, 110, 108, 131};
Surface Loop(162) = {94, 96, 88, 90, 92, 98};
Volume(1) = {161, 162};
Volume(2) = {162};

Physical Volume(1) = {1,2}; // air cavity + controlled air volume
Physical Volume(5) = {2}; // controlled air volume (small)

Surface Loop(163) = {132, 148, 160, 158, 156, 154, 152, 150, 131};
Volume(3) = {163};

Physical Volume(2) = {3}; // porous volume

Physical Surface(3) = {131}; // [porous]-[air cavity] : interface surface

Physical Surface(4) = {148, 150, 152, 156, 158, 160, 132, 154};// porous external surface
