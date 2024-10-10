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
e4   = 0.2;
deg=       30;
angle = deg*Pi/180;


// porous
ee = 0.20;

// size of elements
h = lx1/60;
h2 = ee/3;


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
//Line Loop(87) = {75, 80, -84, -79};
//Plane Surface(88) = {87};
//Line Loop(89) = {76, 85, -81, -80};
//Plane Surface(90) = {89};
//Line Loop(91) = {77, 86, -82, -85};
//Plane Surface(92) = {91};
//Line Loop(93) = {78, 79, -83, -86};
//Plane Surface(94) = {93};
//Line Loop(95) = {75, 76, 77, 78};
//Plane Surface(96) = {95};
//Line Loop(97) = {84, 81, 82, 83};
//Plane Surface(98) = {97};

// wall: corners
Point(15) = {lx3,     ly3  , 0, h/3};
Point(16) = {lx3+l4,     ly3  , 0, h/3};
Point(17) = {lx3,     ly3  , lz4-r4 , h/3};
Point(18) = {lx3+l4,     ly3  , lz4-r4, h/3};
Point(19) = {lx3+r4,     ly3  , lz4-r4 , h/3};
Point(20) = {lx3+l4-r4,     ly3  , lz4-r4, h/3};
Point(21) = {lx3+r4,     ly3  , lz4 , (h2+h)/2};
Point(22) = {lx3+l4-r4,     ly3  , lz4, (h2+h)/2};

// rotate wall
Rotate {{0, 0, 1}, {lx3+l4/2, ly3, 0}, angle} {
  Point{ 18, 20, 21, 19, 17, 16, 15, 22};
}

Line(140) = {15, 16};
Line(141) = {16, 18};
Line(142) = {22, 21};
Line(143) = {17, 15};
Circle(144) = {18, 20, 22};
Circle(145) = {21, 19, 17};
Line Loop(154) = {141, 144, 142, 145, 143, 140};
Plane Surface(155) = {154};

// volume




//Line Loop(113) = {6, 14, -20, -13};
//Ruled Surface(114) = {113};
//Line Loop(115) = {5, 13, -19, -12};
//Ruled Surface(116) = {115};
//Line Loop(117) = {4, 12, -18, -11};
//Ruled Surface(118) = {117};
//Line Loop(119) = {3, 11, -17, -10};
//Ruled Surface(120) = {119};
//Line Loop(121) = {2, 10, -16, -9};
//Ruled Surface(122) = {121};
//Line Loop(123) = {1, 9, -15, -8};
//Ruled Surface(124) = {123};
//Line Loop(125) = {7, 8, -21, -14};
//Ruled Surface(126) = {125};
//Line Loop(127) = {7, 1, 2, 3, 4, 5, 6};
//Line Loop(130) = {21, 15, 16, 17, 18, 19, 20};
//Plane Surface(131) = {130};

// Porous media
//Translate {0, 0, ee} {
//  Duplicata { Surface{131}; }
//}
//Line(256) = {12, 48};
//Line(257) = {5, 40};
//Line(258) = {6, 39};
//Line(259) = {7, 60};
//Line(260) = {14, 56};
//Line(261) = {13, 52};
//Line Loop(262) = {21, 257, -157, -258};
//Plane Surface(263) = {262};
//Line Loop(264) = {258, -163, -259, 20};
//Plane Surface(265) = {264};
//Line Loop(266) = {259, -162, -260, 19};
//Plane Surface(267) = {266};
//Line Loop(268) = {260, -161, -261, 18};
//Plane Surface(269) = {268};
//Line Loop(270) = {17, 261, -160, -256};
//Plane Surface(271) = {270};
//Line Loop(272) = {15, 16, 256, -159, -158, -257};
//Plane Surface(273) = {272};

//Plane Surface(274) = {127};


//Surface Loop(275) = {96, 88, 90, 92, 94, 98};
//Volume(276) = {275};

//Surface Loop(277) = {124, 274, 126, 114, 116, 118, 120, 122, 131};
//Volume(278) = {275, 277};

//Surface Loop(279) = {271, 269, 267, 265, 263, 273, 156, 131};
//Volume(280) = {279};

//Physical Volume(1) = {278,276}; // air cavity + controlled air volume

//Physical Volume(5) = {276}; // controlled air volume (small)

//Physical Volume(2) = {280}; // porous volume

//Physical Surface(3) = {131}; // [porous]-[air cavity] : interface surface

//Physical Surface(4) = {265, 267, 269, 271, 273, 263, 156};// porous external surface

Physical Surface(6) = {155}; // Structure surface

Physical Line(7) = {143, 145, 142, 144, 141};// Structure edge of the free boundary in air

