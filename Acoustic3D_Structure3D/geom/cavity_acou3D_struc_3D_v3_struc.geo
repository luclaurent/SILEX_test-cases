// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.5;

ly2 = 1.0;
lx2 = 3.5;
lx3 = 1.1;

ly4 = 1.2; // structure thickness
R = 1.0; // sphere radius

ly3 = 1;
l4 = 2.5;
lz4 = 1.2;
r4   = 0.5;
e4   = 0.2;
deg=       30;
angle = deg*Pi/180;


// porous
ee = 0.2;

lpy1=0.8;
lpz1=0.6;
lpy2=2.4;
lpz2=1.3;

// size of elements
// air
h = lx1/45;
// air control volume
hc = h/1.5;
// structure
hs = h/2;
// porous
h2 = ee; // /2;


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

// control volume
Point(31) = {lx5-h5/2,     ly5-h5/2  , lz5-h5/2, hc};
Point(32) = {lx5-h5/2,    ly5-h5/2  , lz5+h5/2, hc};
Point(33) = {lx5+h5/2,     ly5-h5/2  , lz5-h5/2, hc};
Point(34) = {lx5+h5/2,    ly5-h5/2  , lz5+h5/2, hc};
Point(35) = {lx5-h5/2,     ly5+h5/2  , lz5-h5/2, hc};
Point(36) = {lx5-h5/2,    ly5+h5/2  , lz5+h5/2, hc};
Point(37) = {lx5+h5/2,     ly5+h5/2  , lz5-h5/2, hc};
Point(38) = {lx5+h5/2,    ly5+h5/2  , lz5+h5/2, hc};
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

// wall: sphere, center (lx3,ly3,0)

Point(100) = {lx3, ly3, 0, hs};


//+ structure surface in contact with air
//Physical Surface(2) = {155, 5155, 5157, 5163, 5156, 5162, 5160, 5161, 5164, 5165, 5158, 5159};


//+
SetFactory("OpenCASCADE");
Sphere(1) = {lx3,ly3, 0, R, -Pi/2, Pi/2, 2*Pi};
//+

//+
Characteristic Length {101, 102} = hs;
//+
Physical Surface(2) = {99};
