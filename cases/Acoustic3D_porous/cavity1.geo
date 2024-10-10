// Parameters: acoustic cavity
lx1 = 0.75;
ly1 = 0.6;
lz1 = 0.4;

lx3 = 0.05;

h = lx3/4;
h2 = 2*h;

Mesh.CharacteristicLengthMax=1.5*h2;
Mesh.ElementOrder = 1;

// Corners
Point(1) = {0,     0  , 0, h};
Point(2) = {lx1,    0  , 0, h2};
Point(3) = {lx1,    ly1 , 0, h2};
Point(4) = {0 ,    ly1 , 0, h};

Point(5) = {0,     0  , lz1, h};
Point(6) = {lx1,    0  , lz1, h2};
Point(7) = {lx1,    ly1 , lz1, h2};
Point(8) = {0 ,    ly1 , lz1, h};

Point(9) = {lx3,    0  , 0, h};
Point(10) = {lx3,    ly1 , 0, h};
Point(11) = {lx3,    0  , lz1, h};
Point(12) = {lx3,    ly1 , lz1, h};


Line(1) = {1, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 10};
Line(5) = {10, 4};
Line(6) = {4, 1};
Line(7) = {9, 10};
Line(8) = {11, 12};
Line(9) = {8, 5};
Line(10) = {5, 11};
Line(11) = {8, 12};
Line(12) = {8, 4};
Line(13) = {12, 10};
Line(14) = {5, 1};
Line(15) = {11, 9};
Line(16) = {7, 3};
Line(17) = {6, 2};
Line(18) = {6, 7};
Line(19) = {7, 12};
Line(20) = {11, 6};


Line Loop(21) = {12, 6, -14, -9};
Plane Surface(22) = {21};
Line Loop(23) = {13, -7, -15, 8};
Plane Surface(24) = {23};
Line Loop(25) = {16, -3, -17, 18};
Plane Surface(26) = {25};
Line Loop(27) = {3, 4, -7, 2};
Plane Surface(28) = {27};
Line Loop(29) = {7, 5, 6, 1};
Plane Surface(30) = {29};
Line Loop(31) = {16, 4, -13, -19};
Plane Surface(32) = {31};
Line Loop(33) = {13, 5, -12, 11};
Plane Surface(34) = {33};
Line Loop(35) = {18, 19, -8, 20};
Plane Surface(36) = {35};
Line Loop(37) = {8, -11, 9, 10};
Plane Surface(38) = {37};
Line Loop(39) = {15, -1, -14, 10};
Plane Surface(40) = {39};
Line Loop(41) = {2, -17, -20, 15};
Plane Surface(42) = {41};


Surface Loop(43) = {26, 32, 28, 42, 36, 24};
Volume(44) = {43};
Surface Loop(45) = {38, 34, 30, 22, 40, 24};
Volume(46) = {45};


Physical Volume(1) = {44}; // air
Physical Volume(2) = {46}; // porous
Physical Surface(3) = {24}; // interface
Physical Surface(4) = {22}; // porous external surface: large one, normal x

Physical Surface(5) = {38, 30}; // 2 porous external surfaces: normal z
Physical Surface(6) = {40, 34}; // 2 porous external surfaces: normal y
