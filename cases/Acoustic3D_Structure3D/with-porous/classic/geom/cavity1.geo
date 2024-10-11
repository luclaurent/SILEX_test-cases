// Parameters: acoustic cavity
lx1 = 7.0;
ly1 = 4.0;
lz1 = 2.8;

ly2 = 3.0;
lx2 = 3.5;

h = lx1/30;

Mesh.CharacteristicLengthMax=h;
Mesh.ElementOrder = 1;

// Corners
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

// lines
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

//surfaces
Line Loop(22) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(23) = {22};
Line Loop(24) = {1, 9, -15, -8};
Plane Surface(25) = {24};
Line Loop(26) = {2, 10, -16, -9};
Plane Surface(27) = {26};
Line Loop(28) = {3, 11, -17, -10};
Plane Surface(29) = {28};
Line Loop(30) = {4, 12, -18, -11};
Plane Surface(31) = {30};
Line Loop(32) = {5, 13, -19, -12};
Plane Surface(33) = {32};
Line Loop(34) = {6, 14, -20, -13};
Plane Surface(35) = {34};
Line Loop(36) = {7, 8, -21, -14};
Plane Surface(37) = {36};
Line Loop(38) = {21, 15, 16, 17, 18, 19, 20};
Plane Surface(39) = {38};

//volume
Surface Loop(40) = {25, 23, 27, 29, 31, 33, 35, 37, 39};
Volume(41) = {40};

//physical groups
Physical Volume(1) = {41};
