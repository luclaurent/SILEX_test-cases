// Parameters: acoustic cavity
lx1 = 1.0;
ly1 = 0.8;
lz1 = 0.6;

lxa = 0.41;
lxashift = 0.1;
th = 0.002;
lza = 0.33;

// size of elements
h =  lx1/20;

Mesh.CharacteristicLengthMax=10*h;
Mesh.ElementOrder = 2;

// Cavity: Corners
Point(1) = {0,     0  , 0, h};
Point(2) = {lx1,    0  , 0, h};
Point(3) = {lx1,    ly1 , 0, h};
Point(4) = {0 ,    ly1 , 0, h};
Point(5) = {0,     0  , lz1, h};
Point(6) = {lx1,    0  , lz1, h};
Point(7) = {lx1,    ly1 , lz1, h};
Point(8) = {0 ,    ly1 , lz1, h};

Point(10) = {lxa-th/2,   0  , 0, h};
Point(11) = {lxa+th/2,   0  , 0, h};
Point(12) = {lxa-th/2,   0  , lza , h};
Point(13) = {lxa+th/2,   0  , lza , h};

Point(20) = {lxa+lxashift-th/2, ly1, 0,    h};
Point(21) = {lxa+lxashift+th/2, ly1, 0,    h};
Point(22) = {lxa+lxashift-th/2, ly1, lza , h};
Point(23) = {lxa+lxashift+th/2, ly1, lza , h};


// Cavity: lines

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};

Line(5) = {1, 4};
Line(6) = {4, 8};
Line(7) = {8, 5};
Line(8) = {8, 7};
Line(9) = {7, 3};

Line(11) = {2, 3};
Line(12) = {6, 7};

//+
Line(13) = {1, 10};
//+
Line(14) = {10, 12};
//+
Line(15) = {12, 13};
//+
Line(16) = {13, 11};
//+
Line(17) = {11, 2};
//+
Line(18) = {4, 20};
//+
Line(19) = {20, 22};
//+
Line(20) = {22, 23};
//+
Line(21) = {23, 21};
//+
Line(22) = {21, 3};
//+
Line(23) = {10, 20};
//+
Line(24) = {12, 22};
//+
Line(25) = {13, 23};
//+
Line(26) = {11, 21};
//+
Curve Loop(1) = {5, 18, -23, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {23, 19, -24, -14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {-15, -25, 20, 24};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {-16, -26, 21, 25};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {-17, -11, 22, 26};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {11, -9, -12, 3};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {-5, -6, -7, 1};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {-18, -19, -20, -21, -22, 9, 8, 6};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {13, 14, 15, 16, 17, -3, -2, -1};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {2, 12, -8, 7};
//+
Plane Surface(10) = {10};
//+
Surface Loop(1) = {7, 1, 8, 2, 3, 9, 4, 5, 6, 10};
//+
Volume(1) = {1};


//+// Fluid volume
Physical Volume(10) = {1};

//+ // Free fluid surface
Physical Surface(30) = {10};

//+// Structure surface
Physical Surface(20) = {7, 1, 2, 3, 4, 5, 6, 9, 8};


