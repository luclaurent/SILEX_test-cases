// Parameters: acoustic cavity
lx1 = 1.0;
ly1 = 0.8;
lz1 = 0.6;


// size of elements
h =  lx1*100000;

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


// Cavity: lines

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 1};
Line(5) = {1, 4};
Line(6) = {4, 8};
Line(7) = {8, 5};
Line(8) = {8, 7};
Line(9) = {7, 3};
Line(10) = {3, 4};
Line(11) = {2, 3};
Line(12) = {6, 7};

Line Loop(17) = {-7, 1, -5, -6};
Plane Surface(18) = {17};
Line Loop(19) = {-2, -3, -4, -1};
Plane Surface(20) = {19};
Line Loop(21) = {3, 11, -9, -12};
Plane Surface(22) = {21};
Line Loop(23) = {-11, -10, 5, 4};
Plane Surface(24) = {23};
Line Loop(25) = {2, 12, -8, 7};
Plane Surface(26) = {25};
Line Loop(27) = {8, 9, 10, 6};
Plane Surface(28) = {27};

//Surface Loop(31) = {18, 26, 20, 22, 24, 28};
//Volume(32) = {31};


//Physical Volume(1) = {32}; // air cavity 

Surface Loop(1025) = {18, 26, 20, 22, 24, 28};

Volume(1027) = {1025};

Physical Volume(10) = {1027,1028}; // air cavity 

// Structure surface
Physical Surface(20) = {20, 22, 24, 28, 18}; 

// Free fluid surface
Physical Surface(30) = {26};
