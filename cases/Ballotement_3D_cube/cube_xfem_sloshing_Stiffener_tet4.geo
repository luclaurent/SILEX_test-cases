// Parameters: acoustic cavity
lx1 = 1.0;
ly1 = 0.8;
lz1 = 0.6;


lxa = 0.41;
// th = 0.01;
lza = 0.33;


// size of elements
h =  lx1/50;

Mesh.CharacteristicLengthMax=10*h;
Mesh.ElementOrder = 1;

// Cavity: Corners
//Point(1) = {0,     0  , 0, h};
//Point(2) = {lx1,    0  , 0, h};
//Point(3) = {lx1,    ly1 , 0, h};
//Point(4) = {0 ,    ly1 , 0, h};
//Point(5) = {0,     0  , lz1, h};
//Point(6) = {lx1,    0  , lz1, h};
//Point(7) = {lx1,    ly1 , lz1, h};
//Point(8) = {0 ,    ly1 , lz1, h};


// stiffener 
Point(10) = {lxa,   0  , 0, h};
//Point(11) = {lxa+th,   0  , 0, h};
Point(12) = {lxa,   0  , lza , h};
//Point(13) = {lxa+th,   0  , lza , h};

Point(20) = {lxa,    ly1, 0, h};
//Point(21) = {lxa+th, ly1, 0, h};
Point(22) = {lxa,    ly1, lza , h};
//Point(23) = {lxa+th, ly1, lza , h};




// Cavity: lines

//Line(1) = {1, 5};
//Line(2) = {5, 6};
//Line(3) = {6, 2};
//Line(4) = {2, 1};
//Line(5) = {1, 4};
//Line(6) = {4, 8};
//Line(7) = {8, 5};
//Line(8) = {8, 7};
//Line(9) = {7, 3};
//Line(10) = {3, 4};
//Line(11) = {2, 3};
//Line(12) = {6, 7};

//+
Line(13) = {10, 12};
//+
Line(14) = {12, 22};
//+
Line(15) = {22, 20};
//+
Line(16) = {20, 10};
//+
Curve Loop(1) = {16, 13, 14, 15};
//+
Plane Surface(1) = {1};


//+ // stiffener surface
Physical Surface(50) = {1};

//+ // stiffener edge in the fluid
Physical Curve(60) = {14};
