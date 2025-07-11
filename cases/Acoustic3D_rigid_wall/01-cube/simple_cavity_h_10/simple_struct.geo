// Parameters: acoustic cavity
lx1 = 0.6001;
ly1 = 1.0;
lz1 = 1.0;


a = 0.66;

// size of elements
h =  lx1/10;
h2 = lz1/10;

//h  = lx1*10;
//h2 = lz1*10;

Mesh.CharacteristicLengthMax=10*h;
Mesh.ElementOrder = 1;

// rectangle: Corners
Point(1) = {lx1,    0  , 0, h};
Point(2) = {lx1,    ly1  , 0, h};
Point(3) = {lx1,    ly1 , lz1, h};
Point(4) = {lx1,    0 , lz1, h};

// Cavity: lines

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(10) = {1,2,3,4};
Plane Surface(11) = {10};

Physical Surface(2) = {11};


