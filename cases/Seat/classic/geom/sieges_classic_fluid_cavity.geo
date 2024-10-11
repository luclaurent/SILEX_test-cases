// Parameters: acoustic cavity
lx = 5.0;
ly = 2.0;
lz = 2.0;
 
bx1  = 2.0;
lys = 0.7;
h1  = 0.6;
a1 = -0.1;
r  = 0.1;
a2 = 0.6;
h2 = 0.8;
h3 = 0.2;
a3 = 0.1;

bx2 = 4.0;
by2 = ly-lys;

//nf = lx/40;
nf = lx/30;

Mesh.CharacteristicLengthMax=nf;

// Corners
Point(1) = {0,     0  , 0, nf};
Point(2) = {lx,    0  , 0, nf};
Point(3) = {lx,    ly , 0, nf};
Point(4) = {0 ,    ly , 0, nf};

Point(5) = {0,     0  , lz, nf};
Point(6) = {lx,    0  , lz, nf};
Point(7) = {lx,    ly , lz, nf};
Point(8) = {0 ,    ly , lz, nf};


// curved surface on top
Point(300)= {0 ,    ly/2 , 0.8*lz, nf};
Point(301)= {lx ,    ly/2 , 0.8*lz, nf};

Circle(1) = {8, 300, 5};
Circle(2) = {7, 301, 6};

Line(3) = {5, 1};
Line(4) = {2, 1};
Line(5) = {2, 3};
Line(6) = {6, 2};
Line(7) = {7, 3};
Line(8) = {7, 8};
Line(9) = {6, 5};
Line(10) = {8, 4};
Line(11) = {4, 1};
Line(12) = {4, 3};
Line Loop(13) = {11, -4, 5, -12};
Plane Surface(14) = {13};
Line Loop(15) = {11, -3, -1, 10};
Plane Surface(16) = {15};
Line Loop(17) = {9, 3, -4, -6};
Plane Surface(18) = {17};
Line Loop(19) = {6, 5, -7, 2};
Plane Surface(20) = {19};
Line Loop(21) = {12, -7, 8, 10};
Plane Surface(22) = {21};
Line Loop(23) = {9, -1, -8, 2};
Ruled Surface(24) = {23};
Surface Loop(25) = {16, 14, 18, 24, 22, 20};
Volume(26) = {25};

Physical Volume(1) = {26};
Physical Surface(2) = {16};
