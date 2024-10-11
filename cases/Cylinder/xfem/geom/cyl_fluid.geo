// Parameters: acoustic cavity
lx = 1;
ly = 1;
lz = 0.03;
rs  = 0.8;
i  = 70;
nf  = lx/i;
is = 70;
ns  = rs/is;
// ns=nf;

// Corners of fluid cavity
Point(1) = {0,  0 , 0, nf};
Point(2) = {lx,  0, 0, nf};
Point(3) = {lx, ly, 0, nf};
Point(4) = {0,  ly, 0, nf};

Point(5) = {0,  0, lz, nf};
Point(6) = {lx, 0, lz, nf};
Point(7) = {lx, ly, lz, nf};
Point(8) = {0,  ly, lz, nf};


// Structure
Point(9)   = {rs,    0, 0, ns};
Point(10)  = {0,   rs, 0, ns};
Point(11)  = {rs,   0, lz, ns};
Point(12)  = {0,   rs, lz, ns};

Circle(13) = {9, 1, 10};
Circle(14) = {11, 5, 12};


Line(15) = {10, 12};
Line(16) = {11, 9};
Line(17) = {6, 2};
Line(18) = {7, 3};
Line(19) = {8, 4};
Line(20) = {5, 1};
Line(21) = {5, 6};
Line(22) = {2, 1};
Line(23) = {1, 4};
Line(24) = {8, 5};
Line(25) = {4, 3};
Line(26) = {7, 8};
Line(27) = {7, 6};
Line(28) = {2, 3};
Line Loop(29) = {24, 20, 23, -19};
Plane Surface(30) = {29};
Line Loop(31) = {19, 25, -18, 26};
Plane Surface(32) = {31};
Line Loop(33) = {18, -28, -17, -27};
Plane Surface(34) = {33};
Line Loop(35) = {22, -20, 21, 17};
Plane Surface(36) = {35};
Line Loop(37) = {22, 23, 25, -28};
Plane Surface(38) = {37};
Line Loop(39) = {27, -21, -24, -26};
Plane Surface(40) = {39};
Line Loop(41) = {15, -14, 16, 13};
Ruled Surface(42) = {41};
Surface Loop(43) = {40, 34, 32, 30, 36, 38};
Volume(44) = {43};
Physical Volume(1) = {44};
// Physical Surface(46) = {42};
