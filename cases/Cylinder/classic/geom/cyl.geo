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

Line(15) = {5, 11};
Line(16) = {11, 6};
Line(17) = {6, 2};
Line(18) = {9, 11};
Line(19) = {9, 2};
Line(20) = {9, 1};
Line(21) = {1, 5};
Line(22) = {5, 12};
Line(23) = {12, 10};
Line(24) = {10, 4};
Line(25) = {4, 8};
Line(26) = {8, 12};
Line(27) = {8, 7};
Line(28) = {7, 6};
Line(29) = {2, 3};
Line(30) = {3, 4};
Line(31) = {10, 1};
Line Loop(32) = {22, 23, 31, 21};
Plane Surface(33) = {32};
Line Loop(34) = {21, 15, -18, 20};
Plane Surface(35) = {34};
Line Loop(36) = {18, 16, 17, -19};
Plane Surface(37) = {36};
Line(38) = {7, 3};
Line Loop(39) = {38, -29, -17, -28};
Plane Surface(40) = {39};
Line Loop(41) = {38, 30, 25, 27};
Plane Surface(42) = {41};
Line Loop(43) = {25, 26, 23, 24};
Plane Surface(44) = {43};
Line Loop(45) = {24, -30, -29, -19, 13};
Plane Surface(46) = {45};
Line Loop(47) = {26, -14, 16, -28, -27};
Plane Surface(48) = {47};
Line Loop(49) = {22, -14, -15};
Plane Surface(50) = {49};
Line Loop(51) = {20, -31, -13};
Plane Surface(52) = {51};
Line Loop(53) = {23, -13, 18, 14};
Ruled Surface(54) = {53};
Surface Loop(55) = {33, 50, 35, 52, 54};
Volume(56) = {55};
Surface Loop(57) = {48, 44, 42, 40, 46, 37, 54};
Volume(58) = {57};

Physical Volume(1) = {56, 58};
Physical Surface(2) = {54};
