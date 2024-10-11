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

//nf1 = lx/45;
nf1=lx/40;
//nf1 = lx/20;
ns = lys/25;
nf2 = nf1/1.5;

lzfin = 1.2*h1+r+r+h2+r;

Mesh.CharacteristicLengthMax=nf1;

// Corners
Point(1) = {0,     0  , 0, nf1};
Point(2) = {lx,    0  , 0, nf1};
Point(3) = {lx,    ly , 0, nf1};
Point(4) = {0 ,    ly , 0, nf1};

Point(5) = {0,     0  , lz, nf1};
Point(6) = {lx,    0  , lz, nf1};
Point(7) = {lx,    ly , lz, nf1};
Point(8) = {0 ,    ly , lz, nf1};

// curved surface on top
Point(300)= {0 ,    ly/2 , 0.8*lz, nf2};
Point(301)= {lx ,    ly/2 , 0.8*lz, nf2};

Point(9) = {0.9*bx1-a1-r-a2-r-a3, 0  , 0, nf2};
Point(10) = {0.9*bx1-a1-r-a2-r-a3, ly  , 0, nf2};
Point(11) = {0.9*bx1-a1-r-a2-r-a3, 0  , lzfin, nf2};
Point(12) = {0.9*bx1-a1-r-a2-r-a3, ly  , lzfin, nf2};

Point(13) = {1.1*bx2, 0  , 0, nf2};
Point(14) = {1.1*bx2, ly  , 0, nf2};
Point(15) = {1.1*bx2, 0  , lzfin, nf2};
Point(16) = {1.1*bx2, ly  , lzfin, nf2};

Circle(58) = {8, 300, 5};
Circle(59) = {7, 301, 6};

Line(1) = {1, 9};
Line(2) = {9, 13};
Line(3) = {13, 2};
Line(4) = {2, 3};
Line(5) = {3, 14};
Line(6) = {14, 13};
Line(7) = {14, 10};
Line(8) = {10, 9};
Line(9) = {1, 4};
Line(10) = {4, 10};
Line(11) = {10, 12};
Line(12) = {12, 11};
Line(13) = {11, 9};
Line(14) = {5, 1};
Line(15) = {8, 5};
Line(16) = {8, 7};
Line(17) = {7, 3};
Line(18) = {7, 6};
Line(19) = {6, 5};
Line(20) = {11, 15};
Line(21) = {15, 13};
Line(22) = {2, 6};
Line(23) = {15, 16};
Line(24) = {16, 12};
Line(25) = {14, 16};
Line(26) = {4, 8};
Line Loop(27) = {14, 9, 26, 15};
Plane Surface(28) = {27};
Line Loop(29) = {26, 16, 17, 5, 25, 24, -11, -10};
Plane Surface(30) = {29};
Line Loop(31) = {24, -11, -7, 25};
Plane Surface(32) = {31};
Line Loop(33) = {17, -4, 22, -18};
Plane Surface(34) = {33};
Line Loop(35) = {16, 18, 19, -15};
Plane Surface(36) = {35};
Line Loop(37) = {12, 13, -8, 11};
Plane Surface(38) = {37};
Line Loop(39) = {21, -6, 25, -23};
Plane Surface(40) = {39};
Line Loop(41) = {6, -2, -8, -7};
Plane Surface(42) = {41};
Line Loop(43) = {8, -1, 9, 10};
Plane Surface(44) = {43};
Line Loop(45) = {6, 3, 4, 5};
Plane Surface(46) = {45};
Line Loop(47) = {22, 19, 14, 1, -13, 20, 21, 3};
Plane Surface(48) = {47};
Line Loop(49) = {21, -2, -13, 20};
Plane Surface(50) = {49};
Surface Loop(51) = {36, 30, 28, 48, 34, 46, 44, 50, 42, 32};
Line Loop(52) = {20, 23, 24, 12};
Plane Surface(53) = {52};

Line Loop(60) = {59, 19, -58, 16};
Ruled Surface(61) = {60};
Line Loop(62) = {58, -15};
Plane Surface(63) = {62};
Line Loop(64) = {59, -18};
Plane Surface(65) = {64};
Surface Loop(66) = {65, 61, 63, 36};
Volume(67) = {66};

Surface Loop(54) = {44, 48, 34, 30, 28, 36, 46, 53, 40, 38};
Volume(55) = {54};
Surface Loop(56) = {50, 42, 32, 38, 53, 40};
Volume(57) = {56};
Physical Volume(3) = {55, 57,67};


