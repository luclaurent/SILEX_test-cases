// Parameters: acoustic cavity
lx1 = 1.0;
ly1 = 1.0;
lz1 = 1.0;


a = 0.66;

// size of elements
h =  lx1/30;
//h2 = lz1/20;

//h  = lx1*10;
//h2 = lz1*10;

Mesh.CharacteristicLengthMax=10*h;
Mesh.ElementOrder = 1;

// Cavity: Corners
Point(1) = {0,     0  , 0, h};
Point(2) = {lx1,    0  , 0, h};
Point(3) = {lx1,    ly1 , 0, h};
Point(4) = {0 ,    ly1 , 0, h};
Point(5) = {0,     0  , lz1, h};
Point(6) = {lx1,    0  , lz1, h};
Point(7) = {lx1,    ly1 , lz1, h};
Point(8) = {0 ,    ly1 , lz1, h};

xc= 0.7;
yc= 0.6;
zc= 0.5;
lxc = 0.2;
lyc = 0.2;
lzc = 0.2;
//hc = lxc/4;
hc = h/1.5;
// Cavity: control volume
Point(101) = {xc,     yc  , zc, hc};
Point(102) = {xc+lxc,    yc  , zc, hc};
Point(103) = {xc+lxc,    yc+lyc , zc, hc};
Point(104) = {xc ,    yc+lyc , zc, hc};
Point(105) = {xc,     yc  , zc+lzc, hc};
Point(106) = {xc+lxc,    yc  , zc+lzc, hc};
Point(107) = {xc+lxc,    yc+lyc , zc+lzc, hc};
Point(108) = {xc ,    yc+lyc , zc+lzc, hc};

Line(101) = {101, 105};
Line(102) = {105, 106};
Line(103) = {106, 102};
Line(104) = {102, 101};
Line(105) = {101, 104};
Line(106) = {104, 108};
Line(107) = {108, 105};
Line(108) = {108, 107};
Line(109) = {107, 103};
Line(1010) = {103, 104};
Line(1011) = {102, 103};
Line(1012) = {106, 107};

// structure
//Point(9) = {a,     0  , 0, h2};
//Point(10) = {a,     0  , lz1, h2};
//Point(11) = {a,    ly1 , lz1, h2};
//Point(12) = {a ,    ly1 , 0, h2};

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

//Line(13) = {10, 11};
//Line(14) = {11, 12};
//Line(15) = {12, 9};
//Line(16) = {9, 10};

Line Loop(17) = {7, -1, 5, 6};
Plane Surface(18) = {17};
Line Loop(19) = {2, 3, 4, 1};
Plane Surface(20) = {19};
Line Loop(21) = {3, 11, -9, -12};
Plane Surface(22) = {21};
Line Loop(23) = {11, 10, -5, -4};
Plane Surface(24) = {23};
Line Loop(25) = {2, 12, -8, 7};
Plane Surface(26) = {25};
Line Loop(27) = {8, 9, 10, 6};
Plane Surface(28) = {27};
//Line Loop(29) = {16, 13, 14, 15};
//Plane Surface(30) = {29};

//Surface Loop(31) = {18, 26, 20, 22, 24, 28};
//Volume(32) = {31};

// surfaces control volume
Line Loop(1013) = {1010, 106, 108, 109};
Plane Surface(1014) = {1013};
Line Loop(1015) = {105, -1010, -1011, 104};
Plane Surface(1016) = {1015};
Line Loop(1017) = {105, 106, 107, -101};
Plane Surface(1018) = {1017};
Line Loop(1019) = {101, 102, 103, 104};
Plane Surface(1020) = {1019};
Line Loop(1021) = {103, 1011, -109, -1012};
Plane Surface(1022) = {1021};
Line Loop(1023) = {1012, -108, 107, 102};
Plane Surface(1024) = {1023};

//Physical Volume(1) = {32}; // air cavity 

//Physical Surface(6) = {30};// Structure surface



Surface Loop(1025) = {18, 26, 20, 22, 24, 28};
Surface Loop(1026) = {1024, 1022, 1020, 1018, 1016, 1014};
Volume(1027) = {1025,1026};
//Volume(1027) = {1025};
Volume(1028) = {1026};

Physical Volume(1) = {1027,1028}; // air cavity 

Physical Volume(5) = {1028}; // control volume 

//+
Physical Surface(6) = {24};
//+


