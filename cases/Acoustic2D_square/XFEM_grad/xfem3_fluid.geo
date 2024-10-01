// Parameters: acoustic cavity
lx = 1.0;
ly = 1.0;

ax = 0.2;
ay = 0.2;
l1 = 0.2;
l2 = 0.2;

l3 = 0.5;
l4 = 0.2;
l5 = 1.0;

// xfem 1
//n1 = lx/3; 

// xfem 2
//n1 = lx/20;

nh=10;
nv=20;
n1=lx/(nv*1.5);

// xfem 3
//n1 = lx/27;

// xfem 4
//n1 = lx/81;

// xfem 5
//n1 = lx/243;

// xfem 6
//n1 = lx/350;

// Corners
Point(1)  = {0,     0  , 0, n1};
Point(2)  = {lx,    0  , 0, n1};
Point(3)  = {lx,    ly , 0, n1};
Point(4)  = {0 ,    ly , 0, n1};

Point(5)  = {ax,    ay    , 0, n1};
Point(6)  = {ax+l1, ay    , 0, n1};
Point(7)  = {ax+l1, ay+l2 , 0, n1};
Point(8)  = {ax ,   ay+l2 , 0, n1};

Point(9)  = {l3,     0  , 0, n1};
Point(10) = {l3+l4,    0  , 0, n1};
Point(11) = {l3+l4,   l5 , 0, n1};
Point(12) = {l3 ,    l5 , 0, n1};


Line(1) = {1, 9};
Line(20) = {9,10};

Line(3) = {4, 12};
Line(30) = {3, 11};
Line(4) = {4, 1};

Line Loop(5) = {4, 1, 20, 2, 3};

Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};


Physical Line(2) = {4};


Line(21) = {10, 2};
Line(22) = {2, 3};
Line(23) = {10, 11};
Line(24) = {11, 12};
Line(25) = {12, 9};

Transfinite Line {20, 24 } = nh Using Progression 1;
Transfinite Line {23, 25 } = nv Using Progression 1;

Line Loop(26) = {8, 9, 6, 7};
Plane Surface(27) = {26};
Physical Surface(5) = {27};

Line Loop(29) = {20, 23, 24, 25};
Plane Surface(30) = {29};
Transfinite Surface {30};
//Line Loop(31) = {21, 22, 3, 4, 1, -25, -24, -23};
//Plane Surface(32) = {26, 31};
//Physical Surface(1) = {32, 30, 27};


Line Loop(31) = {22, 30, -23, 21};
Plane Surface(32) = {31};
Line Loop(33) = {25, -1, -4, 3};
Plane Surface(34) = {26, 33};
Physical Surface(1) = {32, 30, 34, 27};
