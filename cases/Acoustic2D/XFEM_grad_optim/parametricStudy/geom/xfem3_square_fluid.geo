// Parameters: acoustic cavity
lx = 5;
ly = 5;

ax = 3.6;
ay = 3.6;
l1 = 0.8;
l2 = 0.8;

l3 = 0.5;
l4 = 0.2;
l5 = 1.0;

// xfem 1
//n1 = lx/3; 

// xfem 2
//n1 = lx/20;

nh=40;
nv=80;
n1=ly/(nv*1.);

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


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};

Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};


Physical Line(2) = {1};
Physical Line(3) = {4};

Line(22) = {2, 3};


Line Loop(26) = {8, 9, 6, 7};
Plane Surface(27) = {26};
Physical Surface(5) = {27};

Plane Surface(15) = {5,-26};
Physical Surface(1) = {15,27};
