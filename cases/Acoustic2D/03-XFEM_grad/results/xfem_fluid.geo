// Parameters: acoustic cavity
lx = 1.0;
ly = 1.0;

// xfem 1
//n1 = lx/3; 

// xfem 2
n1 = lx/10;


//n1=lx/50;

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

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {4, 1, 2, 3};

Plane Surface(6) = {5};
Physical Surface(1) = {6};
Physical Line(2) = {4};
