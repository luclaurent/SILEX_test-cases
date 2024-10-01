// Parameters: acoustic cavity
lx = 1.0;
ly = 1.0;
 
ax = 0.6;
hy = 0.65;


ep = 0.0004;

n1 = lx/160;
n2 = ep/30;

//n1 = lx/10;
//n2 = ep/10;

// Corners
Point(1)  = {0,     0  , 0, n1};
Point(2)  = {lx,    0  , 0, n1};
Point(3)  = {lx,    ly , 0, n1};
Point(4)  = {0 ,    ly , 0, n1};

Point(9)  = {ax-ep/2,  0  , 0, n1};
Point(10) = {ax+ep/2,  0  , 0, n1};
Point(11) = {ax-ep/2,  hy , 0, n2};
Point(12) = {ax+ep/2,  hy , 0, n2};

Point(17) = {ax,  hy , 0, n2};
Point(19) = {ax,  hy+ep/2 , 0, n2};

Circle(1) = {11, 17, 19};
Circle(2) = {19, 17, 12};

Line(5) = {1, 9};
Line(6) = {9, 11};
Line(7) = {12, 10};
Line(8) = {10, 2};
Line(9) = {2, 3};
Line(10) = {3, 4};
Line(11) = {4, 1};

Line Loop(12) = {11, 5, 6, 1, 2, 7, 8, 9, 10};
Plane Surface(13) = {12};
Physical Surface(1) = {13};
Physical Line(2) = {11};
Physical Line(3) = {1, 2};
