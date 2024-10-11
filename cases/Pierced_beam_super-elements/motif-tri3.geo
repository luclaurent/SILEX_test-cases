h=       10;
a=        0;
b=        0;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {100.0, 0.0, 0.0, h};
Point(3) = {100.0, 100.0, 0.0, h};
Point(4) = {0.0, 100.0, 0.0, h};

Point(5) = {50.0, 50.0, 0.0, h};
Point(6) = {50.0+30.0, 50.0,  0.0, h};
Point(7) = {50.0, 50.0+30.0, 0.0, h};
Point(8) = {50.0-30.0, 50.0,  0.0, h};
Point(9) = {50.0, 50.0-30.0, 0.0, h};

Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Line Loop(9) = {8, 5, 6, 7};
Line Loop(10) = {2, 3, 4, 1};
Plane Surface(11) = {9, 10};

Physical Line(2) = {6};
Physical Line(3) = {8};
Physical Surface(1) = {11};
