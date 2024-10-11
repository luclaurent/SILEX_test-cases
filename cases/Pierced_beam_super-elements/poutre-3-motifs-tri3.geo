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

Point(10) = {200.0, 0.0   , 0.0, h};
Point(11) = {200.0, 100.0 , 0.0, h};
Point(12) = {300.0, 0.0   , 0.0, h};
Point(13) = {300.0, 100.0 , 0.0, h};

Point(15) = {150.0, 50.0, 0.0, h};
Point(16) = {150.0+30.0, 50.0,  0.0, h};
Point(17) = {150.0, 50.0+30.0, 0.0, h};
Point(18) = {150.0-30.0, 50.0,  0.0, h};
Point(19) = {150.0, 50.0-30.0, 0.0, h};

Point(25) = {250.0, 50.0, 0.0, h};
Point(26) = {250.0+30.0, 50.0,  0.0, h};
Point(27) = {250.0, 50.0+30.0, 0.0, h};
Point(28) = {250.0-30.0, 50.0,  0.0, h};
Point(29) = {250.0, 50.0-30.0, 0.0, h};


Circle(1) = {6, 5, 7};
Circle(2) = {7, 5, 8};
Circle(3) = {8, 5, 9};
Circle(4) = {9, 5, 6};

Line(5) = {1, 2};
Line(7) = {3, 4};
Line(8) = {4, 1};


Line(12) = {2, 10};
Line(13) = {10, 12};
Line(14) = {12, 13};
Line(15) = {13, 11};
Line(16) = {11, 3};
Circle(17) = {26, 25, 27};
Circle(18) = {27, 25, 28};
Circle(19) = {28, 25, 29};
Circle(20) = {29, 25, 26};
Circle(21) = {16, 15, 17};
Circle(22) = {17, 15, 18};
Circle(23) = {18, 15, 19};
Circle(24) = {19, 15, 16};
Line Loop(25) = {7, 8, 5, 12, 13, 14, 15, 16};
Line Loop(26) = {2, 3, 4, 1};
Line Loop(27) = {22, 23, 24, 21};
Line Loop(28) = {18, 19, 20, 17};
Plane Surface(29) = {25, 26, 27, 28};


Physical Line(2) = {14};
Physical Line(3) = {8};
Physical Surface(1) = {29};
