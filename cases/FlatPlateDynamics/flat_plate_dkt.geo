Lx=0.5;
Ly=0.5;
cl1 = Lx/50;

xpta = 0.1;
ypta = 0.2;

xptb = 0.3;
yptb = 0.4;

Point(1)  = {0,   0, 0, cl1};
Point(121) = {xpta, 0, 0, cl1};
Point(122) = {xptb, 0, 0, cl1};
Point(2)  = {Lx, 0, 0, cl1};
Point(23)  = {Lx, ypta, 0, cl1};
Point(232)  = {Lx, yptb, 0, cl1};
Point(3)  = {Lx, Ly, 0, cl1};
Point(34)  = {xpta, Ly, 0, cl1};
Point(342)  = {xptb, Ly, 0, cl1};
Point(4)  = {0, Ly, 0, cl1};
Point(14)  = {0, ypta, 0, cl1};
Point(142)  = {0, yptb, 0, cl1};

Point(10)  = {xpta, ypta, 0, cl1};

Point(101)  = {xpta, yptb, 0, cl1};
Point(102)  = {xptb, ypta, 0, cl1};

Point(20)  = {xptb, yptb, 0, cl1};

//+
Line(1) = {1, 121};
//+
Line(2) = {121, 122};
//+
Line(3) = {122, 2};
//+
Line(4) = {2, 23};
//+
Line(5) = {23, 232};
//+
Line(6) = {232, 3};
//+
Line(7) = {3, 342};
//+
Line(8) = {342, 34};
//+
Line(9) = {34, 4};
//+
Line(10) = {4, 142};
//+
Line(11) = {142, 14};
//+
Line(12) = {14, 1};
//+
Line(13) = {121, 10};
//+
Line(14) = {10, 101};
//+
Line(15) = {101, 34};
//+
Line(16) = {122, 102};
//+
Line(17) = {102, 20};
//+
Line(18) = {20, 342};
//+
Line(19) = {14, 10};
//+
Line(20) = {10, 102};
//+
Line(21) = {102, 23};
//+
Line(22) = {142, 101};
//+
Line(23) = {101, 20};
//+
Line(24) = {20, 232};
//+
Curve Loop(1) = {1, 13, -19, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 16, -20, -13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, -21, -16};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {19, 14, -22, 11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {20, 17, -23, -14};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {21, 5, -24, -17};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {24, 6, 7, -18};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {23, 18, 8, -15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {22, 15, 9, 10};
//+
Plane Surface(9) = {9};
//+
Physical Surface(1) = {9, 4, 1, 8, 5, 2, 7, 6, 3};
//+
Physical Curve(2) = {4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3};
