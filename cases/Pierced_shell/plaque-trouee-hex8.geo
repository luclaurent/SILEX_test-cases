cl1 = 2;
ep=cl1;

n = 1+10/cl1;

Point(1) = {0, 50, 0};
Point(2) = {0, 200, 0};
Point(3) = {50, 0, 0};
Point(4) = {100, 0, 0};
Point(5) = {100, 200, 0};
Point(6) = {0, 0, 0};
Point(7) = {0, 50, ep};
Point(8) = {0, 200, ep};
Point(9) = {50, 0, ep};
Point(10) = {100, 0, ep};
Point(11) = {100, 200, ep};
Point(12) = {0, 0, ep};

Point(13) = {0, 100, 0};
Point(14) = {0, 100, ep};
Point(15) = {100, 100, 0};
Point(16) = {100, 100, ep};

Point(17) = {50*Sqrt(2)/2, 50*Sqrt(2)/2, 0};
Point(18) = {50*Sqrt(2)/2, 50*Sqrt(2)/2, ep};

Circle(1) = {17, 6, 1};
Circle(2) = {18, 12, 7};
Circle(3) = {17, 6, 3};
Circle(4) = {18, 12, 9};


Line(5) = {3, 4};

Line(6) = {4, 15};
Line(7) = {15, 5};
Line(8) = {5, 2};
Line(9) = {2, 13};
Line(10) = {13, 1};
Line(11) = {13, 15};
Line(12) = {11, 8};
Line(13) = {8, 14};
Line(14) = {14, 7};
Line(15) = {9, 10};
Line(16) = {10, 16};
Line(17) = {16, 11};
Line(18) = {16, 14};
Line(19) = {16, 18};
Line(20) = {16, 15};
Line(21) = {10, 4};
Line(22) = {9, 3};
Line(23) = {18, 17};
Line(24) = {7, 1};
Line(25) = {14, 13};
Line(26) = {8, 2};
Line(27) = {11, 5};
Line(44) = {17, 15};

Line Loop(30) = {7, 8, 9, 11};
Plane Surface(31) = {30};
Line Loop(32) = {8, -26, -12, 27};
Plane Surface(33) = {32};
Line Loop(34) = {9, -25, -13, 26};
Plane Surface(35) = {34};
Line Loop(36) = {7, -27, -17, 20};
Plane Surface(37) = {36};
Line Loop(38) = {17, 12, 13, -18};
Plane Surface(39) = {38};
Line Loop(40) = {10, -24, -14, 25};
Plane Surface(41) = {40};
Line Loop(42) = {18, 14, -2, -19};
Plane Surface(43) = {42};

Line Loop(45) = {16, 19, 4, 15};
Plane Surface(46) = {45};
Line Loop(47) = {6, -20, -16, 21};
Plane Surface(48) = {47};
Line Loop(49) = {5, -21, -15, 22};
Plane Surface(50) = {49};
Line Loop(51) = {1, -24, -2, 23};
Ruled Surface(52) = {51};
Line Loop(53) = {22, -3, -23, 4};
Ruled Surface(54) = {53};

Line Loop(63) = {44, -20, 19, 23};
Plane Surface(64) = {63};


Line Loop(59) = {6, -44, 3, 5};
Plane Surface(60) = {59};
Line Loop(61) = {44, -11, 10, -1};
Plane Surface(62) = {61};
Line Loop(65) = {18, 25, 11, -20};
Plane Surface(66) = {65};
Surface Loop(67) = {66, 35, 31, 37, 33, 39};
Volume(68) = {67};
Surface Loop(69) = {41, 62, 52, 43, 64, 66};
Volume(70) = {69};
Surface Loop(71) = {60, 48, 46, 54, 50, 64};
Volume(72) = {71};

Transfinite Line {4,5,6,7,8,9,10,1,3,44,15,16,17,12,13,14,2,11,18,19} = n Using Progression 1;
Transfinite Line {22,21,20,27,26,25,24,23} = 2 Using Progression 1;
Transfinite Surface {31,33,35,37,39,41,43,46,48,50,52,54,64,60,62,66};
Recombine Surface {31,33,35,37,39,41,43,46,48,50,52,54,64,60,62,66};
Transfinite Volume {68,70,72};


Physical Surface(1) = {35, 41};
Physical Surface(2) = {50};
Physical Surface(3) = {31, 62, 60};
Physical Surface(4) = {33};
Physical Volume(5) = {70, 72, 68};

