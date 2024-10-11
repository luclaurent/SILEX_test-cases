n=1+10/1;

Point(1) = {0, 50, 0};
Point(2) = {0, 200, 0};
Point(3) = {50, 0, 0};
Point(4) = {100, 0, 0};
Point(5) = {100, 200, 0};
Point(6) = {0, 0, 0};

Point(7) = {0, 100, 0};
Point(8) = {100, 100, 0};

Point(9) = {50*Sqrt(2)/2, 50*Sqrt(2)/2, 0};

Line(1) = {3, 4};
Line(2) = {4, 8};
Line(3) = {8, 5};
Line(4) = {5, 2};
Line(5) = {2, 7};
Line(6) = {7, 1};
Line(7) = {7, 8};
Line(8) = {9, 8};
Circle(9) = {1, 6, 9};
Circle(10) = {9, 6, 3};

Transfinite Line {10, 9, 8, 7, 6, 5, 4, 3, 2, 1} = n Using Progression 1;
Line Loop(11) = {4, 5, 7, 3};
Plane Surface(12) = {11};
Line Loop(13) = {6, 9, 8, -7};
Plane Surface(14) = {13};
Line Loop(15) = {10, 1, 2, -8};
Plane Surface(16) = {15};
Transfinite Surface {12};
Transfinite Surface {14};
Transfinite Surface {16};


Recombine Surface {12, 14, 16};


Physical Line(1) = {6, 5};
Physical Line(2) = {1};
Physical Line(3) = {4};
Physical Surface(4) = {12, 14, 16};
