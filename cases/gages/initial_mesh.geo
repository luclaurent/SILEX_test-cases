n1=       10;
n2=       10;
a=        0;
b=        0;
Point(1) = {0.0, 0.0, 0.0, n1};
Point(2) = {192.0, 0.0, 0.0, n1};
Point(3) = {192.0, 40.0, 0.0, n1};
Point(4) = {0.0, 40.0, 0.0, n1};
Point(5) = {58.5-b, 24.0+a, 0.0, n1};
Point(6) = {58.5-b, 16.0-a, 0.0, n1};
Point(7) = {133.5+b, 24.0+a, 0.0, n1};
Point(8) = {133.5+b, 16.0-a, 0.0, n1};
Point(9) = {47.5-b, 24.0+a, 0.0, n2};
Point(10) = {47.5-b, 16.0-a, 0.0, n2};
Point(11) = {144.5+b, 24.0+a, 0.0, n2};
Point(12) = {144.5+b, 16.0-a, 0.0, n2};
Point(13) = {58.5-b, 0.0, 0.0, n2};
Point(14) = {58.5, 40.0, 0.0, n2};
Point(15) = {133.5+b, 0.0, 0.0, n2};
Point(16) = {133.5+b, 40.0, 0.0, n2};
Line(1) = {1, 13};
Line(2) = {13, 15};
Line(3) = {15, 2};
Line(4) = {2, 3};
Line(5) = {3, 16};
Line(6) = {16, 14};
Line(7) = {14, 4};
Line(8) = {4, 1};
x = Sqrt(11*11-6*6);
Point(17) = {58.5+x-b, 24.0+6.0+a, 0.0, n2};
Point(18) = {58.5+x-b, 16.0-6.0-a, 0.0, n2};
Point(19) = {133.5-x+b, 24.0+6.0+a, 0.0, n2};
Point(20) = {133.5-x+b, 16.0-6.0-a, 0.0, n2};
Circle(9) = {9, 5, 17};
Circle(10) = {10, 6, 18};
Circle(11) = {20, 8, 12};
Circle(12) = {11, 7, 19};
Line(13) = {18, 20};
Line(14) = {12, 11};
Line(15) = {19, 17};
Line(16) = {9, 10};
Line Loop(17) = {7, 8, 1, 2, 3, 4, 5, 6};
Line Loop(18) = {9, -15, -12, -14, -11, -13, -10, -16};
Plane Surface(19) = {17, 18};
Physical Surface(1) = {19};
Physical Line(2) = {8};
Physical Line(3) = {4};
