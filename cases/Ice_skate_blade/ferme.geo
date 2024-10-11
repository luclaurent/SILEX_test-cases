// ferme
h = 0.3;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.75, 0.0, 0.0, h};
Point(3) = {3.5, 0.0, 0.0, h};
Point(4) = {5.25, 0.0, 0.0, h};
Point(5) = {7.0, 0.0, 0.0, h};
Point(6) = {1.75, 0.75, 0.0, h};
Point(7) = {3.5, 1.5, 0.0, h};
Point(8) = {5.25, 0.75, 0.0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {1, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {2, 6};
Line(10) = {4, 8};
Line(11) = {3, 6};
Line(12) = {3, 7};
Line(13) = {3, 8};

Physical Line(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13};
