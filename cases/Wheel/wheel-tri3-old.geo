n = 4;
h = 5;
R1 = 30;
R2 = 90;
R3 = 130;
e  = 20;
beta = 2*Pi/n;
Point(1) = {R1, 0, 0, h};
Point(2) = {R2, 0, 0, h};
Point(3) = {R3, 0, 0, h};
Point(4) = {R1*Cos(Asin(e/(2*R1))), e/2, 0, h};
Point(5) = {R2*Cos(Asin(e/(2*R2))), e/2, 0, h};
Point(6) = {R1*Cos(beta/2), R1*Sin(beta/2), 0, h};
Point(7) = {R2*Cos(beta/2), R2*Sin(beta/2), 0, h};
Point(8) = {R3*Cos(beta/2), R3*Sin(beta/2), 0, h};
Point(9) = {0, 0, 0, h};
Line(1) = {4, 5};
Circle(2) = {4, 9, 6};
Circle(3) = {5, 9, 7};
Circle(4) = {3, 9, 8};
Symmetry {0, 1, 0, 0} {
  Duplicata { Line{2, 1, 3, 4}; }
}
For i In {1:(n-1)}
Rotate {{0, 0, 1}, {0, 0, 0}, beta*i} {
  Duplicata { Line{5, 2, 7, 6, 1, 3, 8, 4}; }
}
EndFor

If (n == 2)
Line Loop(17) = {15, -4, 8, -16};
Line Loop(18) = {11, -3, -1, 2, -9, 12};
Line Loop(19) = {14, -7, -6, 5, -10, 13};
Plane Surface(100) = {17, 18, 19};
EndIf

If (n == 3)
Line Loop(25) = {16, -23, 24, -8, 4, -15};
Line Loop(26) = {14, -19, -20, 17, -10, 13};
Line Loop(27) = {12, 11, -3, -1, 2, -9};
Line Loop(28) = {6, 7, -22, -21, 18, -5};
Plane Surface(100) = {25, 26, 27, 28};
EndIf

If (n == 4)
Line Loop(33) = {16, -23, 24, -31, 32, -8, 4, -15};
Line Loop(34) = {14, -19, -20, 17, -10, 13};
Line Loop(35) = {3, -11, -12, 9, -2, 1};
Line Loop(36) = {6, 7, -30, -29, 26, -5};
Line Loop(37) = {25, -18, 21, 22, -27, -28};
Plane Surface(100) = {33, 34, 35, 36, 37};
EndIf

If (n == 6)
Line Loop(49) = {24, -31, 32, -39, 40, -47, 48, -8, 4, -15, 16, -23};
Line Loop(50) = {22, -27, -28, 25, -18, 21};
Line Loop(51) = {20, 19, -14, -13, 10, -17};
Line Loop(52) = {12, 11, -3, -1, 2, -9};
Line Loop(53) = {6, 7, -46, -45, 42, -5};
Line Loop(54) = {44, 43, -38, -37, 34, -41};
Line Loop(55) = {36, 35, -30, -29, 26, -33};
Plane Surface(100) = {49, 50, 51, 52, 53, 54, 55};
EndIf

If (n == 7)
Line Loop(57) = {31, -24, 23, -16, 15, -4, 8, -56, 55, -48, 47, -40, 39, -32};
Line Loop(58) = {29, 30, -35, -36, 33, -26};
Line Loop(59) = {28, 27, -22, -21, 18, -25};
Line Loop(60) = {20, 19, -14, -13, 10, -17};
Line Loop(61) = {12, 11, -3, -1, 2, -9};
Line Loop(62) = {6, 7, -54, -53, 50, -5};
Line Loop(63) = {52, 51, -46, -45, 42, -49};
Line Loop(64) = {44, 43, -38, -37, 34, -41};
Plane Surface(100) = {57, 58, 59, 60, 61, 62, 63, 64};
EndIf








