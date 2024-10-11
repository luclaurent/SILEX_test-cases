// load the CAD file
Merge "piece_caterpillar.step" ;

// give the elemental length
h=10;

// give the elemental length
Mesh.ElementOrder = 2;

Characteristic Length {30, 31, 32, 62, 40, 61, 33, 29, 41, 60, 59, 16, 7, 8, 42, 17, 2, 58, 6, 18, 27, 28, 57, 36, 14, 54, 38, 12, 63, 65, 64, 39, 37, 13, 66, 55, 11, 35, 56, 15, 50, 53, 26, 19, 3, 1, 4, 5, 43, 44, 23, 34, 20, 9, 10, 24, 22, 49, 45, 47, 25, 51, 21, 48, 52, 46} = h;

Physical Volume(1) = {1};
Physical Surface(2) = {32};
Physical Surface(3) = {22, 26};
Physical Surface(4) = {35, 36};
Physical Surface(5) = {18, 17};

