// load the CAD file
Merge "piston_sans_segments.step" ;

// give the elemental length
h=0.8;
Characteristic Length {10, 5, 13, 20, 1, 2, 3, 4, 12, 6, 11, 15, 14, 17, 16, 21, 24, 7, 26, 8, 25, 23, 22, 18, 9, 19} = h;

// physical group = 1 : symetry surface, y-fixed
Physical Surface(1) = {15, 13};

// physical group = 2 : cylinder surface
Physical Surface(2) = {10};

// physical group = 3 : symetry surface, x-fixed
Physical Surface(3) = {14};

// physical group = 4 : "pressure surface"
Physical Surface(4) = {4, 3};

// physical group = 5 : volume elements
Physical Volume(5) = {1};

