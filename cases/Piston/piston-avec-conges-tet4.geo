// load the CAD file
Merge "piston_avec_conges.step" ;

// give the elemental length


// coefficient raffinement maillage

// taille generale
h=2;
Characteristic Length {56, 53, 52, 55, 11, 43, 42, 48, 49, 41, 51, 10, 31, 30, 36, 44, 45, 46, 14, 37, 29, 39, 47, 50, 40, 19, 18, 24, 32, 33, 17, 34, 25, 6, 27, 35, 38, 28, 20, 21, 22, 23, 26, 5, 13, 57, 7, 58, 15, 16, 4, 2, 1, 3, 54, 12, 60, 62, 61, 59, 9, 8} = h;

// dans les conges, direction radiale
n1=4/(h/5);
Transfinite Line {32, 88, 28, 89, 22, 35, 46, 90, 42, 91, 36, 60, 49, 92, 56, 93, 50, 30, 87, 26, 86, 24, 33, 44, 85, 40, 84, 38, 47, 58, 83, 54, 82, 52} = n1 Using Progression 1;

// dans les conges, direction circonferentielle
n2=80/(h/4);
Transfinite Line {6, 25, 34, 39, 48, 53, 31, 29, 23, 45, 43, 37, 59, 57, 51, 27, 41, 55} = n2 Using Progression 1;

// physical group = 1 : symetry surface, y-fixed
Physical Surface(1) = {25, 27};

// physical group = 2 : cylinder surface
Physical Surface(2) = {22};

// physical group = 3 : symetry surface, x-fixed
Physical Surface(3) = {26};

// physical group = 4 : "pressure surface"
Physical Surface(4) = {17, 18};

// physical group = 5 : volume elements
Physical Volume(5) = {1};


