Merge "propeller.step";

h= 5;
Mesh.ElementOrder = 2;
// 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D)
Mesh.Algorithm3D=7;
// 2D mesh algorithm  (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
// Mesh.Algorithm=6;

Characteristic Length {36, 34, 37, 35, 27, 26, 25, 22, 33, 30, 43, 42, 41, 40, 57, 56, 55, 54, 17, 18, 29, 15, 16, 51, 50, 39, 23, 38, 45, 44, 28, 14, 67, 71, 70, 66, 68, 73, 72, 69, 24, 31, 65, 64, 52, 9, 53, 59, 58, 19, 21, 32, 20, 10, 49, 48, 47, 46, 63, 62, 61, 60, 7, 11, 12, 8, 13, 2, 3, 5, 6, 1, 4} = h;

Physical Surface(1) = {55, 56};
Physical Surface(2) = {31};
Physical Surface(3) = {27, 3};
Physical Volume(4) = {1};
Physical Surface(5) = {45, 39, 48, 49, 47, 43, 38, 44, 40, 42, 46, 36, 28, 41, 37, 53, 56, 55, 31, 50, 54, 51, 29, 20, 26, 18};



