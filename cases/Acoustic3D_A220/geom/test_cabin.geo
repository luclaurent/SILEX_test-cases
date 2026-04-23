//------------------------------------------------------------
// Airbus A220 cabin air volume (no seats) - approximate model
// Units: meters
// Coordinate system:
//   x: longitudinal (nose -> tail)
//   y: lateral (left/right)
//   z: vertical (floor -> ceiling)
//------------------------------------------------------------

SetFactory("OpenCASCADE");

//----------------------
// Global mesh size
//----------------------
lc = 0.15; // characteristic length for points (adjust as needed)

//----------------------
// Basic unit conversion
//----------------------
inch = 0.0254;

//----------------------
// Cabin longitudinal dimensions (from provided drawing)
//----------------------
// Forward area (cockpit to first seat row)
L_front  = 32.0 * inch;

// Forward cabin: 6 rows @ 33 in pitch
L_block1 = 6.0 * 33.0 * inch;

// Door / galley area between cabins
L_door   = 42.0 * inch;

// Aft cabin: 14 rows @ 32 in pitch
L_block2 = 14.0 * 32.0 * inch;

// Total cabin length (can be overridden if you have a more precise value)
L_cabin  = L_front + L_block1 + L_door + L_block2;

//----------------------
// Cross-section dimensions (interior cabin)
//----------------------
// Width at floor level from drawing: 121.9 in (3.10 m)
width_floor_in  = 121.9;
width_floor     = width_floor_in * inch;
halfWidthFloor  = 0.5 * width_floor;

// Cabin interior height from floor to crown: 84 in (2.13 m)
height_cabin_in = 84.0;
height_cabin    = height_cabin_in * inch;

//----------------------
// Circle geometry for inner fuselage
// We construct a circle that passes through:
//   (y, z) = (-halfWidthFloor, 0)
//   (y, z) = ( halfWidthFloor, 0)
//   (y, z) = (0, height_cabin)
// The circle centre (0, zc) and radius R are solved analytically:
//   zc = (h^2 - a^2) / (2*h)
//   R  = sqrt(a^2 + zc^2)
// where h = height_cabin, a = halfWidthFloor.
//----------------------
a  = halfWidthFloor;
h  = height_cabin;
zc = (h*h - a*a) / (2.0*h);
R  = Sqrt(a*a + zc*zc);

//----------------------
// Points (cross-section in the y-z plane at x = 0)
//----------------------
Point(1) = {0.0, -halfWidthFloor, 0.0,          lc}; // left floor corner
Point(2) = {0.0,  halfWidthFloor, 0.0,          lc}; // right floor corner
Point(3) = {0.0,  0.0,            height_cabin, lc}; // ceiling crown
Point(4) = {0.0,  0.0,            zc,           lc}; // circle centre (for arcs)

//----------------------
// Curves: circular fuselage + flat floor
//----------------------
// Left side arc: from left floor corner up to crown
Circle(1) = {1, 4, 3};
// Right side arc: from crown down to right floor corner
Circle(2) = {3, 4, 2};
// Floor between right and left corners
Line(3)   = {2, 1};

// Closed curve loop and surface (cabin cross-section)
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

//----------------------
// Extrude cross-section along x to form 3D cabin air volume
//----------------------
extr[] = Extrude{L_cabin, 0.0, 0.0} {
  Surface{1};
};

// The volume created by the extrusion is extr[1] when using OpenCASCADE
cabinVolume = extr[1];

//----------------------
// Physical groups for use in solvers
//----------------------
// Volume of cabin air
Physical Volume("CabinAir") = {cabinVolume};

// Optionally tag inlet/outlet/bulkhead surfaces for boundary conditions:
frontSurface = extr[0];        // original surface moved to x = L_cabin
backSurface  = 1;              // original cross-section at x = 0 (by construction)

// Fuselage wall (cylindrical shell) surfaces created by the extrusion
// are in extr[2], extr[3], extr[4] ... but numbering may vary;
// you can comment these lines if you prefer to inspect in Gmsh first.
Physical Surface("FrontBulkhead") = {frontSurface};
Physical Surface("BackBulkhead")  = {backSurface};

// You may also group all lateral wall surfaces like this after inspection:
// Physical Surface("CabinWall") = {extr[2], extr[3], extr[4]};