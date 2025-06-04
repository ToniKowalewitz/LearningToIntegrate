// Gmsh project created on Wed Aug 22 16:06:05 2018
// -> gmsh 3.0 or higher required.
SetFactory("OpenCASCADE");

// Create volumetric objects.
// These automatically are given as geometrical volume entities.
// -> Radii given by fourth entry.
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {0, 0, 0, 2, -Pi/2, Pi/2, 2*Pi};
// -> If radii are changed, make sure that box is larger than spheres.
Box(3) = {-2.1, -2.1, 0, 4.2, 4.2, 2.1};

// Use boolean operations to reduce everything to desired halfsphere problem.
BooleanDifference{ Volume{1}; Delete; }{ Volume{3}; }
BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; }
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; }

// Add RX point at [0,0,0].
// -> You can control mesh size by changing last entry.
Point(17) = {0, 0, 0, 1};
Point{17} In Surface{10};

// Add geometrical surface entities.
Line Loop(15) = {19};
Plane Surface(14) = {15};

Line Loop(16) = {22};
Plane Surface(15) = {15, 16};

// Add physical point entity.
Physical Point("RX") = {17};

// Add physical line entity.
Physical Line("TX") = {19};

// Add physical surface entities.
Physical Surface("surface") = {12, 10};
Physical Surface("subsurface") = {11};
Physical Surface("interior") = {9};

// Add physical volume entities.
Physical Volume("innerShell") = {1};
Physical Volume("outerShell") = {2};
