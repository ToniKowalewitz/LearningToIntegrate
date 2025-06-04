// Gmsh project created on Mon Aug 20 09:01:27 2018
SetFactory("OpenCASCADE");

// Define geometries.
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
Line Loop(1) = {1};
Plane Surface(1) = {1};

Circle(2) = {0, 0, 0, 2.5, 0, 2*Pi};
Line Loop(2) = {2};
Plane Surface(2) = {2};

Rectangle(3) = {-11, 0, 0, 22, 11, 0};

BooleanDifference{ Surface{1}; Delete;}{ Surface{3};}
BooleanDifference{ Surface{2}; Delete;}{ Surface{3}; Delete;}
BooleanDifference{ Surface{1}; Delete;}{ Surface{2};}

// Define physical entities.
Physical Line("surface") = {10, 11, 12};
Physical Line("subsurface") = {13};
Physical Line("interior") = {9};

Physical Surface("inner_part") = {2};
Physical Surface("outer_part") = {1};
