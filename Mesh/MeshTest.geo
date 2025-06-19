// Gmsh project created on Sat Jun 07 09:45:24 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.45, 0, 0, 1.0};
//+
Point(6) = {0.55, 0, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("fixed", 1) = {4};
//+
Transfinite Curve {6} = 100 Using Progression 1;
//+
Transfinite Curve {1} = 100 Using Progression 1;
//+
Transfinite Curve {2} = 100 Using Progression 1;
//+
Transfinite Curve {5} = 45 Using Progression 1;
//+
Transfinite Curve {3} = 45 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;