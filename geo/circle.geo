// Gmsh project created on Sun Mar 26 11:44:28 2017
meshSize = 10;
Include "line.dat";
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {-1, 0, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Point(6) = {8, 0, 0, 1.0};
//+
Point(7) = {0, 8, 0, 1.0};
//+
Point(8) = {-8, 0, 0, 1.0};
//+
Point(9) = {0, -8, 0, 1.0};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {6, 1, 7};
//+
Circle(6) = {7, 1, 8};
//+
Circle(7) = {8, 1, 9};
//+
Circle(8) = {9, 1, 6};
//+
Line Loop(9) = {6, 7, 8, 5};
//+
Line Loop(10) = {2, 3, 4, 1};
//+
Plane Surface(11) = {9, 10};
//+
Physical Line(GAMMADIR) = {2, 3, 4, 1};
//+
Physical Line(GAMMAINF) = {6, 7, 5, 8};
//+
Physical Surface(OMEGA) = {11};
//+
Transfinite Line {2, 1, 3, 4} = meshSize Using Progression 1;
//+
Transfinite Line {6, 5, 8, 7} = 8*meshSize Using Progression 1;
