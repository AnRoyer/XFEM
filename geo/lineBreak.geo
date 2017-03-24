Include "Line.dat";

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {2, 1, 0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Transfinite Line {1} = 100 Using Progression 1;
Transfinite Line {2} = 100*Sqrt[2] Using Progression 1;

Physical Point(GAMMADIR) = {1};
Physical Point(GAMMAINF) = {3};
Physical Line(OMEGA) = {1,2};
