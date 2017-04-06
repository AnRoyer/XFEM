Include "Line.dat";

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};

Line(1) = {1,2};
Transfinite Line {1} = 1000 Using Progression 1;

Physical Point(GAMMADIR) = {1};
Physical Point(GAMMAINF) = {2};
Physical Line(OMEGA) = {1};
