
r2 = 200;
r1 = 100;
lz = 750;

espesor = r2-r1 ;

ms = 1.75*espesor ; 

Point(1)  = {0,0,0,ms};

Point(2) = {r1,0,0,ms};
Point(3) = {0,-r1,0,ms};
Point(4) = {-r1,0,0,ms};
Point(5) = {0,r1,0,ms};

Point(6) = {r2,0,0,ms};
Point(7) = { 0,-r2,0,ms};
Point(8) = {-r2,0,0,ms};
Point(9) = {0,r2,0,ms};

Point(11)  = {0,0,lz,ms};

Point(12) = {r1,0,lz,ms};
Point(13) = { 0,-r1,lz,ms};
Point(14) = {-r1,0,lz,ms};
Point(15) = {0, r1,lz,ms};

Point(16) = {r2,0,lz,ms};
Point(17) = { 0,-r2,lz,ms};
Point(18) = {-r2,0,lz,ms};
Point(19) = {0,r2,lz,ms};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

Circle(11) = {12,11,13};
Circle(12) = {13,11,14};
Circle(13) = {14,11,15};
Circle(14) = {15,11,12};

Circle(15) = {16,11,17};
Circle(16) = {17,11,18};
Circle(17) = {18,11,19};
Circle(18) = {19,11,16};

Line(21) = {2,12};
Line(22) = {3,13};
Line(23) = {4,14};
Line(24) = {5,15};

Line(25) = {6,16};
Line(26) = {7,17};
Line(27) = {8,18};
Line(28) = {9,19};


Curve Loop(1) = {-5,25,15,-26};
Curve Loop(2) = {-6,26,16,-27};
Curve Loop(3) = {-7,27,17,-28};
Curve Loop(4) = {-8,28,18,-25};

Curve Loop(5) = {-1,21,11,-22};
Curve Loop(6) = {-2,22,12,-23};
Curve Loop(7) = {-3,23,13,-24};
Curve Loop(8) = {-4,24,14,-21};


Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};

Surface(5) = {-5};
Surface(6) = {-6};
Surface(7) = {-7};
Surface(8) = {-8};

Curve Loop(9) = {5,6,7,8};
Curve Loop(10) = {-4,-3,-2,-1};
//+
Plane Surface(9) = {9,10};

Curve Loop(11) = {-18,-17,-16,-15};
Curve Loop(12) = {11,12,13,14};
//+
Plane Surface(10) = {11,12};

Surface Loop(1) = {1,2,3,4,5,6,7,8,9,10};
Volume(1) = {1};


Physical Point("_node_fixed-ui") = {9,19};
Physical Point("_node_fixed-uj") = {6,16};

Physical Surface("_triangle_fixed-uk") = {9,10};
Physical Surface("_triangle_pressure") = {5,6,7,8};

Physical Volume("mat_tetrahedron_") = {1};
