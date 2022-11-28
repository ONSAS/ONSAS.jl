

// sizes for another case: rInt = 0.05 ; rExt = 0.15 ;
rInt = 0.10 ; rExt = 0.15 ; tz = -0.2 ;

msExt = 0.012*(rExt*2.0*3.14) ;  //
msInt = 0.012*(rInt*2.0*3.14) ;  //

// plane z=0
Point(1)  = {0   , 0, 0};

Point(2)  = { rInt,     0, 0, msInt};
Point(3)  = {    0, -rInt, 0, msInt};
Point(4)  = {-rInt,     0, 0, msInt};
Point(5)  = {    0,  rInt, 0, msInt};

Point(7)  = { rExt,     0, 0, msExt};
Point(8)  = {    0, -rExt, 0, msExt};
Point(9)  = {-rExt,     0, 0, msExt};
Point(10) = {   0,  rExt, 0, msExt};

// plane z=tz
Point(21) = {0   , 0, tz};

Point(22) = { rInt,     0, tz, msInt};
Point(23) = {    0, -rInt, tz, msInt};
Point(24) = {-rInt,     0, tz, msInt};
Point(25) = {    0,  rInt, tz, msInt};

Point(27) = { rExt,     0, tz, msExt};
Point(28) = {    0, -rExt, tz, msExt};
Point(29) = {-rExt,     0, tz, msExt};
Point(30) = {    0,  rExt, tz, msExt};


// plane z=0
Circle(1) = {2, 1, 3}; //
Circle(2) = {3, 1, 4}; //
Circle(3) = {4, 1, 5}; //
Circle(4) = {5, 1, 2}; //

Circle(5) = {7, 1, 10}; //
Circle(6) = {10, 1, 9}; //
Circle(7) = {9, 1, 8}; //
Circle(8) = {8, 1, 7}; //

// plane z=tz
Circle(11) = {22, 21, 23}; //
Circle(12) = {23, 21, 24}; //
Circle(13) = {24, 21, 25}; //
Circle(14) = {25, 21, 22}; //

Circle(15) = {27, 21, 30}; //
Circle(16) = {30, 21, 29}; //
Circle(17) = {29, 21, 28}; //
Circle(18) = {28, 21, 27}; //

// inner lines
Line(21) = {2, 22}; //
Line(22) = {3, 23}; //
Line(23) = {4, 24}; //
Line(24) = {5, 25}; //

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

Line Loop(2) = {11, 12, 13, 14, 15, 16, 17, 18};
Plane Surface(2) = {2};

// inner surface
Line Loop(3) = {11, -22, -1, 21 };
Surface(3) = {3};

Line Loop(4) = {12, -23, -2, 22 };
Surface(4) = {4};

Line Loop(5) = {13, -24, -3, 23 };
Surface(5) = {5};

Line Loop(6) = {14, -21, -4, 24 };
Surface(6) = {6};

// outer lines
Line(25) = {7,  27}; //
Line(26) = {8,  28}; //
Line(27) = {9,  29}; //
Line(28) = {10, 30}; //

// outer surface
Line Loop(7) = {18, -25, -8, 26 };
Surface(7) = {7};

Line Loop(8) = {17, -26, -7, 27 };
Surface(8) = {8};

Line Loop(9) = {16, -27, -6, 28 };
Surface(9) = {9};

Line Loop(10) = {15, -28, -5, 25 };
Surface(10) = {10};


Surface Loop(1) = {1,2,3,4,5,6,7,8,9,10};
Volume(1) = {1};

Physical Volume("ring_solid") = {1};
