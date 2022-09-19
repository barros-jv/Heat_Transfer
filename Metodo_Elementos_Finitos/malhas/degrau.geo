
wall = 0.18;

L1 = 1.0;
L2 = 10.0;

/* 
 *        4                4L                  3
 *         o -------------------------------- o         
 *         |                                  |       
 *         |                                  |       
 *         |                                  |       Y         
 *      L  |                                  | L     ^
 *         |                                  |       |
 *         |                                  |       |
 *         o -------------------------------- o       o -----> X
 *        1                 4L                 2
 * */

Point(1)  = {0.0, 0.5*L1, 0.0,  wall};
Point(2)  = {0.1*L2, 0.5*L1, 0.0,  wall};
Point(3)  = {0.1*L2, 0.0, 0.0,  wall};
Point(4)  = {1*L2, 0.0, 0.0,  wall}; 
Point(5)  = {1*L2, 1*L1, 0.0,  wall};
Point(6)  = {0.0, 1*L1, 0.0,  wall};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

//+
Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

//+ boundary conditions for stream function
Physical Line('inlet') = {6};
Physical Line('outlet') = {4};
Physical Line('paredeInf') = {1,2,3};
Physical Line('paredeSup') = {5};

Physical Surface('surface') = {1};
