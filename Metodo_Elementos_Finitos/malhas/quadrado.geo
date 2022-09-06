
wall = 0.18;

L = 1.0;


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

Point(1)  = {0.0, 0.0, 0.0,  wall}; // p0
Point(2)  = {1*L, 0.0, 0.0,  wall}; // p3
Point(3)  = {1*L, 1*L, 0.0,  wall}; // p4
Point(4)  = {0.0, 1*L, 0.0,  wall}; // p5

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//+ boundary conditions for stream function
Physical Line('inlet') = {4};
Physical Line('outlet') = {2};
Physical Line('paredeInf') = {1};
Physical Line('paredeSup') = {3};

Physical Surface('surface') = {1};
