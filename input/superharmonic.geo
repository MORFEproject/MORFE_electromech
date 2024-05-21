
// global dimensions

L=510; // length of beam
//L=200; // length of beam
H=1.5; // thick. of beams

// mesh parameters
 
//lc1=1.75;
lc1=5;
//lc1=20;

// point coordinates

Point(1) = {-L/2,0,0,lc1};
Point(2) = {0,0,0,lc1};
Point(3) = {L/2,0,0,lc1};
Point(4) = {L/2,H,0,lc1};
Point(5) = {-L/2,H,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};
Line Loop(6) = {1,2,3,4,5};

// surface

Plane Surface(1) = {6};

// physical entities

Physical Line(1) = {5}; // left clamp
Physical Line(2) = {3};  // right clamp
Physical Line(3) = {4};   // force
Physical Surface(4) = {1}; 

// creates second order mesh and saves

Mesh.MshFileVersion=2;
Mesh.ElementOrder=2;
Mesh 2;
Save 'demo2.msh';
