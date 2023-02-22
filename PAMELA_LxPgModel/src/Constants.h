#define ProtonMass 938.2720813 // MeV/c^2
#define AU 149.5978707E9 // m [1 AU in meters]
#define r0 2.0E9 // [??]
#define f0 1.E-8 //

#define Vw 400   // km/s
#define HP 122   // AU
#define DRIFT 0  // drift: ON/OFF 
#define MULT 5

// #define BETAExp 1


// phi = PHIxKSCALE/kscale
#define PHIxKSCALE 2220 

// force local storage of output 
//#define STOREHERE = 0;

// ---- grids ----

const double RMIN = 0;
const double RMAX = HP*AU;
const int    NR = 610; // 122 x 5

const double PMIN= 10.;     // 10 MV
const double PMAX= 200.e+3; // 200 GV
const int    NP  = 500;

const double KMIN= 1.5;
const double KMAX=  22;
const double NK  = 100;

const double XMIN = -4.0;
const double XMAX =  4.0;
const double NX   =  99;





