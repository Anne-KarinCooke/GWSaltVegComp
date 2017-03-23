// Storage arrays for the daily time steps declared and initialized


#include <array>

const int rows=10;
const int cols=10;
const int time=20;
const int deltat=12;


// soil moisture [mm]
double M[rows][cols][time] = {{{0.0}}};

 //  plant biomass density
 double P[rows][cols][time]= {{{0.0}}};
// //h
double h[rows][cols][time]= {{{0.0}}};
// //CM
 double CM[rows][cols][time]= {{{0.0}}};
// // SmI
 double SmI[rows][cols][time]= {{{0.0}}};
// // SmM
double SmM[rows][cols][time]= {{{0.0}}};
// // In
 double In[rows][cols][time]= {{{0.0}}};
// //Svir
 double Svir[rows][cols][time]= {{{0.0}}};
// // flux
double flux[rows][cols][time]= {{{0.0}}};
// //q
double q[rows][cols][time]= {{{0.0}}};
// //runon
double runon[rows][cols][time]= {{{0.0}}};
// //mb
double mb[rows][cols][time]= {{{0.0}}};
//
// //
// // Storage arrays for the SUBdaily time steps declared and initialized
//
// //h sub
 double h_sub[rows][cols][deltat]= {{{0.0}}};
// //P sub
 double P_sub[rows][cols][deltat]= {{{0.0}}};
// //M sub
 double M_sub[rows][cols][deltat]= {{{0.0}}};
// //CM sub
 double CM_sub[rows][cols][deltat]= {{{0.0}}};
// // SmIsub
 double SmI_sub[rows][cols][deltat]= {{{0.0}}};
// // SmMsub
double SmM_sub[rows][cols][deltat]= {{{0.0}}};
// // Isub
 double I_sub[rows][cols][deltat]= {{{0.0}}};
// //Svirsub
double Svir_sub[rows][cols][deltat]= {{{0.0}}};
// // fluxsub
 double flux_sub[rows][cols][deltat]= {{{0.0}}};
// //qsubsub
 double q_sub[rows][cols][deltat]= {{{0.0}}};
// //runonsub
 double runon_sub[rows][cols][deltat]= {{{0.0}}};
// //mbsub
 double mb_sub[rows][cols][deltat]= {{{0.0}}};
// //Gr_sub
double Gr_sub[rows][cols][deltat]= {{{0.0}}};
// //Mo_sub
 double Mo_sub[rows][cols][deltat]= {{{0.0}}};
// //WU_sub
 double WU_sub[rows][cols][deltat]= {{{0.0}}};
//
// // Salt leaching
double L_salt[rows][cols][deltat]= {{{0.0}}};
// //salt rise
double U_salt[rows][cols][deltat]={{{0.0}}};
//
//
// ///   OLDS
 double h_old[rows][cols]={{0.0}};
 double P_old[rows][cols]={{0.0}};
 double M_old[rows][cols]={{0.0}};
 double SmI_old[rows][cols]={{0.0}};
 double CM_old[rows][cols]={{0.0}};
 double Svir_old[rows][cols]={{0.0}};
