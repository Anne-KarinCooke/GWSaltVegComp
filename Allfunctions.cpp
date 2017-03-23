#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
// [[Rcpp::export]]


int i;
int j;
int t;
int tt;

// int rain; /// thats no int

#include <array>
#include <list>


const int rows=10;
const int cols=10;
const int time=20;
const int deltat=12;

#include <array>

// const int rows=10;
// const int cols=10;
// const int time=20;
// const int deltat=12;


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

// #include "storageformats.hpp"
#include "Vegetationfunctions.hpp"
#include "Flux.hpp"
#include "Runoff.hpp"
#include "Infiltration.hpp"

//// BIG OLD BALANCES FUNCTION 

 // int rows;
 // int cols;
 // int time;
 // int deltat;

 
float timeincr = 1/deltat;  // needs to be float and not double for some reason!

 double slope;
// // double runon; 
// double alpha_i;
 double cn;
double Mn;
List Rain;

// 
 List soilpar;
 List vegpar;
List saltpar;
// 
double rn[rows][cols];
double Zras[rows][cols];

// [[Rcpp::export]]
List balances2D(){ 
  
  
  
  for (i==1; i< rows; i++) {
    
    for (j==1; j< cols; j++ ){
      
      for (t==2; t< time; t++){  
        
        for (tt==1; tt< (deltat-1); tt++){  
          
          
          if(tt == 1) {
            h_old[i][j] = h[i][j][t-1];
          } else {
            h_old[i][j] =  h_sub[i][j][tt];
          }
          if(tt == 1) {
            P_old[i][j] = P[i][j][t-1];
          } else {
            P_old[i][j] =  P_sub[i][j][tt];
          }
          if(tt == 1) {
            M_old[i][j]= M[i][j][t-1];
          } else {
            M_old[i][j] =  M_sub[i][j][tt];
          }
          if(tt == 1) {
            SmI_old[i][j]= SmI[i][j][t-1];
          } else {
            SmI_old[i][j] =  SmI_sub[i][j][tt];
          }
          if(tt == 1) {
            CM_old[i][j]= CM[i][j][t-1];
          } else {
            CM_old[i][j] =  CM_sub[i][j][tt];
          }
          if(tt == 1) {
            Svir_old[i][j]= Svir[i][j][t-1];
          } else {
            Svir_old[i][j] = Svir_sub[i][j][tt];
          }
          

          q_sub[i][j][tt+1] = OF(h_old[i][j], cn, Mn, slope[i][j])*timeincr; 
          
          runon_sub[i][j][tt+1] = rn[i][j]*q_sub[i][j][tt]; 
          
          
          if (tt==1){ // not right yet
            h_sub[i][j][tt+1] =  h_old[i][j] + (10*Rain[t]*timeincr) - Infil( h_old[i][j],  P_old[i][j], alpha_i, vegpar["k"], vegpar["W0"])- q_sub[i][j][tt] +runon_sub[i][j][tt];
            return h_sub[i][j][tt+1];
          }
          else {
            h_sub[i][j][tt+1] =  h_old[i][j]  - Infil( h_old[i][j],  P_old[i][j], alpha_i, vegpar["k"], vegpar["W0"])- q_sub[i][j][tt] +runon_sub[i][j][tt];
            return h_sub[i][j][tt+1];
          }
          
          
          if(h_sub[i,j,tt+1]<(soilpar["K_s"]*timeincr)) {
            
            alpha_i = 1.0;
          }
          else {
            
            alpha_i = 1-(h_sub[i][j][tt+1]-soilpar["K_s"]*timeincr)/h_sub[i][j][tt+1];
            
          }
          
          I_sub[i][j][tt+1] = Infil(h_old[i][j], P_old[i][j], alpha_i, vegpar["k"], vegpar["W0"])*timeincr; //// problem like above, how to deal with parameters from list
          
          
          // Water uptake           
          WU_sub[i][j][tt] = WU(M_sub[i][j][tt+1],P_old[i][j], vegpar["gmax"], vegpar["k1"])*timeincr;
          // Growth            
          Gr_sub[i][j][tt] = Gr(Svir_old[i][j], P_old[i][j], vegpar["c"], vegpar["gmax"], vegpar["k1"])*timeincr; 
          //Mortality
          Mo_sub[i][j][tt] = Mo(P_old[i][j], M_old[i][j], Svir_old[i][j],vegpar["d"])*timeincr;
          // Plant biomass balance             
          P_sub[i][j][tt+1] = P_old[i][j] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; /// not sure if this all is ok this way or too close to R
          
          // Water balance before drainage
          M_sub[i][j][tt+1] = M_old[i][j] + I_sub[i][j][tt] - WU_sub[i][j][tt];
          
          // Drainage/Capillary rise (vertical water flux)          
          
          flux_sub[i][j][tt+1] = L_n(M_sub[i][j][tt+1],Zras[i][j],soilpar["b"],soilpar["K_s"],soilpar["psi_s_bar"]);  /// how yo read in Zras,,, change the soilpar and vegpar stuff
          
          // Adjustment for M including flux
          //
          M_sub[i][j][tt+1] = M_sub[i][j][tt+1] +  flux_sub[i][j][tt+1]*timeincr; 
          
          // Salt balance
          
          // salt leaching
          // 
          if(flux_sub[i][j][tt+1] < 0.0 ) {
            L_salt[i][j][tt+1] = saltpar["f"] *CM_sub[i][j][tt+1] * flux_sub[i][j][tt+1]*timeincr;
          } else {
            L_salt[i][j][tt+1] = 0.0;
          }
          
          // salt upflow
          if(flux_sub[i][j][tt+1] > 0.0 ) {
            U_salt[i][j][tt+1] = CM_sub[i][j][tt+1] *saltpar["f"] * flux_sub[i][j][tt+1]*timeincr;
          } else {
            L_salt[i][j][tt+1] = 0.0;
          }
          
          
          // # salt mass coming in with infiltration
          SmI_sub[i][j][tt+1] = SmI_old[i][j] + I_sub[i][j][tt] * saltpar["ConcConst"];
          
          // #salt mass in soil
          SmM_sub[i][j][tt+1] = SmI_sub[i][j][tt+1] + U_salt[i][j][tt+1] - L_salt[i][j][tt+1];
          
          // # salt concentration in soil
          CM_sub[i][j][tt+1] = (SmM_sub[i][j][tt+1]/M_sub[i][j][tt+1])*(1/58.44);         
          
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
          Svir_sub[i][j][tt+1] = soilpar["n"]* vegpar["Zr"]*(pow((soilpar["h1bar"]* pow(10,-1)),(1/soilpar["b"])))*(soilpar["h1bar"]*pow(10,-1)*pow((M_sub[i][j][tt+1]/(soilpar["n"]*vegpar["Zr"])),-soilpar["b"]))+pow((3.6*CM_sub[i][j][tt+1]),(-1/soilpar["b"]))
            
            // in R:  Svir_sub[i,j,tt+1]<-soilpar$n*vegpar$Zr*((soilpar$h1bar*10^-1)^(1/soilpar$b))*
            // ((soilpar$h1bar*10^-1)*(M_sub[i,j,tt+1]/(soilpar$n*vegpar$Zr))^(-soilpar$b)+(3.6*CM_sub[i,j,tt+1]))^(-1/soilpar$b)
            //                         
            // # checking the mass balance
            mb_sub[i][j][tt] = I_sub[i][j][tt] - WU_sub[i][j][tt] + flux_sub[i][j][tt]*timeincr;  
          
        }
        // # Aggregating the substep results to daily values.
        
        P[i][j][t] = P_sub[i][j][deltat];
        M[i][j][t] = M_sub[i][j][deltat];
        h[i][j][t] = h_sub[i][j][deltat];
        CM[i][j][t] = CM_sub[i][j][deltat];
        SmI[i][j][t] = SmI_sub[i][j][deltat];
        SmM[i][j][t] = SmM_sub[i][j][deltat];
        Svir[i][j][t] = Svir_sub[i][j][deltat];
        
        
        for(int tt = 1; tt =< deltat; ++tt)
        {
          double sumI;
          double sumflux;
          double sumq;
          double sumrunon;
          double mb;
          
          sumI += I_sub[i][j][tt];
          sumflux += flux_sub[i][j][tt]; // timeincr thingi?!
          sumq += q_sub[i][j][tt];
          sumrunon += runon_sub[i][j][tt];
          summb += mb_sub[i][j][tt];
          
         // tt++;
        } 
        

        In[i][j][t] = sumI;

        flux[i][j][t]= sumflux; 

        q[i][j][t] = sumq;

        runon[i][j][t] = sumrunon;

        mb[i][j][t] = summb;
        
      }
    }
  }
 
 List out = Rcpp::List::create(Rcpp::Named("P") = P,
                                Rcpp::Named("M") = M, 
                                Rcpp::Named("h") = h, 
                                Rcpp::Named("CM") = CM, 
                                Rcpp::Named("SmM") = SmM, 
                                Rcpp::Named("In") = In, 
                                Rcpp::Named("flux") = flux, 
                                Rcpp::Named("Svir") = Svir, 
                                Rcpp::Named("q") = q,
                                Rcpp::Named("mb") = mb,
                                Rcpp::Named("runon") = runon);
  return out;      
}  


