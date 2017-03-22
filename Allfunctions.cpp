#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
//Sandy Clay Loam
// Rain_function<-function(time){ Rain <- rep(1, time)
//   return(Rain)}
// Rain<-Rain_function(time=time)

#include <list>

int const rows = 10;
int const cols = 10;
int const time = 10;
int const deltat =12;

float timeincr = 1/deltat;  // needs to be float and not double for some reason!

int i;
int j;
int t;
int tt;

// int rain; /// thats no int

#include <array>
// #include <iostream>

#include "storageformats.hpp"
#include "Constants.hpp"
#include "Vegetationfunctions.hpp"

///array matrix, slope matrix

//// BIG OLD BALANCES FUNCTION 

double balances2D(rows, cols, time, deltat, List rain, List soilpar, List vegpar, List saltpar, slope, runon){ // fix rain and slope
  
  
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
          
          double slope; // think about this again
          
          q_sub[i][j][tt+1] = OF(h_old[i][j], soilpar["cn"], soilpar["Mn"], slope)*timeincr; // define slope!!
          
          runon_sub[i][j][tt+1] = rn[i][j]*q_sub[i][j][tt]; /// rn !!! source
          
          
          if (tt==1){ // not right yet
            h_sub[i][j][tt+1] =  h_old[i][j] + (10*rain[t]*timeincr) - Infil( h_old[i][j],  P_old[i][j], soilpar["alpha_i"], vegpar["k"], vegpar["W0"])- q_sub[i][j][tt] +runon_sub[i][j][tt];
            return h_sub[i][j][tt+1];
          }
          else {
            h_sub[i][j][tt+1] =  h_old[i][j]  - Infil( h_old[i][j],  P_old[i][j], soilpar["alpha_i"], vegpar["k"], vegpar["W0"])- q_sub[i][j][tt] +runon_sub[i][j][tt];
            return h_sub[i][j][tt+1];
          }
          
          
          if(h_sub[i,j,tt+1]<(soilpar["K_s"]*timeincr)) {
            
            soilpar["alpha_i"] = 1.0;
          }
          else {
            
            soilpar["alpha_i"] = 1-(h_sub[i][j][tt+1]-soilpar["K_s"]*timeincr)/h_sub[i][j][tt+1];
            
          }
          
          I_sub[i][j][tt+1] = Infil(h_old[i][j], P_old[i][j], soilpar["alpha_i"], vegpar["k"], vegpar["W0"])*timeincr; //// problem like above, how to deal with parameters from list
          
          
          // Water uptake           
          WU_sub[i][j][tt] = WU(M_sub[i][j][tt+1],P_old[i][j],soilpar["gmax"], soilpar["k1"])*timeincr;
          // Growth            
          Gr_sub[i][j][tt] = Gr(Svir_old[i][j], P_old[i][j], soilpar["c"], soilpar["gmax"], soilpar["k1"])*timeincr; 
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
            L_salt[i][j][tt+1] = soilpar["f"] *CM_sub[i][j][tt+1] * flux_sub[i][j][tt+1]*timeincr;
          } else {
            L_salt[i][j][tt+1] = 0.0;
          }
          
          // salt upflow
          if(flux_sub[i][j][tt+1] > 0.0 ) {
            U_salt[i][j][tt+1] = CM_sub[i][j][tt+1] *soilpar["f"] * flux_sub[i][j][tt+1]*timeincr;
          } else {
            L_salt[i][j][tt+1] = 0.0;
          }
          
          
          // # salt mass coming in with infiltration
          SmI_sub[i][j][tt+1] = SmI_old[i][j] + I_sub[i][j][tt] * soilpar["ConcConst"];
          
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
        
        
        double sumI = std::accumulate(I_sub[i][j][1], I_sub[i][j][deltat], 0);
        In[i][j][t] = sumI;
        
        double sumflux = std::accumulate(flux_sub[i][j][1], flux_sub[i][j][deltat], 0);
        flux[i][j][t]= sumflux; 
        
        double sumq = std::accumulate(q_sub[i][j][1], q_sub[i][j][deltat], 0);
        q[i][j][t] = sumq;
        
        double sumrunon = std::accumulate(runon_sub[i][j][1], runon_sub[i][j][deltat], 0);
        runon[i][j][t] = sumq;
        
        double summb = std::accumulate(mb_sub[i][j][1], mb_sub[i][j][deltat], 0);
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
