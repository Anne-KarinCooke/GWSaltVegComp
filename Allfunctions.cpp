#include <Rcpp.h>
using namespace Rcpp;
#include <array>
#include <list>
#include <iostream>
// #include "storageformats.hpp"
#include "Vegetationfunctions.h"
#include "Flux.h"
#include "Runoff.h"
#include "Infiltration.h"


// [[Rcpp::export]]
List balances2D( NumericVector Rain, int alpha_i=1,
                List vegpar, List soilpar, List saltpar,
                int rows=10, int cols=10,  double slope=0.001){

  // all these are part of your function you cannot define them outside
  int i;
  int j;
  int t;
  int tt;
  int t_old;

  int time=t.size();
  const int deltat=12;


  // soil moisture [mm]
  double M[rows][cols][time];

  //  plant biomass density
  double P[rows][cols][time];
  // //h
  double h[rows][cols][time];
  // //CM
  double CM[rows][cols][time];
  // // SmI
  double SmI[rows][cols][time];
  // // SmM
  double SmM[rows][cols][time];
  // // In
  double In[rows][cols][time];
  // //Svir
  double Svir[rows][cols][time];
  // // flux
  double flux[rows][cols][time];
  // //q
  double q[rows][cols][time];
  // //runon
  double runon[rows][cols][time];
  // //mb
  double mb[rows][cols][time];
  //
  // //
  // // Storage arrays for the SUBdaily time steps declared and initialized
  //
  // //h sub
  double h_sub[rows][cols][deltat];
  // //P sub
  double P_sub[rows][cols][deltat];
  // //M sub
  double M_sub[rows][cols][deltat];
  // //CM sub
  double CM_sub[rows][cols][deltat];
  // // SmIsub
  double SmI_sub[rows][cols][deltat];
  // // SmMsub
  double SmM_sub[rows][cols][deltat];
  // // Isub
  double I_sub[rows][cols][deltat];
  // //Svirsub
  double Svir_sub[rows][cols][deltat];
  // // fluxsub
  double flux_sub[rows][cols][deltat];
  // //qsubsub
  double q_sub[rows][cols][deltat];
  // //runonsub
  double runon_sub[rows][cols][deltat];
  // //mbsub
  double mb_sub[rows][cols][deltat];
  // //Gr_sub
  double Gr_sub[rows][cols][deltat];
  // //Mo_sub
  double Mo_sub[rows][cols][deltat];
  // //WU_sub
  double WU_sub[rows][cols][deltat];
  //
  // // Salt leaching
  double L_salt[rows][cols][deltat];
  // //salt rise
  double U_salt[rows][cols][deltat];
  //
  //
  // ///   OLDS

  // double h_old[rows][cols]={{0.0}};
  // double P_old[rows][cols]={{0.0}};
  // double M_old[rows][cols]={{0.0}};
  // double SmI_old[rows][cols]={{0.0}};
  // double CM_old[rows][cols]={{0.0}};
  // double Svir_old[rows][cols]={{0.0}};



  //// BIG OLD BALANCES FUNCTION


  double timeincr = 1/deltat;  // needs to be float and not double for some reason!

  // // double runon;
  // double alpha_i;
   double cn;
  double Mn;

  double rn[rows][cols];
  double Zras[rows][cols];



  for (i==1; i< rows; i++) {

    for (j==1; j< cols; j++ ){

      for (t==2; t< time; t++){

        for (tt==1; tt< (deltat); tt++){


          if(tt == 1) {
            int t_old = t-1;
          } else {
            int t_old = tt;
          }
          // calculation of sub daily runoff and runon
          runon_sub[i][j][tt] = rn[i][j]*q_sub[i][j][t_old];
          q_sub[i][j][tt] = OF(h[i][j][t_old], cn, Mn, slope[i][j])*timeincr;

          double Rain_in;

          if (Rain[t] > 0.0 & tt == 1){
            double Rain_in = 10.0*Rain[t];
          } else {
            double Rain_in = 0.0;
          }
          // calculate water depth on soil
          h_sub[i][j][tt] =  h[i][j][t_old] + Rain_in -
              Infil( h[i][j][t_old],  P[i][j][t_old], alpha_i, vegpar["k"], vegpar["W0"])-
              q_sub[i][j][tt] + runon_sub[i][j][tt];

          // adjust infiltration rate
          if(h_sub[i,j,tt]<(soilpar["K_s"]*timeincr)) {
            alpha_i = 1.0;
          } else {
            alpha_i = 1-(h_sub[i][j][tt]-soilpar["K_s"]*timeincr)/h_sub[i][j][tt];
          }

          I_sub[i][j][tt] = Infil(h[i][j][t_old], P[i][j][t_old],
                                  alpha_i, vegpar["k"], vegpar["W0"])*timeincr; //// problem like above, how to deal with parameters from list


          // Water uptake
          WU_sub[i][j][tt] = WU(M_sub[i][j][t_old],P[i][j][t_old], vegpar["gmax"], vegpar["k1"])*timeincr;
          // Growth
          Gr_sub[i][j][tt] = Gr(Svir[i][j][t_old], P[i][j][t_old], vegpar["c"], vegpar["gmax"], vegpar["k1"])*timeincr;
          //Mortality
          Mo_sub[i][j][tt] = Mo(P[i][j][t_old], M[i][j][t_old], Svir[i][j][t_old],vegpar["d"])*timeincr;
          // Plant biomass balance
          P_sub[i][j][tt] = P[i][j][t_old] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; /// not sure if this all is ok this way or too close to R

          // Water balance before drainage
          M_sub[i][j][tt] = M[i][j][t_old] + I_sub[i][j][tt] - WU_sub[i][j][tt];

          // Drainage/Capillary rise (vertical water flux)

          flux_sub[i][j][tt] = L_n(M_sub[i][j][tt],Zras[i][j],soilpar["b"],soilpar["K_s"],soilpar["psi_s_bar"]);  /// how yo read in Zras,,, change the soilpar and vegpar stuff

          // Adjustment for M including flux
          //
          M_sub[i][j][tt] = M_sub[i][j][tt] +  flux_sub[i][j][tt]*timeincr;

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
          Svir_sub[i][j][tt+1] = soilpar["n"]* vegpar["Zr"]*(pow((soilpar["h1bar"]* pow(10,-1)),(1/soilpar["b"])))*(soilpar["h1bar"]*pow(10,-1)*pow((M_sub[i][j][tt+1]/(soilpar["n"]*vegpar["Zr"])),-soilpar["b"]))+pow((3.6*CM_sub[i][j][tt+1]),(-1/soilpar["b"]));

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

        double sumI;
        double sumflux;
        double sumq;
        double sumrunon;
        double mb;

        for(int tt = 1; tt =< deltat; ++tt)
        {


          sumI += I_sub[i][j][tt];
          sumflux += flux_sub[i][j][tt]; // timeincr thingi?!
          sumq += q_sub[i][j][tt];
          sumrunon += runon_sub[i][j][tt];
          summb += mb_sub[i][j][tt];


        }


        In[i][j][t] = sumI;

        flux[i][j][t]= sumflux;

        q[i][j][t] = sumq;

        runon[i][j][t] = sumrunon;

        mb[i][j][t] = summb;

      }
    }
  }

 List out = Rcpp::List::create(Rcpp::Named("P") = P[rows][cols][time],
                                Rcpp::Named("M") = M[rows][cols][time],
                                Rcpp::Named("h") = h[rows][cols][time],
                                Rcpp::Named("CM") = CM[rows][cols][time],
                                Rcpp::Named("SmM") = SmM[rows][cols][time],
                                Rcpp::Named("In") = In[rows][cols][time],
                                Rcpp::Named("flux") = flux[rows][cols][time],
                                Rcpp::Named("Svir") = Svir[rows][cols][time],
                                Rcpp::Named("q") = q[rows][cols][time],
                                Rcpp::Named("mb") = mb[rows][cols][time],
                                Rcpp::Named("runon") = runon[rows][cols][time]);
  return out;
}


