#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SurfaceWB(const int rows=10, const int cols=10, int i,int j, int tt, int t_old, double alpha_i, double timeincr, double cn, double Mn, double Rain){ // Surface Balance 
  
  // //h sub
  double h_sub[rows][cols][deltat] ={{{0}}};
  // //P sub
  double P_sub[rows][cols][deltat]={{{0}}};
  // //M sub
  double M_sub[rows][cols][deltat]={{{0}}};
  // //qsub
  double q_sub[rows][cols][deltat]={{{0}}};
  // //runonsub
  double runon_sub[rows][cols][deltat]={{{0}}};
  
  
  if(tt == 1) {
    int t_old = t-1;
  } else {
    int t_old = tt;
  }
  // calculation of sub daily runoff and runon
  runon_sub[i][j][tt] = rn[i][j]*q_sub[i][j][t_old]; 
  q_sub[i][j][tt] = OF(h[i][j][t_old], cn, Mn, slope[i][j])*timeincr; 
  
  
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
                          alpha_i, vegpar["k"], vegpar["W0"])*timeincr;
}



/*** R

*/
