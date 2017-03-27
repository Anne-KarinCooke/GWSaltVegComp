#include <Rcpp.h>
using namespace Rcpp;
 // [[Rcpp::export]]
 double OF(double h, double cn, double Mn, double slope){
   
   double qq = (cn/Mn)*(pow(h,1.666667))*sqrt(slope);
   return qq;
 }
                                    

             
                       
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*(P+k*W0)/(P+k);
  return I;
  
  
}
// [[Rcpp::export]]
List SurfaceWB(double alpha_i=1, double cn = 4, double Mn =10, double Rain =2){ // Surface Balance 
 
  
  int i;
  int j;
  int t;
  int tt;
  int t_old;
  const int deltat = 12;
  
  for (tt=1; tt< deltat; tt++){ 
    
    
  const int rows = 10;
  const int cols = 10;

  const int time = 10;
  double slope = 0.001;
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  List vegpar = (Rcpp::List::create(Rcpp::Named("k") = k,
                 Rcpp::List::create(Rcpp::Named("W0") = W0)));
  double K_s = 3.51; // cm/day
  List soilpar= Rcpp::List::create(Rcpp::Named("K_s") = K_s);

  
  // //h sub
  double h_sub[rows][cols][deltat] ={{{0}}};

  // //qsub
  double q_sub[rows][cols][deltat]={{{0}}};
  // //runonsub
  double runon_sub[rows][cols][deltat]={{{0}}};
  //  Isub
  double I_sub[rows][cols][deltat]={{{0}}};
  
  double h[rows][cols][time]={{{0}}};
  double P[rows][cols][time]={{{0}}};
  
  double timeincr = 1/12;

  double rn[rows][cols];

    
  if(tt == 1) {
    int t_old = t-1;
  } else {
    int t_old = tt;
  }
  // calculation of sub daily runoff and runon
  runon_sub[i][j][tt] = rn[i][j]*q_sub[i][j][t_old]; 
  q_sub[i][j][tt] = OF(h[i][j][t_old], cn, Mn, slope)*timeincr; 
  
  double Rain_in;
    
  if (Rain > 0.0 & tt == 1){ 
    Rain_in = 10.0*Rain; 
  } else {
   Rain_in = 0.0;
  }
  // calculate water depth on soil
  h_sub[i][j][tt] =  h[i][j][t_old] + Rain_in - 
    Infil( h[i][j][t_old],  P[i][j][t_old], alpha_i, vegpar["k"], vegpar["W0"])- 
    q_sub[i][j][tt] + runon_sub[i][j][tt];
  
  // adjust infiltration rate
  if(h_sub[i][j][tt] < (K_s*timeincr)) {
    alpha_i = 1.0;
  } else {
    alpha_i = 1-(h_sub[i][j][tt] - (K_s*timeincr))/h_sub[i][j][tt];
  }
  
  I_sub[i][j][tt] = Infil(h[i][j][t_old], P[i][j][t_old], 
                          alpha_i, vegpar["k"], vegpar["W0"])*timeincr;
  

}
  return(Rcpp::List::create(Rcpp::Named("I_sub") = I_sub,
                            Rcpp::Named("h_sub") = h_sub));
}

/*** R

*/
