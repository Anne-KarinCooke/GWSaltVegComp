#include <Rcpp.h>
using namespace Rcpp;


//vertical water flux function (capillary rise and drainage), eq from Salvucci 1993
// [[Rcpp::export]]
double L_n(double M, double Z, double n, double Zr, double b, double hb,double K_s, double psi_s_bar){
  

  double s=M/(n*Zr); // extract n and Zr from list and define them
  double psi = hb*pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));// define hb and b
  double m = 2 + 3/b;
  double qf = (pow((Z/hb),(-m))-(pow((psi/hb),-m)))/((1+pow((psi/hb),-m))+(m-1)*pow((Z/hb),-m));
  double flux = K_s * qf;
  
  return flux;
                                                       
}


                                       
                                       
                                       /*** R
                                       L_n()
                                       
                                       */                                           