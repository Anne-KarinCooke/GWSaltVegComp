#include <Rcpp.h>
using namespace Rcpp;

// Infiltration function Infil
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){

  double I=alpha_i*h*((P+k*W0)/(P+k));
  return I;
}

/*** R
Infil(3.0,2.0,0.5,0.3,0.7)
#> Infil(3.0,2.0,0.5,0.3,0.7)
#[1] 1.441304
*/


