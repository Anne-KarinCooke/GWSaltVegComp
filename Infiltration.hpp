#include <Rcpp.h>
using namespace Rcpp;

// Infiltration function Infil

double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*(P+k*W0)/(P+k);
  return I;
}

