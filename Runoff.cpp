#include <Rcpp.h>
using namespace Rcpp;
// Overland flow RUNOFF, kinematic wave approach, as used in from Saco et al 2013

double OF(double h, double cn, double Mn, double slope){
  
  double qq = (cn/Mn)*(pow(h,1.666667))*(pow(slope,0.5));
  return qq;
}

  
  