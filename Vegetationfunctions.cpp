#include <Rcpp.h>
using namespace Rcpp;

// VEGETATION FUNCTIONS

//Plant water uptake function WU
// [[Rcpp::export]]
double WU(double M, double P, double gmax, double k1 ) {  /// not quite happy with the list item calling...list par
  
  double Wu = gmax*(M/(M+k1))*P;   // function is called WU, output is called Wu, cannot be the same
  return Wu;
  
}


//Plant Growth function Gr
// [[Rcpp::export]]
double Gr(double M, double P, double c, double gmax, double k1){
  
  double Gro = c*WU(M,P,gmax,k1);
  return Gro;
}


// Plant mortality function Mo
// [[Rcpp::export]]
double Mo(double P, double M, double Svir, double d ){
  
  double Mort=P*(d*(M/Svir));
  
  //List out = Rcpp::List::create(Rcpp::Named("Mort") = Mort);
  return Mort;
                         
}  

/*** R
WU(3.0,2.0,0.5,0.3)
Gr(3.0,2.0,0.3,0.5,0.3)
Mo(2.0,3.0,0.7,0.2)

# > WU(3.0,2.0,0.5,0.3)
#[1] 0.9090909

# > Gr(3.0,2.0,0.3,0.5,0.3)
# [1] 0.2727273

#> Mo(2.0,3.0,0.7,0.2)
#[1] 1.714286
*/
