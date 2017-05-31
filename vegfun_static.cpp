#include <Rcpp.h>
using namespace Rcpp;


// very close to soilfun.cpp, but for veg

// ########################################
// [[Rcpp::export]]
List  Veg_cpp() {
  
// default

   double k = 12.0;//Saco et al, 2013
   double W0 = 0.2;//Saco et al, 2013
   double gmax = 0.05;//Saco et al, 2013
   double c = 10.0;//Saco et al, 2013
   double k1 = 5.0;//Saco et al, 2013
   double d = 0.24;//Saco et al, 2013 //fraction of plant mortality



 
 return(Rcpp::List::create(Rcpp::Named("k") = k,
                            Rcpp::Named("W0") = W0,
                            Rcpp::Named("gmax") = gmax,
                            Rcpp::Named("c") = c,
                            Rcpp::Named("k1") = k1,
                            Rcpp::Named("d") = d));

}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

//  /*** R
// Soil_cpp("S Clay Loam")
// */