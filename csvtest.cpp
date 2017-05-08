#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]] 
mat write_flowdirTable() {
  // mat::fixed<3,9> flowdirTable;
  arma::mat flowdirTable(10,10, fill::randu);
  return flowdirTable;
  
}

// [[Rcpp::export]] 
mat bla(){
mat B = write_flowdirTable();
//Rcpp::Rcout << "B:" << std::endl << B << std::endl;
B.save("B.txt", arma::raw_ascii);
system("R CMD BATCH GeoTiff.R");
// load from disk
arma::mat flowdir_new;
flowdir_new.load("new_flowdir.txt");
arma::mat slp_matrix_new;
slp_matrix_new.load("slp_matrix.txt");

return flowdir_new;
}


/*** R
bla()
*/


