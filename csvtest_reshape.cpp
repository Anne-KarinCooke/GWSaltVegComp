#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]] 
mat write_elev() {
  // mat::fixed<3,9> flowdirTable;
  arma::mat elev(10,10, fill::randu);
  return elev;
  
}

// [[Rcpp::export]] 
List bla(){
mat B = write_elev();
//Rcpp::Rcout << "B:" << std::endl << B << std::endl;
B.save("B.txt", arma::raw_ascii);
system("R CMD BATCH GeoTiff.R");
// load from disk
arma::mat flowdir_new;
flowdir_new.load("new_flowdir.txt");
arma::mat slp_matrix_new;



slp_matrix_new.load("slp_matrix.txt");
mat slp_matrix_new_ten = reshape(slp_matrix_new, 10, 10);

mat flowdir_new_ten(10,10, fill::zeros);

vec vec1 = flowdir_new.row(1);
vec vec2 = flowdir_new.row(3);
vec vec3 = flowdir_new.row(5);
vec vec4 = flowdir_new.row(7);
vec vec5 = flowdir_new.row(9);
vec vec6 = flowdir_new.row(11);
vec vec7 = flowdir_new.row(13);
vec vec8 = flowdir_new.row(15);
vec vec9 = flowdir_new.row(17);
vec vec10 = flowdir_new.row(19);


  flowdir_new_ten(0,6::10)
    
flowdir_new.row(1) = 

for(int i=0; i<5; i++){
  
  for(int j=1; j<17; j++){

    
    flowdir_new_ten(0,6+i) = flowdir_new(1,0+i);
    
    flowdir_new_ten(1,6+i) = flowdir_new(0,1+i);
    
    if ( j % 2 == 0){
    flowdir_new_ten(j,6+i) = flowdir_new(j+1,1+i);
    }
 }
}

 // extract a column vector


 
List output = List::create(Rcpp::Named("slp_matrix_new_ten") = slp_matrix_new_ten,
                           Rcpp::Named("flowdir_new_ten") =  flowdir_new_ten);
 

return output;
}



/*** R
bla()

*/


