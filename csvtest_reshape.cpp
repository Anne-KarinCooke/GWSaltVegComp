#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]] 
mat write_elev(int rows, int cols) {
  // mat::fixed<3,9> flowdirTable;
  arma::mat elev(rows,cols, fill::randu);
  return elev;
  
}

// [[Rcpp::export]] 
List bla( int rows, int cols){
mat B = write_elev(rows, cols);
//Rcpp::Rcout << "B:" << std::endl << B << std::endl;
B.save("B.txt", arma::raw_ascii);
system("R CMD BATCH GeoTiff.R");
// load from disk

arma::mat flowdir_new;
flowdir_new.load("new_flowdir.txt");


arma::mat slp_matrix_new;
slp_matrix_new.load("slp_matrix.txt");


mat flowdir_new_ten(rows,cols, fill::zeros);
mat slp_matrix_new_ten(rows,cols, fill::zeros);


for(int i=0; i<(cols/2); i++){
  flowdir_new_ten(0,(cols/2)+i) = flowdir_new(1,0+i);
  flowdir_new_ten(0,0+i) = flowdir_new(0,0+i);
  
  for(int j=1; j<cols; j++){
  for(int k=1; (k<rows) & (k%2!=0) ; k++){
   
   for(int l=1; (l<rows) & (l%2==0) ; l++){
 



    flowdir_new_ten(j,(cols/2)+i) = flowdir_new(k,0+i); //1,3,5


    flowdir_new_ten(j,0+i) = flowdir_new(l,0+i); //2,4,6

       
    }
  }
 }
}



 
List output = List::create(Rcpp::Named("flowdir") =  flowdir_new_ten,
                           Rcpp::Named("slope") =  slp_matrix_new_ten);
                           
 

return output;
}



/*** R
rows=10
cols=10
results <- bla(rows, cols)

results$flowdir
results$slope

*/


