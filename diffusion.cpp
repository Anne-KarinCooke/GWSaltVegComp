#include <cmath> 
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]] 
mat write_DiffdirectionTable() {

  arma::mat DiffdirTable(2,8, fill::zeros);
  
  DiffdirTable(0,0) = 0.0;
  DiffdirTable(0,1) = -1.0;
  DiffdirTable(0,2) = -1.0;
  DiffdirTable(0,3) = -1.0;
  DiffdirTable(0,4) = 0.0;
  DiffdirTable(0,5) = 1.0;
  DiffdirTable(0,6) = 1.0;
  DiffdirTable(0,7) = 1.0;

  
  DiffdirTable(1,0) = 1.0;
  DiffdirTable(1,1) = 1.0;
  DiffdirTable(1,2) = 0.0;
  DiffdirTable(1,3) = -1.0;
  DiffdirTable(1,4) = -1.0;
  DiffdirTable(1,5) = -1.0;
  DiffdirTable(1,6) = 0.0;
  DiffdirTable(1,7) = 1.0;
  
  return DiffdirTable;
  
}

// // Gradient calcualtion
// // [[Rcpp::export]]
// mat diffGrad(int ro, int co, mat DiffdirTable, mat Msub, double dx){
//   
//   mat grad(ro, co,fill::zeros); 
//   
//   int a;
//   int x;
//   int y;
//   int ii;
//   int jj;
//   
//   for (ii=1; ii< (ro-1); ii++) {
//      
//     for (jj=1; jj< (co-1); jj++ ){    
//       
//       for (a=0; a < 7; a++) {
// 
//         x = DiffdirTable(0,a);
//         y = DiffdirTable(1,a);
//         
//         //destination(ii+x,jj+y) = qq(ii,jj);
//         grad(ii+x,jj+y) = std::abs(Msub(ii+x,jj+y)-Msub(ii,jj))*dx;
//       }
//     
//     }
//   }
//   
//   return grad;
//   
// }

// Diffusion
// [[Rcpp::export]]
mat Diffusion(int ro, int co, mat DiffdirTable, mat grad, mat Msub, double dt, double D){
  
  
  mat diffloss(ro, co,fill::zeros); 
  mat diffgain(ro, co,fill::zeros);
  
  int a;
  int x;
  int y;
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    
    for (jj=1; jj< (co-1); jj++ ){    
      
      for (a=0; a < 7; a++) {
  
        x = DiffdirTable(0,a);
        y = DiffdirTable(1,a);
        
        diffloss(ii,jj) += -D * grad(ii+x,jj+y)*Msub(ii,jj)*dt;
 
        diffgain(ii+x,jj+y) = -D * grad(ii+x,jj+y)*Msub(ii,jj)*dt; 

      }
    }
  }
  
  return diffloss;
  
}

// diffusion of M
double D; // d can be Dp or Dm 
double dx = extent/rows; // grid cell spacing
double dt = timeincr * pow(dx,2)/D; // time step
double Dp = 0.3;

//diffGradstore = diffGrad(rows, cols, write_DiffdirectionTable(),  M_sub.slice(tt), dx); // matrix
//mat slope exists
diffstore = Diffusion(rows,cols, write_DiffdirectionTable(), slope, M_sub.slice(tt),dt, Dp); //matrix

//maybe instead of concentratin gradient just elevation gradient

M_sub(i,j,tt+1) = M_sub(i,j,tt+1) +  (timeincr * flux_sub(i,j,tt));
M_sub(i,j,tt+1) = M_sub(i,j,tt+1) + diffstore(i,j);
  
  
// distance matrix for plant interference function
// [[Rcpp::export]]
mat Distances(int ro, int co, double kk, double ll){

  mat dist(ro, co,fill::zeros); 
  
  for (ii=1; ii< (ro-1); ii++) {
    
    for (jj=1; jj< (co-1); jj++ ){  
      
      dist(ii,jj) = sqrt(pow((ii-kk),2) + pow((ii-ll),2));
      
    }
  }
  return dist;
}

//Distances(rows,cols,i,j);

/*** R

*/
