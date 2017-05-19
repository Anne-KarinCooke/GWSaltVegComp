#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]

  arma::cube P_subA = arma::zeros(rows, cols, deltat);
  arma::cube P_subB = arma::zeros(rows, cols, deltat);
  arma::cube P_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Gr_subA = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subB = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Mo_subA = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subB = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube WU_subA = arma::zeros(rows, cols, deltat);
  arma::cube WU_subB = arma::zeros(rows, cols, deltat);
  arma::cube WU_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube P_A = arma::zeros(rows, cols, time);
  arma::cube P_B = arma::zeros(rows, cols, time);
  arma::cube P_C = arma::zeros(rows, cols, time);
  
  arma::cube WU_A = arma::zeros(rows, cols, time);
  arma::cube WU_B = arma::zeros(rows, cols, time);
  arma::cube WU_C = arma::zeros(rows, cols, time);
  
  // of all species together

  arma::cube P_sub = arma::zeros(rows, cols, deltat);
  arma::cube P = arma::zeros(rows, cols, time);
  
  arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  arma::cube WU = arma::zeros(rows, cols, time);


// Alternative Growth function including the carrying capacity         
// [[Rcpp::export]]
double Gr(double M, double P, double c, double gmax, double k1, double P0, double sigmaP, double conc){
  
  double Gro = c*WU(M,P,gmax,k1)*(P0*exp(-sigmaP*conc)-P);
  return Gro;
}
// P0 needs to be defined, sigmaP is different for different species
Gr_sub(i,j,tt)= timeincr * Gr(Svir_sub(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_in, P0, sigmaP, CM_sub(i,j,tt));



// dPe^(t∙Sa)
  
  Mo_sub(i,j,tt) = timeincr * Mo(P_sub(i,j,tt), M_sub(i,j,tt+1), Svir_sub(i,j,tt),d_in); 

// [[Rcpp::export]]
double Mo(double P, double M, double Svir, double d, double dur){
  
  double Mort=P*(d*exp(CM_sub(i,j,tt)*dur);
                   
                   return Mort;
                   
}

// Duration of salinity exposure

// with RollRcpp?

double dur(){
  CM(i,j,t) - CM(i,j,)
}

double length = 90;

for(int i = 0; i < length; i++)
{
  
  sum = vector1[length + j - 1] + vector1[length + j - 2];
  a = sum / length;
  
  if CM(i,j,t) - CM(i,j,t-90) = 
  
        
//### Species A Halophyte

double sigmaPA;
          
//### Species B salt-TOLERANT NON-Halophyte

double sigmaPB;
          

//### Species C salt-SENSITIVE NON-Halophyte

double sigmaPC;
          // highert than sigmaPB

        
        // SEED transport with runoff
                  
  double c20;


// this has to be added
c2=c20*exp(-sigmaP*CM_sub(i,j,tt)) 
  
  //already in main model
  if((q_sub(i,j,tt) > c1) & (q_sub(i,j,tt)<c2))  
  {
    qsd_sub(i,j,tt) = (1.0/c1) * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
    //qsd_sub(i,j,tt) = c1 * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
  }
  
  if((q_sub(i,j,tt)*c1) > c2){ //(Saco, 2007)
    qsd_sub(i,j,tt) = c2 * P_sub(i,j,tt);
  }
                  

 // distance matrix for plant interference function
 // [[Rcpp::export]]
 mat Distances(int ro, int co, double kk, double ll){
   
   mat dist(ro, co,fill::zeros); 
   
   for (ii=1; ii< (ro-1); ii++) {
     // 
     for (jj=1; jj< (co-1); jj++ ){  
       
       dist(ii,jj) = sqrt(pow((ii-kk),2) + pow((ii-ll),2));
       
     }
   }
   return dist;
 }
          

// Interference function 
// [[Rcpp::export]]
double interference(int ro, int co, double kk, double ll, mat Psub, double L, double pi = 3.141593){
  
  mat interf(ro, co,fill::zeros);
  mat w(ro, co,fill::zeros);
  

     w = (1/(2*pi*(L*L)))*exp(-(abs(Distances(ro, co, kk, ll))*abs(Distances(ro, co, kk, ll)))/(2*(L*L)));
     interf = w*Psub;
     
     for (ii=1; ii< (ro-1); ii++) {
       for (jj=1; jj< (co-1); jj++){  
         
     double sum += interf(ii,jj);

           }
        }
  return sum;

      }

  
  double intf = interference(rows, cols, i, jj, q_sub.slice(tt), double L, double pi = 3.141593);
  
  P_sub(i,j,tt+1) = P_sub(i,j,tt) + Gr_sub(i,j,tt)- Mo_sub(i,j,tt) - qsd_sub(i,j,tt) + runonsd_sub(i,j,tt) - intf;


/*** R

*/
