#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*((P+k*W0)/(P+k));
  return I;
}


// [[Rcpp::export]]
List SurfaceSoilSaltWBGRID(double alpha_i, double cn, double Mn, double Rain, double Zras
                           ){ // Surface Balance
  
  
  int i = 0;
  int j= 0;
  int t = 1;
  int tt = 0;
  // int t_old = 0;
  int deltat = 12;
  
  float timeincr = 1/deltat;
  
  int rows = 10;
  int cols = 10;
  
  int time = 100;
  double k_in = 12.0;
  double W0_in = 0.2;
  // //h sub
  double h_sub[rows][cols][deltat];
  double h[rows][cols][time];
  h[0][0][0]=50.0;
  
  for (i=0; i< rows; i++) {

       for (j=0; j< cols; j++ ){
  
        for (t = 1; t< time; t++){
          
          for (tt = 0; tt< (deltat); tt++){
            
            if(tt == 0) {
              h_sub[i][j][tt] = h[i][j][t-1];

      }
      
      double Rain_in;
      
      if ((Rain > 0.0) & (tt == 0)){
        Rain_in = 10.0*Rain;
      } else {
        Rain_in = 0.0;
      }
      
      
      // Rain_in 
      // calculate water depth on soil
      h_sub[i][j][tt+1] =  h_sub[i][j][tt] + Rain_in
        - (Infil(h_sub[i][j][tt], 10.0, alpha_i=1.0, k_in, W0_in)*0.833333); // - q_sub[i][j][tt] + runon_sub[i][j][tt];
     
     Rcpp::Rcout << h_sub[i][j][tt];
    }
  
    
    h[i][j][t] = h_sub[i][j][deltat];
    
  }
       }
  }
      
      List out(Rcpp::List::create(Rcpp::Named("h") = h[rows][cols][time]));
      return out;
  }
  
  // 
  // //[[Rcpp::export]]
  // List Grid_run(){
  //   int i;
  //   int j;
  //   int rows = 10;
  //   int cols = 10;
  //   // // for testing remove later
  //   double alpha_i;
  //   double cn;
  //   double Mn;
  //   double Rain;
  //   double Zras;
  //   
  //   for (i=0; i< rows; i++) {
  //     for (j=0; j< cols; j++ ){
  //       
  //       SurfaceSoilSaltWBGRID(alpha_i =1.0, cn=0.01, Mn=0.04, Rain=10.0, Zras=1000.0);
  //       
  //     }
  //   }
  //   List Out(Rcpp::List::create(Rcpp::Named("results") = SurfaceSoilSaltWBGRID(alpha_i=1.0, cn = 0.01,
  //                                           Mn =0.004, Rain =10.0, Zras = 1000.0)));
  //                                           
  //   return Out;
  // }

  
/*** R
SurfaceSoilSaltWBGRID(alpha_i =1.0, cn=0.01, Mn=0.04, Rain=10.0, Zras=1000.0)
# Grid_run()
*/