#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*(P+k*W0)/(P+k);
  return I;
}
// [[Rcpp::export]]
double OF(double h, double cn, double Mn, double slope){
  
  double qq = (cn/Mn)*(pow(h,1.666667))*sqrt(slope);
  return qq;
}

// [[Rcpp::export]]
double L_n(double M, double Z, double n, double Zr, double b, double hb, double K_s, double psi_s_bar){
  
  double s=M/(n*Zr); //
  double psi = hb*pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));// define hb and b
  double m = 2 + 3/b;
  double qf = (pow((Z/hb),(-m))-(pow((psi/hb),-m)))/((1+pow((psi/hb),-m))+(m-1)*pow((Z/hb),-m)); // Salvucci 1993

  
  double flux = K_s * qf;
  
  return flux;
  
}


// qf <-((Z/hb)^(-m)-(psi/hb)^(-m))/(1+(psi/hb)^(-m)+(m-1)*(Z/hb)^(-m))
  // q<- (s^(2*soilpar$b+3))-(1+((3/2)/(m)))*(hb/(Z))^m 

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
  
  return Mort;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 Infil(20.0,10.0, 1.0,12,0.2)*(1/12) 



# I_sub[i][j][tt] = Infil(h[i][j][t_old], P[i][j][t_old],
#                         alpha_i, vegpar["k"], vegpar["W0"])*timeincr;

OF(10, 0.01, 0.04, 0.001)

# double OF(double h, double cn, double Mn, double slope){
# 
L_n(200.0, 8000.0, 0.34, 400.0, 13.2, 150, 35.5, -0.0015)

WU(200, 110, 0.05, 5 )
Mo(20, 100, 90, 0.25)

*/
