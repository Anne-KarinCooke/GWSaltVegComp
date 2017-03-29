#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SaltBalance_cpp (double M_sub, double flux_sub, double SmI_old ){
  
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  List vegpar = (Rcpp::List::create(Rcpp::Named("k") = k,
                                    Rcpp::List::create(Rcpp::Named("W0") = W0)));
  double K_s = 3.51*10; // mm/day
  double psi_sh = -10.0;
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double b = 13.48; // neurotheta LMC
  //  double nvg = 1.089;
  // double  avg  = 0.0591;
  double s_fc = 0.364/n; // Field capacity
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = psi_s_bar*(-1e4);
  double spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
  double h1bar = -psi_s_bar;
  double Zr = 400.0;
  
  List soilpar= Rcpp::List::create(Rcpp::Named("n") = n,
                                   Rcpp::Named("b") = b
                                   );
  
  double ConcConst = 0.1;
  double CMgw = 0.1;
  double f = 1;
  
  int deltat = 12;
  float timeincr = 1/deltat;
  
  List saltpar= Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                                   Rcpp::Named("f") = f,
                                   Rcpp::Named("CMgw") = CMgw);

  // //CM sub
  double CM_sub;
  // // SmIsub
  double SmI_sub;
  // // SmMsub
  double SmM_sub;

  double I_sub = 10;

  // //Svirsub
  double Svir_sub;
  // // Salt leaching
  double L_salt;
  // //salt rise
  double U_salt;
  
  if(flux_sub < 0.0 ) {
    L_salt = f * CM_sub * flux_sub*timeincr;
  } else {
    L_salt = 0.0;
  }
  
  // salt upflow
  if(flux_sub > 0.0 ) {
    U_salt = CMgw * f * flux_sub*timeincr;
  } else {
    L_salt = 0.0;
  }
  
  
  // # salt mass coming in with infiltration
  SmI_sub = SmI_old + I_sub * ConcConst;
  
  // #salt mass in soil
  SmM_sub = SmI_sub + U_salt - L_salt;
  
  // # salt concentration in soil
  CM_sub = (SmM_sub/M_sub)*(1/58.44);         
  
  // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
  Svir_sub = n* Zr*(pow((h1bar * pow(10,-1)),(1/b)))*(h1bar*pow(10,-1)*pow((M_sub/(n*Zr)),-b))+pow((3.6*CM_sub),(-1/b));
  
  
  return(Rcpp::List::create(Rcpp::Named("SmM_sub") = SmM_sub,
                            Rcpp::Named("SmI_sub") = SmI_sub,
                            Rcpp::Named("Svir_sub") = Svir_sub));
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
SaltBalance_cpp(100,-15,0.001)
*/
