#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List SaltBalance_cpp2D(double M_sub, double flux_sub, double SmI_old ){
  
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
  double CMgw = 1.0;
  double f = 1.0;
  
  int deltat = 12;
  float timeincr = 1/deltat;
  
  List saltpar= Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                                   Rcpp::Named("f") = f,
                                   Rcpp::Named("CMgw") = CMgw);
  
  // //CM sub
  double CM_sub[deltat];
  // // SmIsub
  double SmI_sub[deltat];
  // // SmMsub
  double SmM_sub[deltat];
  
  double I_sub = 2.0;
  
  // //Svirsub
  double Svir_sub[deltat];
  // // Salt leaching
  double L_salt[deltat];
  // //salt rise
  double U_salt[deltat];
  
  int tt;
  for (tt==1; tt< (deltat); tt++){
    
  if(flux_sub < 0.0 ) {
    L_salt[tt] = f * CM_sub[tt] * flux_sub*timeincr;
  } else {
    L_salt[tt] = 0.0;
  }
  
  // salt upflow
  if(flux_sub > 0.0 ) {
    U_salt[tt] = CMgw * f * flux_sub*timeincr;
  } else {
    L_salt[tt] = 0.0;
  }
  
  
  // # salt mass coming in with infiltration
  SmI_sub[tt] = SmI_old + I_sub * ConcConst;
  
  // #salt mass in soil
  SmM_sub[tt] = SmI_sub[tt] + U_salt[tt] - L_salt[tt];
  
  // # salt concentration in soil
  CM_sub[tt] = (SmM_sub[tt]/M_sub)*(1/58.44);         
  
  // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
  Svir_sub[tt] = n* Zr*(pow((h1bar * pow(10,-1)),(1/b)))*(h1bar*pow(10,-1)*pow((M_sub/(n*Zr)),-b))+pow((3.6*CM_sub[tt]),(-1/b));
  
  }
  return(Rcpp::List::create(Rcpp::Named("SmM_sub") = SmM_sub[deltat],
                            Rcpp::Named("SmI_sub") = SmI_sub[deltat],
                            Rcpp::Named("Svir_sub") = Svir_sub[deltat]));
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
SaltBalance_cpp2D(100,100,10)
*/