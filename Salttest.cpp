#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SaltBalance_cpp (flux, saltpar, timeincr, CM_sub, ){
  
  
  // //CM sub
  double CM_sub[rows][cols][deltat];
  // // SmIsub
  double SmI_sub[rows][cols][deltat];
  // // SmMsub
  double SmM_sub[rows][cols][deltat];
  // // Isub
  double I_sub[rows][cols][deltat];
  // //Svirsub
  double Svir_sub[rows][cols][deltat];
  // // Salt leaching
  double L_salt[rows][cols][deltat];
  // //salt rise
  double U_salt[rows][cols][deltat];
  
  if(flux_sub[i][j][tt] < 0.0 ) {
    L_salt[i][j][tt] = saltpar["f"] * CM_sub[i][j][tt] * flux_sub[i][j][tt]*timeincr;
  } else {
    L_salt[i][j][tt] = 0.0;
  }
  
  // salt upflow
  if(flux_sub[i][j][tt] > 0.0 ) {
    U_salt[i][j][tt] = CM_sub[i][j][tt] *saltpar["f"] * flux_sub[i][j][tt]*timeincr;
  } else {
    L_salt[i][j][tt] = 0.0;
  }
  
  
  // # salt mass coming in with infiltration
  SmI_sub[i][j][tt] = SmI_old[i][j] + I_sub[i][j][tt] * saltpar["ConcConst"];
  
  // #salt mass in soil
  SmM_sub[i][j][tt] = SmI_sub[i][j][tt] + U_salt[i][j][tt] - L_salt[i][j][tt];
  
  // # salt concentration in soil
  CM_sub[i][j][tt] = (SmM_sub[i][j][tt]/M_sub[i][j][tt])*(1/58.44);         
  
  // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
  Svir_sub[i][j][tt] = soilpar["n"]* vegpar["Zr"]*(pow((soilpar["h1bar"]* pow(10,-1)),(1/soilpar["b"])))*(soilpar["h1bar"]*pow(10,-1)*pow((M_sub[i][j][tt]/(soilpar["n"]*vegpar["Zr"])),-soilpar["b"]))+pow((3.6*CM_sub[i][j][tt]),(-1/soilpar["b"]));
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
