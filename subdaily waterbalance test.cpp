#include <Rcpp.h>
using namespace Rcpp;
// SOil Water Balance Function
// [[Rcpp::export]]
WB_cpp(tt, t_old, vegpar, soilpar){
  
  int deltat = 12;
  
  NumericVector flux_sub <deltat>;
  // //mbsub
  double mb_sub[deltat];
  // //Gr_sub
  double Gr_sub[deltat];
  // //Mo_sub
  double Mo_sub[deltat];
  // //WU_sub
  double WU_sub[deltat];
  
  // Water uptake           
  WU_sub[tt] = WU(M_sub[t_old],P[t_old], vegpar["gmax"], vegpar["k1"])*timeincr;
  // Growth            
  Gr_sub[tt] = Gr(Svir[t_old], P[t_old], vegpar["c"], vegpar["gmax"], vegpar["k1"])*timeincr; 
  //Mortality
  Mo_sub[tt] = Mo(P[t_old], M[t_old], Svir[t_old],vegpar["d"])*timeincr;
  // Plant biomass balance             
  P_sub[tt] = P[t_old] + Gr_sub[tt]- Mo_sub[tt]; 
  
  // Water balance before drainage
  M_sub[tt] = M[t_old] + I_sub[tt] - WU_sub[tt];
  
  // Drainage/Capillary rise (vertical water flux)          
  
  flux_sub[tt] = L_n(M_sub[tt],Zras,soilpar["b"],soilpar["K_s"],soilpar["psi_s_bar"]);  
  
  // Adjustment for M including flux
  //
  M_sub[tt] = M_sub[tt] +  flux_sub[tt]*timeincr; 
  
  
  
}
/*** R

*/
