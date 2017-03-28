#include <Rcpp.h>
using namespace Rcpp;
// SOil Water Balance Function

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
double L_n(double M, double Z, double n, double Zr, double b, double K_s, double psi_s_bar){
  
  double hb = psi_s_bar*10e5;
  double s=M/(n*Zr); // extract n and Zr from list and define them
  double psi = hb*pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));// define hb and b
  double m = 2 + 3/b;
  double qf = pow((Z/hb),(-m))-(pow((psi/hb),-m)/(1+pow((psi/hb),-m)))+(m-1)*pow((Z/hb),-m);
  double flux = K_s * qf;
  
  return flux;
  
}

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




// [[Rcpp::export]]
List WB_cpp(NumericVector Rain, int alpha_i=1, double I_sub = 3.0,
                             double slope=0.001, double Zras=3000.0){
  
  int i;
  int j;

  
   int time=5;
   int deltat=6;
   int rows=5;
   int cols=5;
  float timeincr = 1/deltat;  //
  
  // soil moisture [mm]
  double M[rows][cols][time];
  
  //  plant biomass density
  double P[rows][cols][time];
  // //h
  // double h[rows][cols][time];
  

  // //Svir
  double Svir[rows][cols][time];
  // // flux
  // double flux[rows][cols][time];

  // //h sub
  // double h_sub[rows][cols][deltat];
  // //P sub
  double P_sub[rows][cols][deltat];
  // //M sub
  double M_sub[rows][cols][deltat];
 
  // double I_sub[rows][cols][deltat];
  // //Svirsub
  // double Svir_sub[rows][cols][deltat];
  // // fluxsub
  double flux_sub[rows][cols][deltat];
 
  double Gr_sub[rows][cols][deltat];
  // //Mo_sub
  double Mo_sub[rows][cols][deltat];
  // //WU_sub
  double WU_sub[rows][cols][deltat];
  //
  // // Salt leaching
  // double L_salt[rows][cols][deltat];
  // // //salt rise
  // double U_salt[rows][cols][deltat];
  //

  // double Zras[rows][cols];
  // double cn = 0.4;
  // double Mn = 10.0;
  // default
  double Zr = 400.0; //mm, Grass
  double gmax = 0.05;//Saco et al, 2013
  double c = 0.10;//Saco et al, 2013
  double k1 = 5.0;//Saco et al, 2013
  double d = 0.24;//Saco et al, 2013 //fraction of plant mortality
  
  
  
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  List vegpar = (Rcpp::List::create(Rcpp::Named("k") = k,
                                    Rcpp::Named("Zr") = Zr,
                                    Rcpp::Named("gmax") = gmax,
                                    Rcpp::Named("c") = c,
                                    Rcpp::Named("k1") = k1,
                                    Rcpp::Named("d") = d,
                                    Rcpp::List::create(Rcpp::Named("W0") = W0)));
  
  // double psi_sh = -10.0;
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51*10; // mm/day
  double b = 13.48; // neurotheta LMC
  //  double nvg = 1.089;
  // double  avg  = 0.0591;
  double s_fc = 0.364/n; // Field capacity
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = psi_s_bar*(-1e4);
  // double spec_y = 0.054; //Johnson 1967 says 0.05, specific yield.
  // double h1bar = -psi_s_bar;
  
  List soilpar= Rcpp::List::create(Rcpp::Named("n") = n,
                                   Rcpp::Named("b") = b,
                                   Rcpp::Named("K_s") = K_s,
                                   Rcpp::Named("hb") = hb,
                                   Rcpp::Named("psi_s_bar") = psi_s_bar,
                                   Rcpp::Named("s_fc") = s_fc);
                              
  
  
  double ConcConst = 0.1;
  double CMgw = 1.0;
  double f = 1.0;
  
  List saltpar= Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                                   Rcpp::Named("f") = f,
                                   Rcpp::Named("CMgw") = CMgw);
  
  
  P[1][1][1]=10.0;
  M[1][1][1]=10.0;
  Svir[1][1][1]=0.0;
  
  
  int t_old;
  int t;
  int tt;
  
  if(tt == 1) {
    t_old = t-1;
  } else {
     t_old = tt;
  }
  
  // double Rain_in;
  // 
  // if (Rain[t] > 0.0 & tt == 1){
  //   double Rain_in = 10.0*Rain[t];
  // } else {
  //   double Rain_in = 0.0;
  // }
  
  
  for (i=1; i<= rows; i++) {
    
    for (j=1; j<= cols; j++ ){
      
      for (t=1; t<= time; t++){
        
        for (tt=1; tt<= time; tt++){

          
          // Water uptake
          WU_sub[i][j][tt] = WU(M_sub[i][j][t_old],P[i][j][t_old], vegpar["gmax"], vegpar["k1"])*timeincr;
          // Growth
          Gr_sub[i][j][tt] = Gr(Svir[i][j][t_old], P[i][j][t_old], vegpar["c"], vegpar["gmax"], vegpar["k1"])*timeincr;
          //Mortality
          Mo_sub[i][j][tt] = Mo(P[i][j][t_old], M[i][j][t_old], Svir[i][j][t_old],vegpar["d"])*timeincr;
          // Plant biomass balance
          P_sub[i][j][tt] = P[i][j][t_old] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; /// not sure if this all is ok this way or too close to R
          
          // Water balance before drainage
          M_sub[i][j][tt] = M[i][j][t_old] + I_sub - WU_sub[i][j][tt];
          
          // Drainage/Capillary rise (vertical water flux)
          
          flux_sub[i][j][tt] = L_n(M_sub[i][j][tt],Zras,soilpar["n"],vegpar["Zr"],soilpar["b"],soilpar["K_s"],soilpar["psi_s_bar"]);  
          
          
          // Adjustment for M including flux
          //
          M_sub[i][j][tt] = M_sub[i][j][tt] +  flux_sub[i][j][tt]*timeincr;
  
        }
        
        P[i][j][t] = P_sub[i][j][deltat];
        M[i][j][t] = M_sub[i][j][deltat];
      }
    }
  }
  List out = Rcpp::List::create(Rcpp::Named("P") = P[rows][cols][time],
                                Rcpp::Named("M") = M[rows][cols][time]);
  
  return out;
  
}
/*** R
 Rain <- rep(1, 10)
 WB_cpp(Rain)
*/