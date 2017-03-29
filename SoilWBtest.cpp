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
  
  double s=M/(n*Zr); // extract n and Zr from list and define them
  double psi = hb*pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));// define hb and b
  double m = 2 + 3/b;
  double qf = (pow((Z/hb),(-m))-(pow((psi/hb),-m)))/((1+pow((psi/hb),-m))+(m-1)*pow((Z/hb),-m));
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
List SurfaceWB(double alpha_i=1.0, double cn = 0.4, double Mn =10.0, double Rain =10.0){ // Surface Balance
  
  
  int i = 0;
  int j= 0;
  int t = 1;
  int tt = 0;
  int t_old = 0;
  int deltat = 4.0;
  
  float timeincr = 1/deltat;
  
  int rows = 2.0;
  int cols = 2.0;
  
  int time = 3.0;
  double slope = 0.001;
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  List vegpar = (Rcpp::List::create(Rcpp::Named("k") = k,
                                    Rcpp::Named("W0") = W0));
  double K_s = 3.51*10.0; // cm/day
  List soilpar= Rcpp::List::create(Rcpp::Named("K_s") = K_s);
  
  
  // //h sub
  double h_sub[rows][cols][deltat];
  
  // //qsub
  double q_sub[rows][cols][deltat];
  // //runonsub
  double runon_sub[rows][cols][deltat];
  //  Isub
  double I_sub[rows][cols][deltat];
  
  double h[rows][cols][time];
  double P[rows][cols][time];
  double In[rows][cols][time];
  
  P[0][0][0]=0.0;
  P[0][0][1]=10.0;
  // M[1][1][1]=10.0;
  h[0][0][0]=10.0;
  In[0][0][0]=0.0;
  
  
  double rn[rows][cols];
  
  //    for (i = 0; i< rows; i++) {
  
  //      for (j = 0; j< cols; j++ ){
  
           for (t = 0; t< time; t++){
  
              for (tt = 0; tt< (deltat); tt++){
    
    if(tt == 0) {
      t_old = t-1;
    } else {
      t_old = tt;
    }
    
    
    // calculation of sub daily runoff and runon
    runon_sub[i][j][tt] = rn[i][j]*q_sub[i][j][t_old];
    q_sub[i][j][tt] = OF(h[i][j][t_old], cn, Mn, slope)*timeincr;
    
    
    double Rain_in;
    
    if ((Rain > 0.0) & (tt == 0)){
      Rain_in = 10.0*Rain;
    } else {
      Rain_in = 0.0;
    }
    
    
    
    // calculate water depth on soil
    h_sub[i][j][tt] =  h[i][j][t_old] + Rain_in 
    - Infil(h[i][j][t_old], P[i][j][t_old], alpha_i, vegpar["k"], vegpar["W0"])
    -q_sub[i][j][tt] + runon_sub[i][j][tt];
    
    //     // adjust infiltration rate
        if(h_sub[i][j][tt] < (K_s*timeincr)) {
          alpha_i = 1.0;
        } else {
          alpha_i = 1-(h_sub[i][j][tt] - (K_s*timeincr))/h_sub[i][j][tt];
        }
    // 
        I_sub[i][j][tt] = Infil(h[i][j][t_old], P[i][j][t_old],
                                alpha_i, vegpar["k"], vegpar["W0"])*timeincr;

    // 
  
  // 
      // # Aggregating the substep results to daily values.


      h[i][j][t] = h_sub[i][j][deltat];

      double sumI = 0.0;

      for(int tt = 0; tt < deltat; ++tt)
      {

        sumI += I_sub[i][j][tt];

      }


      In[i][j][t] = sumI;
  // 
  // //    }
  // //      }
  // //    }
  // 
  //     return(Rcpp::List::create(Rcpp::Named("I_sub") = I_sub[rows][cols][deltat],
  //                               Rcpp::Named("I") = In[rows][cols][time],
  //                               Rcpp::Named("P") = P[rows][cols][time],
  //                               Rcpp::Named("h_sub") = h_sub[rows][cols][deltat]));
  
              }}
  
  return(Rcpp::List::create(Rcpp::Named("h_sub") = h_sub[rows][cols][deltat],
                            Rcpp::Named("q_sub") = q_sub[rows][cols][deltat],
                            Rcpp::Named("I_sub") = I_sub[rows][cols][deltat],
                            Rcpp::Named("runon_sub") = runon_sub[rows][cols][deltat]));
  
}

/*** R
SurfaceWB()
  */