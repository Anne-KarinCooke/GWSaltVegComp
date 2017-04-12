#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*((P+k*W0)/(P+k));
  return I;
}
// [[Rcpp::export]]
double OF(double h, double cn, double Mn, double slope){
  
  double qq = (cn/Mn)*(pow(h,1.666667))*sqrt(slope);
  return qq;
}
// [[Rcpp::export]]
double WU(double M, double P, double gmax, double k1 ) {  /// not quite happy with the list item calling...list par
  
  double Wu = gmax*(M/(M+k1))*P;   // function is called WU, output is called Wu, cannot be the same
  return Wu;
}

// [[Rcpp::export]]
double Gr(double M, double P, double c, double gmax, double k1){
  
  double Gro = c*WU(M,P,gmax,k1);
  return Gro;
}

// [[Rcpp::export]]
double Mo(double P, double M, double Svir, double d ){
  
  double Mort=P*(d*(M/Svir));
  
  return Mort;
  
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
// [[Rcpp::export]]
List veg_simple() {
  
  double Zr = 400.0; //mm, Grass
  double gmax = 0.05;//Saco et al, 2013
  double c = 10.0;//Saco et al, 2013
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
                                    Rcpp::Named("W0") = W0));
  return(vegpar);
}

List vegpar = veg_simple();
double k_in = vegpar["k"];
double Zr_in = vegpar["Zr"];
double gmax_in = vegpar["gmax"];
double c_in = vegpar["c"];
double k1_in = vegpar["k1"];
double d_in = vegpar["d"];
double W0_in = vegpar["W0"];

// Soil function
// [[Rcpp::export]]
List soil_simple() {
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51*10.0; // mm/day
  double b = 13.48; // neurotheta LMC
  
  double psi_s_bar = -1.5e-3; // This is the bubbling pressure
  double hb = -psi_s_bar*(10e5); //mm
  double h1bar = -psi_s_bar;

  List soilpar= Rcpp::List::create(Rcpp::Named("n") = n,
                                   Rcpp::Named("b") = b,
                                   Rcpp::Named("K_s") = K_s,
                                   Rcpp::Named("hb") = hb,
                                   Rcpp::Named("psi_s_bar") = psi_s_bar,
                                   Rcpp::Named("h1bar") = h1bar);
  
  return(soilpar);
  
}

List soilpar = soil_simple();
double n_in = soilpar["n"];
double b_in = soilpar["b"];
double K_s_in = soilpar["K_s"];
double hb_in = soilpar["hb"];
double psi_s_bar_in = soilpar["psi_s_bar"];
double h1bar_in = soilpar["h1bar"];


// salt function
// [[Rcpp::export]]
List salt_simple() {
  double ConcConst = 0.0;
  double CMgw = 0.0;
  double f = 1.0;
  
  List saltpar= Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                                   Rcpp::Named("f") = f,
                                   Rcpp::Named("CMgw") = CMgw);
  return(saltpar);
}

List saltpar = salt_simple();
double ConcConst_in = saltpar["ConcConst"];
double f_in = saltpar["f"];
double CMgw_in = saltpar["CMgw"];

// [[Rcpp::export]]
List SurfaceSoilSaltWBGRID(double alpha_i, double cn, double Mn, double Rain, double slope, double Zras,List soilpar, List vegpar, List saltpar){ // Surface Balance , 
  
  
  int i = 0;
  int j= 0;
  int t = 1;
  int tt = 0;
  // int t_old = 0;
  int deltat = 12;
  
  // float timeincr = 1/deltat;
  
  int rows = 10;
  int cols = 10;
  
  int time = 100;
  
  // //h sub
  double h_sub[rows][cols][deltat];
  double q_sub[rows][cols][deltat];
  double runon_sub[rows][cols][deltat];
  double I_sub[rows][cols][deltat];
  double P_sub[rows][cols][deltat];
  double WU_sub[rows][cols][deltat];
  double Gr_sub[rows][cols][deltat];
  double Mo_sub[rows][cols][deltat];
  double M_sub[rows][cols][deltat];
  double flux_sub[rows][cols][deltat];
  //double mb_sub[rows][cols][deltat];
  // 
  double h[rows][cols][time];
  double q[rows][cols][time];
  double runon[rows][cols][time];
  double In[rows][cols][time];
  double P[rows][cols][time];
  double Wu[rows][cols][time];
  double M[rows][cols][time];
  double flux[rows][cols][time];
  
  h[0][0][0] = 10.0;
  P[0][0][0] = 20.0;
  M[0][0][0] = 30.0;
  
  for (i=0; i< rows; i++) {
    
    for (j=0; j< cols; j++ ){
      
      for (t = 1; t< time; t++){
        
        for (tt = 0; tt< (deltat); tt++){
          
          if(tt == 0) {
            h_sub[i][j][tt] = h[i][j][t-1];
            P_sub[i][j][tt] = P[i][j][t-1];
            M_sub[i][j][tt] = M[i][j][t-1];
          }
          
          double Rain_in;
          
          if ((Rain > 0.0) & (tt == 0)){
            Rain_in = 10.0*Rain;
          } else {
            Rain_in = 0.0;
          }
          

          // adjust infiltration rate
          if(h_sub[i][j][tt] < (K_s_in * 0.833333)) {
            alpha_i = 1.0;
          } else {
            alpha_i = 1-((h_sub[i][j][tt] - (K_s_in * 0.833333))/h_sub[i][j][tt]);
          }
          
          q_sub[i][j][tt] = OF(h_sub[i][j][tt], cn, Mn, slope) * 0.833333;
          runon_sub[i][j][tt] = q_sub[i][j][tt] * 0.5;
         
         
          
          // calculate water depth on soil
          h_sub[i][j][tt+1] =  h_sub[i][j][tt] + Rain_in
            - (Infil(h_sub[i][j][tt], 20.0, alpha_i, k_in, W0_in) * 0.833333) - q_sub[i][j][tt]  + runon_sub[i][j][tt]; //
          

          I_sub[i][j][tt] = Infil(h_sub[i][j][tt], P_sub[i][j][tt], alpha_i, k_in, W0_in) * 0.833333; 
          
          
         
           WU_sub[i][j][tt] = WU(M_sub[i][j][tt], P_sub[i][j][tt], gmax_in, k1_in) * 0.833333; //Svir_sub[i][j][tt], for M
        
           M_sub[i][j][tt+1] = M_sub[i][j][tt] + I_sub[i][j][tt] - WU_sub[i][j][tt];
           
        // 
           Gr_sub[i][j][tt] = Gr(M_sub[i][j][tt], P_sub[i][j][tt], c_in, gmax_in, k1_in) * 0.833333;//Svir_sub[i][j][tt],
        //  // // // //Mortality
           Mo_sub[i][j][tt] = Mo(P_sub[i][j][tt], M_sub[i][j][tt], M_sub[i][j][tt],d_in) * 0.833333;  // Svir_sub[i][j][tt] for second M
        //  // // //
        //  // // // // Plant biomass balance
           P_sub[i][j][tt+1] = P_sub[i][j][tt] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt];
           
       
           flux_sub[i][j][tt] = L_n(M_sub[i][j][tt],Zras,n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);  
       
           M_sub[i][j][tt+1] = M_sub[i][j][tt] +  (flux_sub[i][j][tt] * 0.833333); // M_sub[i][j][tt]  or M_sub[tt+1]
           Rcpp::Rcout <<  M_sub[i][j][tt] << std::endl;
        }
        
        h[i][j][t] = h_sub[i][j][deltat];
        Wu[i][j][t] = WU_sub[i][j][deltat];
        P[i][j][t] = P_sub[i][j][deltat];
        M[i][j][t] = M_sub[i][j][deltat];
        
        
        double sumI = 0.0;
        double sumq = 0.0;
        double sumrunon = 0.0;
        double sumflux = 0.0;
        
        for(int tt = 0; tt < deltat; tt++)
        {
          sumI += I_sub[i][j][tt];
          sumrunon += runon_sub[i][j][tt];
          sumq += q_sub[i][j][tt];
          sumflux += flux_sub[i][j][tt] * 0.83333;
        }
        
        q[i][j][t] = sumq;
        runon[i][j][t] = sumrunon;
        In[i][j][t] = sumI;
        flux[i][j][t] = sumflux;
        
      }
    }
  }
  Rcpp::Rcout << "time is" <<  time << std::endl;
  
  List out(Rcpp::List::create(Rcpp::Named("h") = h[rows][cols][time],
                              Rcpp::Named("q") = q[rows][cols][time],
                              Rcpp::Named("In") = In[rows][cols][time],
                              Rcpp::Named("runon") = runon[rows][cols][time],
                              Rcpp::Named("Wu") = Wu[rows][cols][time],
                              Rcpp::Named("P") = P[rows][cols][time],
                              Rcpp::Named("flux") = flux[rows][cols][time],
                              Rcpp::Named("M") = );
                           
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
soilpar_in <- soil_simple()
vegpar_in <- veg_simple()
saltpar_in <- salt_simple()

SurfaceSoilSaltWBGRID(alpha_i =1.0, cn=0.01, Mn=0.04, Rain=1.0, slope=0.001,Zras=1000.0, soilpar=soilpar_in, vegpar=vegpar_in,saltpar=saltpar_in)
# Grid_run()
  */