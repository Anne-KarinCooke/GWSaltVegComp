#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
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

// vegetation function
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
double b_in = b_in;
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
List SurfaceSoilSaltWBGRID(double alpha_i, double cn, double Mn, double Rain, double Zras,
			List soilpar, List vegpar, List saltpar){ // Surface Balance
  
  
  int i = 0;
  int j= 0;
  int t = 0;
  int tt = 0;
  // int t_old = 0;
  int deltat = 12;
  
   float timeincr = 1/deltat;
  
  int rows = 10;
  int cols = 10;
  
  int time = 1000;
  double slope = 0.001;
    
   
  // //h sub
  double h_sub[rows][cols][deltat];
  // //qsub
  double q_sub[rows][cols][deltat];
  // //runonsub
  double runon_sub[rows][cols][deltat];
  //  Isub
  double I_sub[rows][cols][deltat];
  // // fluxsub
  double flux_sub[rows][cols][deltat];
  // //mbsub
  double mb_sub[rows][cols][deltat];
  // //Gr_sub
  double Gr_sub[rows][cols][deltat];
  // //Mo_sub
  double Mo_sub[rows][cols][deltat];
  // //WU_sub
  double WU_sub[rows][cols][deltat];
  // //P sub
  double P_sub[rows][cols][deltat];
  // //M sub
  double M_sub[rows][cols][deltat];
  // //CM sub
  double CM_sub[rows][cols][deltat];
  // // SmIsub
  double SmI_sub[rows][cols][deltat];
  // // SmMsub
  double SmM_sub[rows][cols][deltat];
  // // Salt leaching
  double L_salt[rows][cols][deltat];
  // //salt rise
  double U_salt[rows][cols][deltat];
  // //Svirsub
  double Svir_sub[rows][cols][deltat];
  
  double h[rows][cols][time];
  double P[rows][cols][time];
  double In[rows][cols][time];
  double M[rows][cols][time];
  // //CM
  double CM[rows][cols][time];
  // // SmI
  double SmI[rows][cols][time];
  // // SmM
  double SmM[rows][cols][time];
  // // flux
  double flux[rows][cols][time];
  // //q
  double q[rows][cols][time];
  // //runon
  double runon[rows][cols][time];
  // //Svir
  double Svir[rows][cols][time];
  double mb[rows][cols][time];
  
  
  // P_sub[0][0][0] =10.0;
  // M_sub[0][0][0] =10.0;
  // h_sub[0][0][0] =10.0;
  // Svir_sub[0][0][0] = 9.0;
  
  P[0][0][0]=100.0;
  M[0][0][0]=100.0;
  h[0][0][0]=50.0;

  // In[0][0][0]=0.0;
  Svir[0][0][0]=100.0;
  CM[0][0][0]=0.0;
  SmI[0][0][0]=0.0;
  SmM[0][0][0]=0.0;


  
 // double rn[rows][cols];
  
  // for (i=0; i< rows; i++) {
  // 
  //      for (j=0; j< cols; j++ ){

          
      
      
  
      for (t = 1; t< time; t++){
        
        for (tt = 0; tt< (deltat); tt++){
          
          if(tt == 0) {
            // q_sub[i][j][tt] = q[i][j][t-1];
            // runon_sub[i][j][tt] = runon[i][j][t-1];
            // I_sub[i][j][tt] = In[i][j][t-1];
            flux_sub[i][j][tt] = flux[i][j][t-1];
            h_sub[i][j][tt] = h[i][j][t-1];
            P_sub[i][j][tt] = P[i][j][t-1];
            M_sub[i][j][tt] = M[i][j][t-1];
            SmI_sub[i][j][tt] = SmI[i][j][t-1];
            SmM_sub[i][j][tt] = SmM[i][j][t-1];
            Svir_sub[i][j][tt] = Svir[i][j][t-1];
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
        - (Infil(h_sub[i][j][tt], P_sub[i][j][tt], alpha_i=1.0, k_in, W0_in)*0.83333) - q_sub[i][j][tt] + runon_sub[i][j][tt];

         // calculation of sub daily runoff and runon
         //rn[i][j]
         
         q_sub[i][j][tt] = OF(h_sub[i][j][tt], cn, Mn, slope) * timeincr;
       
       // Rcpp::Rcout << OF(h_sub[i][j][tt], cn=0.01, Mn=0.04, slope=0.001) * 0.8333;
       Rcpp::Rcout <<q_sub[i][j][tt];
       
       runon_sub[i][j][tt] = q_sub[i][j][tt];
        // adjust infiltration rate
        if(h_sub[i][j][tt] < (K_s_in*timeincr)) {
          alpha_i = 1.0;
        } else {
          alpha_i = 1-((h_sub[i][j][tt] - (K_s_in*timeincr))/h_sub[i][j][tt]);
        }

        I_sub[i][j][tt] = Infil(h_sub[i][j][tt], P_sub[i][j][tt],
                                alpha_i=1.0, k_in=12.0, W0_in=0.2)*0.83333;
   

 //         
        // Water uptake
        WU_sub[i][j][tt] = WU(Svir_sub[i][j][tt],P_sub[i][j][tt], gmax_in, k1_in)*timeincr;
        
        // Growth
        Gr_sub[i][j][tt] = Gr(Svir_sub[i][j][tt], P_sub[i][j][tt], c_in, gmax_in, k1_in)*timeincr;
        //Mortality
        Mo_sub[i][j][tt] = Mo(P_sub[i][j][tt], M_sub[i][j][tt], Svir_sub[i][j][tt],d_in)*timeincr;
        
        // Plant biomass balance
        P_sub[i][j][tt] = P_sub[i][j][tt] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; 
        
        
        // Water balance before drainage
        M_sub[i][j][tt] = M_sub[i][j][tt] + I_sub[i][j][tt] - WU_sub[i][j][tt];
        
        
        // Drainage/Capillary rise (vertical water flux)
        
        flux_sub[i][j][tt] = L_n(M_sub[i][j][tt],Zras,n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);  
        
    // Adjustment for M including flux
        //
        M_sub[i][j][tt] = M_sub[i][j][tt] +  flux_sub[i][j][tt]*timeincr;
        
        // salt leaching
        //
        if(flux_sub[i][j][tt] < 0.0 ) {
          L_salt[i][j][tt] = f_in * CM_sub[i][j][tt] * flux_sub[i][j][tt]*timeincr;
        } else {
          L_salt[i][j][tt] = 0.0;
        }
        
        // salt upflow
        if(flux_sub[i][j][tt] > 0.0 ) {
          U_salt[i][j][tt] = CMgw_in * flux_sub[i][j][tt]*timeincr;
        } else {
          U_salt[i][j][tt] = 0.0;
        }
        
        
        // # salt mass coming in with infiltration
        SmI_sub[i][j][tt] = SmI_sub[i][j][tt] + I_sub[i][j][tt] * ConcConst_in;
        
        // #salt mass in soil
        SmM_sub[i][j][tt] = SmI_sub[i][j][tt] + U_salt[i][j][tt] - L_salt[i][j][tt];
        
        //  salt concentration in soil
         CM_sub[i][j][tt] = (SmM_sub[i][j][tt]/M_sub[i][j][tt])*(1.0/58.44);
         
        // 
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
          Svir_sub[i][j][tt] =  n_in* Zr_in *(pow((h1bar_in * 10.0E-1),
			(1/b_in)))*(h1bar_in * 10.0E-1*pow((M_sub[i][j][tt]/(n_in*Zr_in)),-b_in))+pow((3.6*CM_sub[i][j][tt]),(-1.0/b_in));
          
          
          
          
         // # checking the mass balance
          mb_sub[i][j][tt] = I_sub[i][j][tt] - WU_sub[i][j][tt] + flux_sub[i][j][tt]*timeincr;


        // 
        // # Aggregating the substep results to daily values.
    }
    
    P[i][j][t] = P_sub[i][j][deltat];
    M[i][j][t] = M_sub[i][j][deltat];
    h[i][j][t] = h_sub[i][j][deltat];
    CM[i][j][t] = CM_sub[i][j][deltat];
    SmI[i][j][t] = SmI_sub[i][j][deltat];
    SmM[i][j][t] = SmM_sub[i][j][deltat];
    Svir[i][j][t] = Svir_sub[i][j][deltat];
    
    
    
    double sumI = 0.0;
    double sumflux = 0.0;
    double sumrunon = 0.0;
    double summb = 0.0;
    double sumq = 0.0;
    
    for(int tt = 0; tt < deltat; tt++)
    {
      
      
      sumI += I_sub[i][j][tt];
      sumflux += flux_sub[i][j][tt]*timeincr;
      sumq += q_sub[i][j][tt];
      sumrunon += runon_sub[i][j][tt];
      summb += mb_sub[i][j][tt];
      }
    
      In[i][j][t] = sumI;
      flux[i][j][t] = sumflux;
      q[i][j][t] = sumq;
      runon[i][j][t] = sumrunon;
      mb[i][j][t] = summb;
      

          }
         
  //      }
  // }

List out(Rcpp::List::create(Rcpp::Named("P") = P[rows][cols][time],
                            Rcpp::Named("h") = h[rows][cols][time],
                            Rcpp::Named("M") = M[rows][cols][time],
                            Rcpp::Named("Svir") = Svir[rows][cols][time],
                            Rcpp::Named("q") = q[rows][cols][time],
                            Rcpp::Named("mb") = mb[rows][cols][time],
                            Rcpp::Named("In") = In[rows][cols][time],
                            Rcpp::Named("CM") = CM[rows][cols][time],
                            Rcpp::Named("flux") = flux[rows][cols][time],
                            Rcpp::Named("SmI") = SmI[rows][cols][time],
                            Rcpp::Named("runon") = runon[rows][cols][time],
                            Rcpp::Named("SmM") = SmM[rows][cols][time]));
                            
   
  return out;
  
   
}

//[[Rcpp::export]]
List Grid_run(List vegpar, List soilpar, List saltpar){
  int i;
   int j;
   int rows = 10;
     int cols = 10;
// // for testing remove later
 	double alpha_i;
	double cn;
 	double Mn;
 	double Rain;
	double Zras;

    for (i=0; i< rows; i++) {
         for (j=0; j< cols; j++ ){

          SurfaceSoilSaltWBGRID(alpha_i=1.0, cn = 0.01,
 				       Mn =0.004, Rain =10.0, Zras = 1000.0,
 					vegpar=vegpar,soilpar=soilpar,saltpar=saltpar);

         }
    }
    List Out(Rcpp::List::create(Rcpp::Named("results") = SurfaceSoilSaltWBGRID(alpha_i=1.0, cn = 0.01,
 									Mn =0.004, Rain =10.0, Zras = 1000.0,
									vegpar=vegpar,soilpar=soilpar,saltpar=saltpar)));
    return Out;
 }



/*** R

 
 #  cols=5
 #  rows=5
 # 
 #  Store <- list()
 #  sub_store <- list()
 # 
 soilpar_in <- soil_simple()
 vegpar_in <- veg_simple()
 saltpar_in <- salt_simple()
# 
#   for (i in 1:(rows-1)) {
# 
#     for (j in 1:(cols-1)){
# 
#       sub_store[[j]] <-data.frame(SurfaceSoilSaltWBGRID(alpha_i=1.0, cn = 0.001,
#   									Mn =0.04, Rain =10.0, Zras = 1000.0,
#   									vegpar=vegpar_in,
#   									soilpar=soilpar_in,
#   									saltpar=saltpar_in))
# 
#     }
#     Store[[i]] <- sub_store
#   }
#  Store
#  df<- as.data.frame(Store)
#  df
# 
# 

 Grid_run(vegpar=vegpar_in,soilpar=soilpar_in,saltpar=saltpar_in)


*/