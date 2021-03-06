#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

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
  double c = 5.0;//Saco et al, 2013
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

// List vegpar = veg_simple();
// double k_in = vegpar["k"];
// double Zr_in = vegpar["Zr"];
// double gmax_in = vegpar["gmax"];
// double c_in = vegpar["c"];
// double k1_in = vegpar["k1"];
// double d_in = vegpar["d"];
// double W0_in = vegpar["W0"];

// Soil function
// [[Rcpp::export]]
List soil_simple() {
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51*10.0; // mm/day
  double b = 13.48; // neurotheta LMC
  
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = -psi_s_bar*(10E5); //mm
  double h1bar = -psi_s_bar;

 
  
  List soilpar= Rcpp::List::create(Rcpp::Named("n") = n,
                                   Rcpp::Named("b") = b,
                                   Rcpp::Named("K_s") = K_s,
                                   Rcpp::Named("hb") = hb,
                                   Rcpp::Named("psi_s_bar") = psi_s_bar,
                                   Rcpp::Named("h1bar") = h1bar);
  
  return(soilpar);
  
}

// List soilpar = soil_simple();
// double n_in = soilpar["n"];
// double b_in = soilpar["b"];
// double K_s_in = soilpar["K_s"];
// double hb_in = soilpar["hb"];
// double psi_s_bar_in = soilpar["psi_s_bar"];
// double h1bar_in = soilpar["h1bar"];


// salt function
// [[Rcpp::export]]
List salt_simple() {
  double ConcConst = 0.1;
  double CMgw = 0.0;
  double f = 1.0;
  
  List saltpar= Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                                   Rcpp::Named("f") = f,
                                   Rcpp::Named("CMgw") = CMgw);
  return(saltpar);
}

// List saltpar = salt_simple();
// double ConcConst_in = saltpar["ConcConst"];
// double f_in = saltpar["f"];
// double CMgw_in = saltpar["CMgw"];

// [[Rcpp::export]]
List SurfaceSoilSaltWBGRID(Rcpp::List soilpar, Rcpp::List vegpar, Rcpp::List saltpar,
                           Rcpp::List dims,
                           double alpha_i, double cn, double Mn, double Rain, 
                           double slope, double Zras){ // Surface Balance , 
  
  
  int i = 0;
  int j= 0;
  int t = 1;
  int tt = 0;
  // int t_old = 0;
  int deltat = 12;
  
  float timeincr = 1.0/deltat;


  
  double n_in = soilpar["n"];
  double b_in = soilpar["b"];
  double K_s_in = soilpar["K_s"];
  double hb_in = soilpar["hb"];
  double psi_s_bar_in = soilpar["psi_s_bar"];
  double h1bar_in = soilpar["h1bar"];
  

  double k_in = vegpar["k"];
  double Zr_in = vegpar["Zr"];
  double gmax_in = vegpar["gmax"];
  double c_in = vegpar["c"];
  double k1_in = vegpar["k1"];
  double d_in = vegpar["d"];
  double W0_in = vegpar["W0"];
  
  double ConcConst_in = saltpar["ConcConst"];
  double f_in = saltpar["f"];
  double CMgw_in = saltpar["CMgw"];
  
  int rows = dims["rows"];
  int cols = dims["cols"];
  int time = dims["time"];
  
  // //h sub
  

  arma::cube h_sub = arma::zeros(rows, cols, deltat);
  arma::cube q_sub= arma::zeros(rows, cols, deltat);
  arma::cube runon_sub= arma::zeros(rows, cols, deltat);
  arma::cube I_sub= arma::zeros(rows, cols, deltat);
  arma::cube P_sub= arma::zeros(rows, cols, deltat);
  arma::cube WU_sub= arma::zeros(rows, cols, deltat);
  arma::cube M_sub= arma::zeros(rows, cols, deltat);
  arma::cube flux_sub= arma::zeros(rows, cols, deltat);
  arma::cube CM_sub= arma::zeros(rows, cols, deltat);
  arma::cube SmI_sub= arma::zeros(rows, cols, deltat);
  arma::cube SmM_sub= arma::zeros(rows, cols, deltat);
  arma::cube Svir_sub= arma::zeros(rows, cols, deltat);
  arma::cube mb_sub= arma::zeros(rows, cols, deltat);
  
  arma::cube Gr_sub(rows, cols, deltat);
  arma::cube Mo_sub(rows, cols, deltat);
  arma::cube U_salt(rows, cols, deltat);
  arma::cube L_salt(rows, cols, deltat);
=======
  int time = 100;
  
  // //h sub
  
  
  arma::cube h_sub = arma::zeros(rows, cols, deltat);
  arma::cube q_sub = arma::zeros(rows, cols, deltat);
  arma::cube runon_sub = arma::zeros(rows, cols, deltat);
  arma::cube I_sub = arma::zeros(rows, cols, deltat);
  arma::cube P_sub = arma::zeros(rows, cols, deltat);
  arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  arma::cube M_sub = arma::zeros(rows, cols, deltat);
  arma::cube flux_sub = arma::zeros(rows, cols, deltat);
  arma::cube CM_sub = arma::zeros(rows, cols, deltat);
  arma::cube SmI_sub = arma::zeros(rows, cols, deltat);
  arma::cube SmM_sub = arma::zeros(rows, cols, deltat);
  arma::cube Svir_sub = arma::zeros(rows, cols, deltat);
  arma::cube mb_sub = arma::zeros(rows, cols, deltat);
  
  arma::cube Gr_sub(rows, cols, deltat);
  arma::cube Mo_sub(rows, cols, deltat);
  arma::cube U_salt(rows, cols, deltat);
  arma::cube L_salt(rows, cols, deltat);


  arma::cube h(rows, cols, time);
  arma::cube q(rows, cols, time);
  arma::cube runon(rows, cols, time);
  arma::cube In(rows, cols, time);
  arma::cube P(rows, cols, time);
  arma::cube Wu(rows, cols, time);
  arma::cube M(rows, cols, time);
  arma::cube flux(rows, cols, time);
  arma::cube CM(rows, cols, time);
  arma::cube SmI(rows, cols, time);
  arma::cube SmM(rows, cols, time);
  arma::cube Svir(rows, cols, time);
  arma::cube mb(rows, cols, time);

  
  // double h[rows][cols][time];
  // double q[rows][cols][time];
  // double runon[rows][cols][time];
  // double In[rows][cols][time];
  // double P[rows][cols][time];
  // double Wu[rows][cols][time];
  // double M[rows][cols][time];
  // double flux[rows][cols][time];
  // double CM[rows][cols][time];
  // double SmI[rows][cols][time];
  // double SmM[rows][cols][time];
  // double Svir[rows][cols][time];
  // double mb[rows][cols][time];


  h.fill(20.0);
  P.fill(20.0);
  M.fill(30.0);
  
  Svir.fill(30.0);
  CM.fill(0.0);
  SmI.fill(0.0);
  SmM.fill(0.0);
  
  
  // h(0,0,0) = 10.0;
  // P(0,0,0) = 20.0;
  // M(0,0,0) = 30.0;
  // 
  // Svir(0,0,0) = 30.0;
  // CM(0,0,0) = 0.0;
  // SmI(0,0,0) = 0.0;
  // SmM(0,0,0) = 0.0;
  

  
  for (i=0; i< rows; i++) {
    
    for (j=0; j< cols; j++ ){
      

  //initialise cubes at 0
      h(i,j,0)=10.0;
      P(i,j,0)=20.0;
      M(i,j,0)=30.0;
      Svir(i,j,0) = 30.0;
      CM(i,j,0) = 0.0;
      SmI(i,j,0) = 0.0;
      SmM(i,j,0) = 0.0;
      
      
      for (t = 1; t< time; t++){

        
        for (tt = 0; tt< (deltat-1); tt++){
          
         
          if(tt == 0) {
            h_sub(i,j,tt) = h(i,j,t-1);
            P_sub(i,j,tt) = P(i,j,t-1);
            M_sub(i,j,tt) = M(i,j,t-1);
            CM_sub(i,j,tt) = CM(i,j,t-1);
            SmI_sub(i,j,tt) = SmI(i,j,t-1);
            SmM_sub(i,j,tt) = SmM(i,j,t-1);
            Svir_sub(i,j,tt) = Svir(i,j,t-1);
          }
          
          double Rain_in;
          
          if ((Rain > 0.0) & (tt == 0)){
            Rain_in = 10.0*Rain;
          } else {
            Rain_in = 0.0;
          }
          
          
          // adjust infiltration rate
          if(h_sub(i,j,tt) < (timeincr * K_s_in)) {
            alpha_i = 1.0;
          } else {

            alpha_i = 1.0-((h_sub(i,j,tt) - (timeincr * K_s_in))/h_sub(i,j,tt));
          }
          
          q_sub(i,j,tt) = timeincr * OF(h_sub(i,j,tt), cn, Mn, slope);
          runon_sub(i,j,tt) = q_sub(i,j,tt) * 0.5;

     

          
          q_sub(i,j,tt) = timeincr * OF(h_sub(i,j,tt), cn, Mn, slope) ;
          runon_sub(i,j,tt) = q_sub(i,j,tt) * 0.5;
          
           
          
          // calculate water depth on soil
          h_sub(i,j,tt+1) =  h_sub(i,j,tt) + Rain_in

            - (timeincr * Infil(h_sub(i,j,tt), P_sub(i,j,tt), alpha_i, k_in, W0_in)) - q_sub(i,j,tt)  + runon_sub(i,j,tt); //

          
          I_sub(i,j,tt) = timeincr * Infil(h_sub(i,j,tt), P_sub(i,j,tt), alpha_i, k_in, W0_in); 
           
          
           
          WU_sub(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_sub(i,j,tt), gmax_in, k1_in); 
           
          
           
          M_sub(i,j,tt+1) = M_sub(i,j,tt) + I_sub(i,j,tt) - WU_sub(i,j,tt+1);
          
          // 
          Gr_sub(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_in);
          //  // // // //Mortality

          Mo_sub(i,j,tt) = timeincr * Mo(P_sub(i,j,tt), M_sub(i,j,tt+1), Svir_sub(i,j,tt),d_in);  

          //  // // //
          
          //  // // // // Plant biomass balance
          P_sub(i,j,tt+1) = P_sub(i,j,tt) + Gr_sub(i,j,tt)- Mo_sub(i,j,tt);
          
          Rcpp::Rcout <<  P_sub(i,j,tt);
          
          flux_sub(i,j,tt) = L_n(M_sub(i,j,tt+1),Zras,n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);  
          

          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) +  (timeincr * flux_sub(i,j,tt));
        
       

          
          // salt leaching
          //
          if(flux_sub(i,j,tt) < 0.0 ){
            L_salt(i,j,tt) = f_in * CM_sub(i,j,tt) * (timeincr * flux_sub(i,j,tt));
          } else {
            L_salt(i,j,tt) = 0.0;
          }
          
          // salt upflow
          if(flux_sub(i,j,tt) > 0.0 ) {
            U_salt(i,j,tt) = CMgw_in * (timeincr * flux_sub(i,j,tt));
          } else {
            U_salt(i,j,tt) = 0.0;
          }
          


          // # salt mass coming in with infiltration
          SmI_sub(i,j,tt+1) = SmI_sub(i,j,tt) + (I_sub(i,j,tt) * ConcConst_in);
          

          // # salt mass coming in with infiltration
          SmI_sub(i,j,tt+1) = SmI_sub(i,j,tt) + (I_sub(i,j,tt) * ConcConst_in);
          


          // #salt mass in soil
          SmM_sub(i,j,tt+1) = SmI_sub(i,j,tt) + U_salt(i,j,tt) - L_salt(i,j,tt);
          
          //  salt concentration in soil

          CM_sub(i,j,tt) = (SmM_sub(i,j,tt+1)/M_sub(i,j,tt+1))*(1.0/58.44);

          //
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
          Svir_sub(i,j,tt) =  n_in* Zr_in *(pow((h1bar_in * 10.0E-1),(1.0/b_in)))*(h1bar_in * 10.0E-1*pow((M_sub(i,j,tt)/(n_in*Zr_in)),-b_in))+pow((3.6*CM_sub(i,j,tt+1)),(-1.0/b_in));
          if ((i == 0) & (j == 0) & (tt==0)){
            Rcpp::Rcout <<  Svir_sub(i,j,tt) << std::endl;
          }
          


          CM_sub(i,j,tt+1) = (SmM_sub(i,j,tt+1)/M_sub(i,j,tt+1))*(1.0/58.44);

          
          // # checking the mass balance
          mb_sub(i,j,tt) = I_sub(i,j,tt) - WU_sub(i,j,tt) + (flux_sub(i,j,tt) * timeincr);
          

        }
        
        h(i,j,t) = h_sub(i,j,(deltat-2));
        Wu(i,j,t) = WU_sub(i,j,(deltat-2));
        P(i,j,t) = P_sub(i,j,(deltat-2));
        M(i,j,t) = M_sub(i,j,(deltat-2));
        CM(i,j,t) = CM_sub(i,j,(deltat-2));
        SmI(i,j,t) = SmI_sub(i,j,(deltat-2));
        SmM(i,j,t) = SmM_sub(i,j,(deltat-2));
        Svir(i,j,t) = Svir_sub(i,j,(deltat-2));
        

        double sumI = 0.0;
        double sumq = 0.0;
        double sumrunon = 0.0;
        double sumflux = 0.0;
        double summb = 0.0;
        
        for(int tt = 0; tt < deltat; tt++)
        {
          sumI += I_sub(i,j,tt);
          sumrunon += runon_sub(i,j,tt);
          sumq += q_sub(i,j,tt);
          sumflux += flux_sub(i,j,tt) * timeincr;
          summb += mb_sub(i,j,tt);
        }
        
        q(i,j,t) = sumq;
        runon(i,j,t) = sumrunon;
        In(i,j,t) = sumI;
        flux(i,j,t) = sumflux;
        mb(i,j,t) = summb;
        
      }
    }
  }
  
  // List out(Rcpp::List::create(Rcpp::Named("h") = h[rows][cols][time],
  //                             Rcpp::Named("q") = q[rows][cols][time],
  //                             Rcpp::Named("In") = In[rows][cols][time],
  //                             Rcpp::Named("runon") = runon[rows][cols][time],
  //                             Rcpp::Named("Wu") = Wu[rows][cols][time],
  //                             Rcpp::Named("P") = P[rows][cols][time],
  //                             Rcpp::Named("flux") = flux[rows][cols][time],
  //                             Rcpp::Named("M") = M[rows][cols][time],
  //                             Rcpp::Named("SmM") = SmM[rows][cols][time],
  //                             Rcpp::Named("CM") = CM[rows][cols][time],
  //                             Rcpp::Named("mb") = mb[rows][cols][time],
  //                             Rcpp::Named("Svir") = Svir[rows][cols][time],
  //                             Rcpp::Named("SmI") = SmI[rows][cols][time]));
  List fields;
  arma::field<arma::cube> f1(13);
  f1( 0 ) = h;
  f1( 1 ) = q;
  f1( 2 ) = In;
  f1( 3 ) = runon;
  f1( 4 ) = Wu;
  f1( 5 ) = P;
  f1( 6 ) = flux;
  f1( 7 ) = M;
  f1( 8 ) = SmM;
  f1( 9 ) = CM;
  f1( 10 ) = mb;
  f1( 11 ) = Svir;
  f1( 12 ) = SmI;
  fields["field<cube>"] = f1;
  
  
  List output = List::create(_["fields : field<cube>"] = f1);
  
  return output;
  

}

/*** R
# soilpar_in <- soil_simple()
#   vegpar_in <- veg_simple()
#   saltpar_in <- salt_simple()
#   

 result<- SurfaceSoilSaltWBGRID(soilpar=soil_simple(), vegpar=veg_simple(),
                                saltpar = salt_simple(), dims = list(rows=10,cols=10,time=100),
                                alpha_i =1.0, cn=0.01, Mn=0.04, Rain=1.0, slope=0.001,Zras=1000.0)
result$fields[[1]][1:10,1:10,2]
str(result$fields)


# SurfaceSoilSaltWBGRID(alpha_i =1.0, cn=0.01, Mn=0.04, Rain=1.0, slope=0.001,Zras=1000.0, soilpar=soilpar_in, vegpar=vegpar_in,saltpar=saltpar_in) 
 
# hdata <- as.data.frame(result$fields[1])
# write.table(hdata, "C:/Users/acoo7451/Desktop/hdata.txt", sep="\t")
 
 # Grid_run()
  */