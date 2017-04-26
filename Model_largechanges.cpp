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

arma::mat flowdirTable(3,9, fill::zeros);

flowdirTable(0,0) = 0.0;
flowdirTable(0,1) = 0.7853982;
flowdirTable(0,2) = 1.570796;
flowdirTable(0,3) = 2.356195;
flowdirTable(0,4) = 3.141593; 
flowdirTable(0,5) = 3.926991;
flowdirTable(0,6) = 4.712389;
flowdirTable(0,7) = 5.497788;
flowdirTable(0,8) = 6.283186;

flowdirTable(1,0) = 0.0;
flowdirTable(1,1) = -1.0;
flowdirTable(1,2) = -1.0;
flowdirTable(1,3) = -1.0;
flowdirTable(1,4) = 0.0;
flowdirTable(1,5) = 1.0;
flowdirTable(1,6) = 1.0;
flowdirTable(1,7) = 1.0;
flowdirTable(1,8) = 0.0;


flowdirTable(2,0) = 1.0;
flowdirTable(2,1) = 1.0;
flowdirTable(2,2) = 0.0;
flowdirTable(2,3) = -1.0;
flowdirTable(2,4) = -1.0;
flowdirTable(2,5) = -1.0;
flowdirTable(2,6) = 0.0;
flowdirTable(2,7) = 1.0;
flowdirTable(2,8) = 1.0;

// [[Rcpp::export]]
mat Surface(int ro, int co, mat flowdir, mat flowdirTable, mat qq){
  //flowdirTable
  
  const double pi = 3.141593; 
  
  mat destination(ro, co,fill::zeros); 
  
  //double number = (flowdir(i,j)/(pi/4));
  
  int a;
  int x;
  int y;
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    // 
    for (jj=1; jj< (co-1); jj++ ){    
      
      for (a=0; a < 8; a++) {
        
        if (flowdir(ii,jj) == flowdirTable(0,a))
          
        {  
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) += qq(ii,jj);
        }
        
        if ((flowdir(ii,jj) >= flowdirTable(0,a)) && (flowdir(ii,jj) <= flowdirTable(0,a+1))) 
        { 
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) += qq(ii,jj) * (1.0 - ((flowdir(ii,jj) - flowdirTable(0,a+1))/(pi/4)));
          
        }
        
      }
      
     }
  }
  
  return destination;
  
}

  // [[Rcpp::export]]
  List  Veg_cpp(std::string stype) {
    
    // default
    double Zr = 400.0; //mm, Grass
    double k = 12.0;//Saco et al, 2013
    double W0 = 0.2;//Saco et al, 2013
    double gmax = 0.05;//Saco et al, 2013
    double c = 10.0;//Saco et al, 2013
    double k1 = 5.0;//Saco et al, 2013
    double d = 0.24;//Saco et al, 2013 //fraction of plant mortality
    
    
    
    if (stype == "Fantasy Tree") {
      
      Zr = 400.0; 
      k = 12.0;
      W0 = 0.2;
      gmax = 0.05;
      c = 10.0;
      k1 = 5.0;
      d = 0.24;
      
    }
    
    return(Rcpp::List::create(Rcpp::Named("Zr") = Zr,
                              Rcpp::Named("k") = k,
                              Rcpp::Named("W0") = W0,
                              Rcpp::Named("gmax") = gmax,
                              Rcpp::Named("c") = c,
                              Rcpp::Named("k1") = k1,
                              Rcpp::Named("d") = d));
    
  }



// [[Rcpp::export]]
List  Soil_cpp(std::string stype) {
  double psi_sh = -10.0;
  
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51; // cm/day
  double b = 13.48; // neurotheta LMC
  //  double nvg = 1.089;
  // double  avg  = 0.0591;
  double s_fc = 0.364/n; // Field capacity
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = psi_s_bar*(-1e4);
  double spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
  
  
  if (stype == "L Med Clay Stony") {
    // Medium Light Clay
    n = 0.318; // porosity
    // more soil variables for evaporation & losses
    K_s = 3.51; // cm/day
    b = 13.48; // neurotheta LMC
    // nvg = 1.089;
    // avg = 0.0591;
    s_fc = 0.264/n; // Field capacity
    psi_s_bar = -1.5E-3; // This is the bubbling pressure
    hb = psi_s_bar*(-1e4);
    spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
    
  }
  
  if (stype == "S Clay Loam") {
    // Sandy Clay Loam
    n = 0.367; // porosity
    // more soil variables for evaporation & losses
    K_s = 52.08; // cm/day
    b = 6.4069; // neurotheta sandy clay loam
    // avg = 0.0521;
    // nvg = 1.237;
    s_fc = 0.2677/n; // Field capacity
    psi_s_bar = -1.2E-3;
    hb = psi_s_bar*(-1e4);
    spec_y = 0.07;  //difference Fc and por Johnson 1967 says 0.07 
  }
  
  if (stype == "Loamy Sand") {
    // Loamy sand
    n = 0.37; // porosity
    // more soil variables for evaporation & losses
    K_s = 175.3; // cm/day
    b = 4.5206;
    // avg = 0.0641;
    // nvg = 1.344;
    s_fc = 0.2098/n; // Field capacity
    
    psi_s_bar = -0.66E-3; // This is the bubbling pressure
    spec_y = 0.17;  // changed to 0.1 to increase rise in gw, not difference por and fc
  }
  
  if (stype == "H Clay") {
    // Medium Heavy Clay
    n = 0.4473; // porosity
    // more soil variables for evaporation & losses
    K_s = 2.82; // cm/day
    b = 16.1501; // neurotheta Medium heavy clay
    // avg = 0.0613;
    // nvg = 1.086;
    s_fc = 0.3936/n; // Field capacity
    
    psi_s_bar = -1.4e-3;
    
    hb = psi_s_bar*-1.0e5;
    spec_y = 0.05;  // difference por and fc
  }
  
  if (stype == "M Clay") {
    // Medium Clay
    n =0.4391; // porosity
    // more soil variables for evaporation & losses
    K_s = 6.04; // cm/day
    b =  13.5127;
    // avg = 0.0507;
    //doouble nvg = 1.088;
    s_fc = 0.3818/n; // Field capacity
    
    psi_s_bar = -1.75E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e5;
    spec_y = 0.05; // difference por and fc
    
    
  }
  
  if (stype == "C Sand") {
    // Coarse Sand
    n = 0.368; // porosity
    // more soil variables for evaporation & losses
    K_s = 182.68; // cm/day
    b = 4.1152;
    // avg = 0.0712;
    // nvg = 1.392;
    s_fc = 0.1895/n; // Field capacity
    
    psi_s_bar = -0.61E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e5; //mm
    spec_y = 0.27;  // difference por and fc
  }
  
  
  // Other derived parameters
  double s_h = pow(psi_sh/psi_s_bar,-1/b); // soil moisture at hygroscopic point
  double beta = 2*b+4; //page 714 Laio et al., 2001a
  
  // Define parameters for Eagleson function
  // Beta parameters
  double  beta1 = 2+3/b;
  // alpha parameter
  double a1 = 1+(3/2)/(beta1-1);
  
  return(Rcpp::List::create(Rcpp::Named("n") = n,
                            Rcpp::Named("K_s") = K_s,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("hb") = hb,
                            Rcpp::Named("psi_s_bar") = psi_s_bar,
                            Rcpp::Named("s_fc") = s_fc,
                            Rcpp::Named("s_h") = s_h,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("beta1") = beta1,
                            Rcpp::Named("a1") = a1,
                            Rcpp::Named("spec_y") = spec_y));
}


// [[Rcpp::export]]
List Salt_cpp(std::string stype) {
  
  // default, NO SALT from nowhere
  double ConcConst = 0.0; //ConcConst is the concentration of the salt in the infiltrating water in g/l
  double CMgw = 0.0; //CMgw is the goundwater salt concentration  in g/l
  double f = 1; //f is the soil salt leaching efficiency (whether some salt is retained)
  
  
  if (stype == "Groundwater") {
    
    ConcConst = 0.0;
    CMgw = 0.01;
    f = 1;
    
  }
  
  if (stype == "Rain") {
    
    ConcConst = 0.01;
    CMgw = 0.0;
    f = 1;
    
  }
  
  if (stype == "Both") {
    
    ConcConst = 0.01;
    CMgw = 0.01;
    f = 1;
    
  }
  
  if (stype == "None") {
    
    ConcConst = 0.0;
    CMgw = 0.0;
    f = 1;
    
  }
  
  
  return(Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                            Rcpp::Named("CMgw") = CMgw,
                            Rcpp::Named("f") = f));
  
}

// [[Rcpp::export]]
  List SurfaceSoilSaltWBGRID(Rcpp::List soilpar, Rcpp::List vegpar, Rcpp::List saltpar,
                             Rcpp::List dims, NumericVector Rain,
                             double alpha_i, double cn, double Mn, 
                             NumericMatrix slope, NumericMatrix Zras, mat flowdir){
    
    double ConcConst_in = saltpar["ConcConst"];
    double f_in = saltpar["f"];
    double CMgw_in = saltpar["CMgw"];
    
    double n_in = soilpar["n"];
    double b_in = soilpar["b"];
    double K_s_in = soilpar["K_s"];
    double hb_in = soilpar["hb"];
    double psi_s_bar_in = soilpar["psi_s_bar"];
    double h1bar_in = soilpar["h1bar"];
    double s_fc_in = soilpar["s_fc"];
  
    K_s_in = K_s_in*10.0;
  
  double k_in = vegpar["k"];
  double Zr_in = vegpar["Zr"];
  double gmax_in = vegpar["gmax"];
  double c_in = vegpar["c"];
  double k1_in = vegpar["k1"];
  double d_in = vegpar["d"];
  double W0_in = vegpar["W0"];
  
  int i = 0;
  int j= 0;
  int t = 1;
  int tt = 0;
  // int t_old = 0;
  int deltat = 12;
  
  
  float timeincr = 1.0/deltat;
  
  // # qsd has direction of runoff/runon
  // # magnitude of qsd calculates as:
  double c1 = 2.25; // [1/mm] Saco 2013
  double c2 = timeincr * 0.0002 ;  //[m/d] tranformed to [mm/deltat] ; Saco 2013
  
  int rows = dims["rows"];
  int cols = dims["cols"];
  int time = dims["time"];
  
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
  
  arma::cube Gr_sub = arma::zeros(rows, cols, deltat);
  arma::cube Mo_sub = arma::zeros(rows, cols, deltat);
  arma::cube U_salt = arma::zeros(rows, cols, deltat);
  arma::cube L_salt = arma::zeros(rows, cols, deltat);
  
  // advective seed transport with 
  arma::cube qsd_sub = arma::zeros(rows, cols, deltat);
  arma::cube runonsd_sub = arma::zeros(rows, cols, deltat);
  
  // salt mass in h
  arma::cube Smh_sub = arma::zeros(rows, cols, deltat);
  // salt concentration in h
  arma::cube Ch_sub = arma::zeros(rows, cols, deltat);
  
  // Seepage
  arma::cube seep_sub = arma::zeros(rows, cols, deltat);

  arma::cube h = arma::zeros(rows, cols, time);
  arma::cube q = arma::zeros(rows, cols, time);
  arma::cube runon = arma::zeros(rows, cols, time);
  arma::cube In = arma::zeros(rows, cols, time);
  arma::cube P = arma::zeros(rows, cols, time);
  arma::cube Wu = arma::zeros(rows, cols, time);
  arma::cube M = arma::zeros(rows, cols, time);
  arma::cube flux = arma::zeros(rows, cols, time);
  arma::cube CM = arma::zeros(rows, cols, time);
  arma::cube SmI = arma::zeros(rows, cols, time);
  arma::cube SmM = arma::zeros(rows, cols, time);
  arma::cube Svir = arma::zeros(rows, cols, time);
  arma::cube mb = arma::zeros(rows, cols, time);

  
  // advective seed transport with 
  arma::cube qsd = arma::zeros(rows, cols, time);
  arma::cube runonsd = arma::zeros(rows, cols, time);
  
  // Seepage
  arma::cube seep = arma::zeros(rows, cols, time);
  
  // salt mass in h
  arma::cube Smh = arma::zeros(rows, cols, time);
  // salt concentration in h
  arma::cube Ch = arma::zeros(rows, cols, time);

  
  h.fill(20.0);
  P.fill(20.0);
  M.fill(30.0);

  Svir.fill(30.0);
  CM.fill(0.0);
  SmI.fill(0.0);
  SmM.fill(0.0);
  

  
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
      Smh(i,j,0) = 0.0;
      Ch(i,j,0) = 0.0;
      seep(i,j,0) = 0.0;
      
      for (t = 1; t< (time); t++){
        
        for (tt = 0; tt< (deltat-1); tt++){
          
         
          if(tt == 0) {
            h_sub(i,j,tt) = h(i,j,t-1);
            P_sub(i,j,tt) = P(i,j,t-1);
            M_sub(i,j,tt) = M(i,j,t-1);
            CM_sub(i,j,tt) = CM(i,j,t-1);
            SmI_sub(i,j,tt) = SmI(i,j,t-1);
            SmM_sub(i,j,tt) = SmM(i,j,t-1);
            Svir_sub(i,j,tt) = Svir(i,j,t-1);
            Smh_sub(i,j,tt) = Smh(i,j,t-1);
            Ch_sub(i,j,tt) = Ch(i,j,t-1);
            seep_sub(i,j,tt) = seep(i,j,t-1);
          }
          
          double Rain_in;
          
          if ((Rain(t) > 0.0) & (tt == 0)){
            Rain_in = 10.0 * Rain(t);
          } else {
            Rain_in = 0.0;
          }
          
          
          // adjust infiltration rate
          if(h_sub(i,j,tt) < (timeincr * K_s_in)) {
            alpha_i = 1.0;
          } else {
            alpha_i = 1.0-((h_sub(i,j,tt) - (timeincr * K_s_in))/h_sub(i,j,tt));
          }
          
           
          
          q_sub(i,j,tt) = timeincr * OF(h_sub(i,j,tt), cn, Mn, slope(i,j));
          
        
          
          mat runon_store(rows, cols, fill::ones);
          runon_store = Surface(rows, cols, flowdir, flowdirTable, q_sub.slice(tt));
          runon_sub(i,j,tt+1) = runon_store(i,j);  
          
        
          // calculate water depth on soil
          h_sub(i,j,tt+1) =  h_sub(i,j,tt) + Rain_in
            - (timeincr * Infil(h_sub(i,j,tt),P_sub(i,j,tt), alpha_i, k_in, W0_in)) - q_sub(i,j,tt) + seep_sub(i,j,tt) + runon_sub(i,j,tt) ;
           
           
          
          I_sub(i,j,tt) = timeincr * Infil(h_sub(i,j,tt), P_sub(i,j,tt), alpha_i, k_in, W0_in); 
           
          
           
          WU_sub(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_sub(i,j,tt), gmax_in, k1_in); 
           
          
           
          M_sub(i,j,tt+1) = M_sub(i,j,tt) + I_sub(i,j,tt) - WU_sub(i,j,tt+1);
          
          // 
          Gr_sub(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_in);
          //  // // // //Mortality
          Mo_sub(i,j,tt) = timeincr * Mo(P_sub(i,j,tt), M_sub(i,j,tt+1), Svir_sub(i,j,tt),d_in);  
          //  // // //
          //  // // // // Plant biomass balance
          

          if((q_sub(i,j,tt) > c1) & (q_sub(i,j,tt)<c2))
          {
            qsd_sub(i,j,tt) = (1.0/c1) * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
            //qsd_sub(i,j,tt) = c1 * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
          }
          
          if((q_sub(i,j,tt)*c1) > c2){ //(Saco, 2007)
            qsd_sub(i,j,tt) = c2 * P_sub(i,j,tt);
          }


       
          mat runonsd_store(rows, cols, fill::ones);
          runonsd_store = Surface(rows, cols, flowdir, flowdirTable,qsd_sub.slice(tt));
          runonsd_sub(i,j,tt+1) = runonsd_store(i,j); 
          
              
          P_sub(i,j,tt+1) = P_sub(i,j,tt) + Gr_sub(i,j,tt)- Mo_sub(i,j,tt) - qsd_sub(i,j,tt) + runonsd_sub(i,j,tt);
          
          
          
          flux_sub(i,j,tt) = L_n(M_sub(i,j,tt+1),Zras(i,j),n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);  
          
          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) +  (timeincr * flux_sub(i,j,tt));
        
        // Seepage
        
        if(M_sub(i,j,tt+1) > (s_fc_in * n_in * Zr_in)){
          
          seep_sub(i,j,tt+1) = M_sub(i,j,tt+1)-(s_fc_in * n_in * Zr_in);
          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) - seep_sub(i,j,tt+1);
        } 
          
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
          SmI_sub(i,j,tt+1) = SmI_sub(i,j,tt) + (I_sub(i,j,tt) * Ch_sub(i,j,tt));
          

          // #salt mass in soil
          SmM_sub(i,j,tt+1) = SmI_sub(i,j,tt+1) + U_salt(i,j,tt) - L_salt(i,j,tt) - (CM_sub(i,j,tt) * seep_sub(i,j,tt+1));

          //  salt concentration in soil
          CM_sub(i,j,tt+1) = SmM_sub(i,j,tt+1)/M_sub(i,j,tt+1);
          // Rcpp::Rcout <<  CM_sub(i,j,tt);
          
          
          
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
      
          Svir_sub(i,j,tt+1) = n_in * Zr_in * (pow((h1bar_in * 0.1),(1.0/b_in))) * 
          pow((h1bar_in * 0.1 * pow((M_sub(i,j,tt+1)/(n_in * Zr_in)),-b_in))+(3.6 * (CM_sub(i,j,tt+1)*(1.0/58.44))),(-1.0/b_in));
           
 
          // # checking the mass balance
          mb_sub(i,j,tt) = I_sub(i,j,tt) - WU_sub(i,j,tt) + (flux_sub(i,j,tt) * timeincr);
          
         
          
          Smh_sub(i,j,tt+1) = Smh_sub(i,j,tt) + (ConcConst_in * Rain_in) + (CM_sub(i,j,tt+1) * seep_sub(i,j,tt+1)) - (q_sub(i,j,tt)*Ch_sub(i,j,tt)) +  (runon_sub(i,j,tt)*Ch_sub(i,j,tt));
          
          Ch_sub(i,j,tt+1) = (Smh_sub(i,j,tt+1)/h_sub(i,j,tt+1));
          
     
         
          
          
                                            }
        
        h(i,j,t) = h_sub(i,j,(deltat-2));
        Wu(i,j,t) = WU_sub(i,j,(deltat-2));
        P(i,j,t) = P_sub(i,j,(deltat-2));
        M(i,j,t) = M_sub(i,j,(deltat-2));
        CM(i,j,t) = CM_sub(i,j,(deltat-2));
        SmI(i,j,t) = SmI_sub(i,j,(deltat-2));
        SmM(i,j,t) = SmM_sub(i,j,(deltat-2));
        Svir(i,j,t) = Svir_sub(i,j,(deltat-2));
        
        Ch(i,j,t) = Ch_sub(i,j,(deltat-2));
        Smh(i,j,t) = Smh_sub(i,j,(deltat-2));
        seep(i,j,t) = seep_sub(i,j,(deltat-2));

        

        double sumI = 0.0;
        double sumq = 0.0;
       double sumrunon = 0.0;
        double sumflux = 0.0;
        double summb = 0.0;
        double sumqsd = 0.0;
        double sumrunonsd = 0.0;
        
        for(int tt = 0; tt < deltat; tt++)
        {
          sumI += I_sub(i,j,tt);
          sumrunon += runon_sub(i,j,tt);
          sumq += q_sub(i,j,tt);
          sumflux += flux_sub(i,j,tt) * timeincr;
          summb += mb_sub(i,j,tt);
          
          sumqsd += qsd_sub(i,j,tt);
          sumrunonsd += runonsd_sub(i,j,tt);
        }
        
        q(i,j,t) = sumq;
        runon(i,j,t) = sumrunon;
        In(i,j,t) = sumI;
        flux(i,j,t) = sumflux;
        mb(i,j,t) = summb;
        qsd(i,j,t) = sumqsd;
        runonsd(i,j,t) = sumrunonsd;
        
      }
    }
  }
  

  List fields;
  arma::field<arma::cube> f1(18);
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
  f1( 13 ) = Smh;
  f1( 14 ) = Ch;
  f1( 15 ) = qsd;
  f1( 16 ) = runonsd;
  f1( 17 ) = seep;
  fields["field<cube>"] = f1;
  
  
  List output = List::create(_["fields : field<cube>"] = f1);
  
  return output;
  

}

/*** R
result <-SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,
                               saltpar = saltpar1, dims = list(rows=rows,cols=cols,time=time),
                               alpha_i =1.0, cn=0.01, Mn=0.04, Rain=Rain, slope=slp_matrix,Zras=Zras_matrix, flowdir = flowdir)
result$fields[[4]]
  */
