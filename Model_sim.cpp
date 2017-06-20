#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
mat sub1( mat x, uword e) {
  x.shed_col(e-1);
  x.shed_row(e-1);
  return x;
}

// Infiltration
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*((P+k*W0)/(P+k));
  return I;
}
// Overland flow/Runoff
// [[Rcpp::export]]
double OF(double h, double cn, double Mn, double slope){
  
  double qq = (cn/Mn) * (pow(h,1.666667)) * sqrt(slope);
  return qq;
}
// plant water uptake
// [[Rcpp::export]]
double WU(double M, double P, double gmax, double k1 ) {
  
  double Wu = gmax*(M/(M+k1))*P;
  return Wu;
}
////Alternative Growth function including the carrying capacity
// [[Rcpp::export]]
double Gr(double M, double P, double P_all, double c, double gmax, double k1, double P0) { //, double sigmaP, double CM)
  
  // double Gro = c*WU(M,P,gmax,k1)*(P0*exp(-sigmaP*CM)-P);
  double Gro = c * WU(M,P,gmax,k1) * (1 - (P_all/P0)); // P_all is the biomass of all species together in a cell
  return Gro;
}
// // // simple growth
// // [[Rcpp::export]]
// double Gr(double M, double P, double c, double gmax, double k1){
// 
//   double Gro = c*WU(M,P,gmax,k1);
//   return Gro;
// }

//plant mortality function
// [[Rcpp::export]]
double Mo(double P, double d, double sigmaP, double CM){
  
  double Mort=P * (d * exp(CM * sigmaP));
  
  return Mort;
  
}
// Capillary rise and drainage
// [[Rcpp::export]]
double L_n(double M, double Z, double n, double Zr, double b, double hb, double K_s, double psi_s_bar){
  
  double s=M/(n*Zr); 
  double psi = hb * pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));
  double m = 2 + 3/b;
  double qf = (pow((Z/hb),(-m))-(pow((psi/hb),-m)))/((1+pow((psi/hb),-m))+(m-1)*pow((Z/hb),-m));
  double flux = K_s * qf;
  
  return flux;
  
}
// Runon
// [[Rcpp::export]]
mat Surface(int ro, int co, mat flowdir, mat flowdirTable, mat qq, mat filler){
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
          
          
          destination(ii+x,jj+y) = qq(ii,jj);
        }
        
        if ((flowdir(ii,jj) >= flowdirTable(0,a)) && (flowdir(ii,jj) <= flowdirTable(0,a+1))) 
        { 
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) = qq(ii,jj) * filler(ii,jj)* (1.0 - ((flowdir(ii,jj) - flowdirTable(0,a+1))/(pi/4)));
          
        }
        
      }
      
    }
  }
  
  return destination;
  
}
//lateral subsurface water flow
// [[Rcpp::export]]
mat Subsurface(int ro, int co, mat flowdir, mat flowdirTable, mat M, mat filler, double Dm, double timeincr){
  //flowdirTable
  
  const double pi = 3.141593; 
  
  mat destination(ro, co,fill::zeros); 
  
  
  
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
          
          
          destination(ii+x,jj+y) = timeincr * Dm * M(ii,jj);
        }
        
        if ((flowdir(ii,jj) >= flowdirTable(0,a)) && (flowdir(ii,jj) <= flowdirTable(0,a+1))) 
        { 
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) = timeincr * Dm * M(ii,jj) * filler(ii,jj)* (1.0 - ((flowdir(ii,jj) - flowdirTable(0,a+1))/(pi/4)));
          
        }
        
      }
      
    }
  }
  
  return destination;
  
}
// simply the distances from all other cell to one cell
// [[Rcpp::export]]
mat Distances(int ro, int co, double kk, double ll, double dx){  ///kk and ll are the i and j of the certain cell around with the distance matrix is calculated
  
  mat dist(ro, co,fill::zeros); 
  
  int ii;
  int jj;
  
  
  for (ii=0; ii< (ro); ii++) {
    for (jj=0; jj< (co); jj++ ){  
      
      dist(ii,jj) = sqrt(pow((ii-kk),2) + pow((jj-ll),2))*dx; //*dx for output in m
      
    }
  }
  return dist;
} 
// plant Interference function (competition/dacilitation)
// [[Rcpp::export]]
double interference(int ro, int co, double kk, double ll, mat Psub, double dx, double b1, double b2, double q1, double q2, double pi = 3.141593){  // doubel L
  
  mat interf(ro, co,fill::zeros);
  mat w(ro, co,fill::zeros);
  
  int ii;
  int jj;
  
  
  for (ii=0; ii< (ro); ii++) {
    for (jj=0; jj< (co); jj++){
      
      // w = (1/(2*pi*(L*L))) * exp(-(Distances(ro, co, kk, ll, dx) * Distances(ro, co, kk, ll, dx)))/(2*(L*L));
      
      mat dist = Distances(ro, co, kk, ll, dx);
      
      
      w(ii,jj) = b1 *exp(-pow((dist(ii,jj)/q1),2)) - b2 *exp(-pow((dist(ii,jj)/q2),2));
      
    }
  }
  double sum =0.0;
  for (ii=0; ii< (ro); ii++) {
    for (jj=0; jj< (co); jj++){
      
      interf(ii,jj) = w(ii,jj) * Psub(ii,jj);
      
      sum = arma::accu(interf);
    }
  }
  return sum;
}

// function for seed diffusion directions (all directions)
// [[Rcpp::export]] 
mat write_DiffdirectionTable() {  // this works
  
  arma::mat DiffdirTable(2,8, fill::zeros);
  
  DiffdirTable(0,0) = 0.0;
  DiffdirTable(0,1) = -1.0;
  DiffdirTable(0,2) = -1.0;
  DiffdirTable(0,3) = -1.0;
  DiffdirTable(0,4) = 0.0;
  DiffdirTable(0,5) = 1.0;
  DiffdirTable(0,6) = 1.0;
  DiffdirTable(0,7) = 1.0;
  
  
  DiffdirTable(1,0) = 1.0;
  DiffdirTable(1,1) = 1.0;
  DiffdirTable(1,2) = 0.0;
  DiffdirTable(1,3) = -1.0;
  DiffdirTable(1,4) = -1.0;
  DiffdirTable(1,5) = -1.0;
  DiffdirTable(1,6) = 0.0;
  DiffdirTable(1,7) = 1.0;
  
  return DiffdirTable;
  
}

// // seed Diffusion
// [[Rcpp::export]] 
mat seedDiffusionLoss(int ro, int co, mat Medium , double Dp, double timeincr){
  
  mat diffloss(ro, co,fill::zeros);
  
  
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    
    for (jj=1; jj< (co-1); jj++ ){
      
      
      diffloss(ii,jj) = Dp * Medium(ii,jj) * timeincr; 
      
      
    }
  }
  
  
  
  return diffloss;
  
}
// [[Rcpp::export]] 
mat seedDiffusionGain(int ro, int co, mat DiffdirTable, mat Medium , double Dp, double timeincr){
  
  mat diffgain(ro, co,fill::zeros);
  float divide = 1.0/8.0;
  
  int a;
  int x;
  int y;
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    for (jj=1; jj< (co-1); jj++ ){
      for (a=0; a < 7; a++) {
        
        x = DiffdirTable(0,a);
        y = DiffdirTable(1,a);
        
        diffgain(ii+x,jj+y) = Dp  * Medium(ii,jj) *timeincr * divide;
      }
    }
  }
  return diffgain;
}

// vegetation parameter list function 
// [[Rcpp::export]]
List  Veg_cpp() {
  
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  double gmax = 0.05;//Saco et al, 2013
  double c = 10.0;//Saco et al, 2013
  double k1 = 5.0;//Saco et al, 2013
  double d = 0.24;//Saco et al, 2013 //fraction of plant mortality
  
  
  
  return(Rcpp::List::create(Rcpp::Named("k") = k,
                            Rcpp::Named("W0") = W0,
                            Rcpp::Named("gmax") = gmax,
                            Rcpp::Named("c") = c,
                            Rcpp::Named("k1") = k1,
                            Rcpp::Named("d") = d));
  
}  
// salinity parameter list function 
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
// soil parameter list function 
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
    //double nvg = 1.088;
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
// takes elevation, calls and applies TauDEM, outputs flowdir and slope for every cell
// [[Rcpp::export]]
int call_Taudem(arma::mat B){
  
  B.save("B.mat",arma::raw_ascii); 
  system("R CMD BATCH GeoTiff.R");
  return 0;
}
// [[Rcpp::export]]
mat flowdir_load(){
  arma::mat flowdir_new;
  flowdir_new.load("new_flowdir.mat");
  mat flo = sub1(flowdir_new,1);
  
  return flo;
}
// [[Rcpp::export]]
mat slope_load(){
  
  arma::mat slp_matrix_new;
  slp_matrix_new.load("slp_matrix.mat");
  mat slp = sub1(slp_matrix_new,1);
  
  return slp;
}

// translates flowdir from TauDEM (0 to 2pi) into the indeces of the corresponding cells
// [[Rcpp::export]] 
mat write_flowdirTable() {
  
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
  
  return flowdirTable;
  
}
// Creating a elevation matrix (DEM) not used anymore
/////////[[Rcpp::export]]
// mat write_elev(int rows, int cols, double gslp, double ext){  // has to be rewritten to account for statistical propertis of suface elevation heterogeneity
//   
//   int iii;
//   mat elev(rows, cols, fill::randn);
//   
//   elev.row(1) = elev.row(1) + (gslp * ext);
//   
//   for (iii=0; iii< rows; iii++) {
//     
//     elev.row(iii) = elev.row(iii)+(gslp * (ext-((ext/rows) * (iii-1))));
//   }
//   return elev;
// }
// Just a function to find a condition and produce a submatrix (used in following function)
// [[Rcpp::export]]
arma::vec matrix_sub(arma::mat y, double B) {
  // arma::mat Z = M * M.t();                        
  arma::vec v = y.elem( find( y > B )); // Find IDs & obtain values
  return v;
}
// calculates difference in height from one cell to the next highest cell of its surrounding cells
// [[Rcpp::export]]
mat diff_next_highest(int rows, int cols, mat elev){
  
  arma::mat diff_next_highest_cell(rows, cols, fill::zeros);
  
  int i;
  int j;
  for (i=1; i< (rows-1); i++) {
    
    for (j=1; j< (cols-1); j++ ){
      
      
      mat y = elev.submat(i-1, j-1, i+1, j+1); //D8
      
      
      vec rest = matrix_sub(y,elev(i,j));
      
      NumericVector rest1 = as<NumericVector>(wrap(rest));
      
      double minimum = min(rest1);
      
      diff_next_highest_cell(i,j) =  minimum- elev(i,j);
      
      
    }
  }
  
  return diff_next_highest_cell;
}
// adding slope to elevation grid
// [[Rcpp::export]]
mat add_Slope_to_Elev(mat elev, int rows, double gslp, double ext){  
  
  int iii;
  elev.row(0) = elev.row(0) + (gslp * ext);
  
  for (iii=1; iii< rows; iii++) {
    
    elev.row(iii) = elev.row(iii) + (gslp * (ext-((ext/rows) * (iii-1))));
  }
  return elev;
}
// [[Rcpp::export]]
mat Z_matrix(mat elev, int ro, int co, double Z){
  
  mat m(ro, co, fill::zeros);
  int ii;
  int jj;
  
  for (ii = 0; ii< ro; ii++){
    
    for (jj = 0; jj< co; jj++){
      
      m(ii,jj) = elev(ii,jj) * 1000.0 + Z;
      
    }
  }
  return m;
  
}
//gaussian function for halophytes
// [[Rcpp::export]]
double gauss(double CM, double mean, double var, double pi = 3.141593){
  
  double gauss = (1/(var*sqrt(2*pi)))*exp(pow((CM-mean),2)/(2*(var*var)));
  return gauss;
}

// the core of the model
// [[Rcpp::export]]
List SurfaceSoilSaltWBGRID(Rcpp::List soilpar, Rcpp::List vegpar,
                           Rcpp::List dims, NumericVector Rain,
                           Rcpp::List fixedInput, Rcpp::List simInput, mat elev, mat Zras
) {
  
  
  // spatial and temporal dimensions of the domain
  int rows = dims["rows"];
  int cols = dims["cols"];
  // duration of simulation
  int time = dims["time"];
  // groundwater depth
  double Z = dims["Z"];
  // extension of the domain
  double ext = dims["ext"];
  // temporal discretization, subdaily time steps
  int deltat = fixedInput["deltat"];
  
  // hillslope
  double gslp = simInput["gslp"]; 
  
  // salinity
  
  // salt concentration rainfall/irrigation
  double ConcConst_in = simInput["ConcConst"];
  // salt concentration groundwater 
  double CMgw_in = simInput["CMgw"];
  // salt leaching efficiency
  double f_in = fixedInput["f"];
  
  
  // soil parameters
  // porosity
  double n_in = soilpar["n"];
  // campbells b
  double b_in = soilpar["b"];
  /// saturated hudraulic conductivity
  double K_s_in = soilpar["K_s"];
  K_s_in = K_s_in*10.0;  // corrected to be in mm/d
  // bubbling pressure
  double hb_in = soilpar["hb"];
  // matrix potential
  double psi_s_bar_in = soilpar["psi_s_bar"];
  double h1bar_in = soilpar["h1bar"];
  // field capacity
  double s_fc_in = soilpar["s_fc"];
  // soil moisture diffusivity
  double Dm = fixedInput["Dm"]; 
  
  // infiltration rate
  double alpha_i = fixedInput["alpha_i"]; 
  
  // Kinematic wave (runoff) paramters
  double cn = fixedInput["cn"];  // conversion factor
  double Mn = fixedInput["Mn"]; // Mannings n
  
  // vegetation parameters
  //root depth
  double Zr_inA = simInput["ZrA"];
  double Zr_inB = simInput["ZrB"];
  double Zr_inC = simInput["ZrC"];
  
  double k_in = vegpar["k"];
  double gmax_in = vegpar["gmax"];
  // growth conversion factor
  double c_in = vegpar["c"];
  
  double W0_in = vegpar["W0"];
  
  // half saturation constant
  double k1_inA = simInput["k1A"];
  double k1_inB = simInput["k1B"];
  double k1_inC = simInput["k1C"];
  
  // mortality factor, half life?
  double d_inA = simInput["dA"];
  double d_inB = simInput["dB"];
  double d_inC = simInput["dC"];
  
  
  // SPECIES A interference parameters, competition and facilitation
  // double b1A = simInput["b1A"];
  // double b2A = simInput["b2A"];
  double q1A = simInput["q1A"];
  double q2A = q1A/0.6;
  
  // SPECIES B interference parameters, competition and facilitation
  // double b1B = simInput["b1B"];
  // double b2B = simInput["b2B"];
  double q1B = simInput["q1B"];
  double q2B = q1B/0.6;
  
  // SPECIES C interference parameters, competition and facilitation
  // double b1C = simInput["b1C"];
  // double b2C = simInput["b2C"];
  double q1C = simInput["q1C"];
  double q2C = q1C/0.6;
  
  
  double zeta = fixedInput["zeta"]; // relative importance of non-local effects vs local effects (impact of interference)
  
  // ecosystem carrying capacity of plant biomass density
  double P0 =  fixedInput["P0"];
  // plant sensitivity to salinity species A
  double sigmaPA = simInput["sigmaPA"];
  // plant sensitivity to salinity species B
  double sigmaPB = simInput["sigmaPB"];
  // plant sensitivity to salinity  species C
  double sigmaPC = simInput["sigmaPC"];
  
  
  // seed dispersal and diffusion parameters
  double c1 = fixedInput["c1"]; // [1/mm] 
  double c02 = fixedInput["c02"];  //[m/d] 
  //seed diffusivity
  double Dp = fixedInput["Dp"]; 
  
  
  
  int i = 0;
  int j = 0;
  int t = 1;
  int tt = 0;
  
  
  // temporal discretization
  float timeincr = 1.0/deltat;
  // cell size
  double dx = ext/rows; 
  //double dts = timeincr * ((dx*dx)/Dp);
  
  
  // old:
  // // Creating the DEM
  // mat elev = write_elev(rows, cols, gslp, ext);
  // // elev.save("B.txt", arma::raw_ascii);
  // // system("R CMD BATCH GeoTiff.R");
  // // load from disk
  // arma::mat flowdir_new;
  // flowdir_new.load("flowdir_new.txt");
  // //mat flowdir_new = sub1(flowdir,1);
  // arma::mat slope;
  // slope.load("slp_matrix.txt");
  // //mat slope = sub1(slop,1);
  
  
  
  //// Creating the DEM
  // mat elev = write_elev(rows, cols, gslp, ext);
  // matrix with elevation differences
  mat diff = diff_next_highest(rows, cols, elev);
  // flow directions
  int bla = call_Taudem(elev);
  
  mat flowdir_new = flowdir_load();
  //slope of each cell
  mat slope = slope_load();
  
  
  /// storage cubes (arrays) for the subdaily time steps tt:
  
  /// HYDROLOGY
  // overland flow depth
  arma::cube h_sub = arma::zeros(rows, cols, deltat);
  // runoff
  arma::cube q_sub = arma::zeros(rows, cols, deltat);
  //runon
  arma::cube runon_sub = arma::zeros(rows, cols, deltat);
  // subsurface water runon
  arma::cube Subsrunon_sub = arma::zeros(rows, cols, deltat);
  // infiltration
  arma::cube I_sub = arma::zeros(rows, cols, deltat);
  //water mass balance
  arma::cube mb_sub = arma::zeros(rows, cols, deltat);
  //salt mass balance
  arma::cube mb_salt_sub = arma::zeros(rows, cols, deltat);
  //soil moisture
  arma::cube M_sub = arma::zeros(rows, cols, deltat);
  // vertical water flux, capillary rise and drainage
  arma::cube flux_sub = arma::zeros(rows, cols, deltat);
  // Seepage
  arma::cube seep_sub = arma::zeros(rows, cols, deltat);
  
  
  /// VEGETATION
  /// plant biomass density
  // arma::cube P_sub = arma::zeros(rows, cols, deltat);
  // // plant water uptake
  // arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  // plant growth
  arma::cube Gr_sub = arma::zeros(rows, cols, deltat);
  // plant mortality
  arma::cube Mo_sub = arma::zeros(rows, cols, deltat);
  
  // SPECIES A
  // advective seed transport with runoff 
  arma::cube qsd_subA = arma::zeros(rows, cols, deltat);
  // advective seed transport with runon
  arma::cube runonsd_subA = arma::zeros(rows, cols, deltat);
  //SPECIES B
  // advective seed transport with runoff 
  arma::cube qsd_subB = arma::zeros(rows, cols, deltat);
  // advective seed transport with runon
  arma::cube runonsd_subB = arma::zeros(rows, cols, deltat);
  // SPECIES C
  // advective seed transport with runoff 
  arma::cube qsd_subC = arma::zeros(rows, cols, deltat);
  // advective seed transport with runon
  arma::cube runonsd_subC = arma::zeros(rows, cols, deltat);
  
  
  /// SALT
  // soil solution salt concentration
  arma::cube CM_sub = arma::zeros(rows, cols, deltat);
  // salt mass in infiltrating water
  arma::cube SmI_sub = arma::zeros(rows, cols, deltat);
  //soil solution salt mass
  arma::cube SmM_sub = arma::zeros(rows, cols, deltat);
  // virtual saturation
  arma::cube Svir_sub = arma::zeros(rows, cols, deltat);
  // salt transport with subsurface lateral flow
  arma::cube runonSubsSalt_sub = arma::zeros(rows, cols, deltat);
  // salt transport with surface runon
  arma::cube salt_runon_sub = arma::zeros(rows, cols, deltat);
  // salt rising up with capillary rise
  arma::cube U_salt = arma::zeros(rows, cols, deltat);
  // salt leached due to drainage
  arma::cube L_salt = arma::zeros(rows, cols, deltat);
  // salt mass in h
  arma::cube Smh_sub = arma::zeros(rows, cols, deltat);
  // salt concentration in h
  arma::cube Ch_sub = arma::zeros(rows, cols, deltat);
  
  
  
  /// storage cubes (arrays for the daily timesteps t:
  
  arma::cube h = arma::zeros(rows, cols, time);
  arma::cube q = arma::zeros(rows, cols, time);
  arma::cube runon = arma::zeros(rows, cols, time);
  arma::cube In = arma::zeros(rows, cols, time);
  // arma::cube P = arma::zeros(rows, cols, time);
  // arma::cube Wu = arma::zeros(rows, cols, time);
  arma::cube M = arma::zeros(rows, cols, time);
  arma::cube flux = arma::zeros(rows, cols, time);
  arma::cube CM = arma::zeros(rows, cols, time);
  arma::cube SmI = arma::zeros(rows, cols, time);
  arma::cube SmM = arma::zeros(rows, cols, time);
  arma::cube Svir = arma::zeros(rows, cols, time);
  arma::cube mb = arma::zeros(rows, cols, time);
  arma::cube mb_salt = arma::zeros(rows, cols, time);
  
  arma::cube seep = arma::zeros(rows, cols, time);
  arma::cube Smh = arma::zeros(rows, cols, time);
  arma::cube Ch = arma::zeros(rows, cols, time);
  arma::cube runonSubsSalt = arma::zeros(rows, cols, time);
  arma::cube salt_runon = arma::zeros(rows, cols, time);
  arma::cube Subsrunon = arma::zeros(rows, cols, time);
  
  // different species
  arma::cube P_subA = arma::zeros(rows, cols, deltat);
  arma::cube P_subB = arma::zeros(rows, cols, deltat);
  arma::cube P_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Gr_subA = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subB = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Mo_subA = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subB = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube WU_subA = arma::zeros(rows, cols, deltat);
  arma::cube WU_subB = arma::zeros(rows, cols, deltat);
  arma::cube WU_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube qsdA = arma::zeros(rows, cols, time);
  arma::cube runonsdA = arma::zeros(rows, cols, time);
  arma::cube qsdB = arma::zeros(rows, cols, time);
  arma::cube runonsdB = arma::zeros(rows, cols, time);
  arma::cube qsdC = arma::zeros(rows, cols, time);
  arma::cube runonsdC = arma::zeros(rows, cols, time);
  
  arma::cube P_A = arma::zeros(rows, cols, time);
  arma::cube P_B = arma::zeros(rows, cols, time);
  arma::cube P_C = arma::zeros(rows, cols, time);
  
  arma::cube WU_A = arma::zeros(rows, cols, time);
  arma::cube WU_B = arma::zeros(rows, cols, time);
  arma::cube WU_C = arma::zeros(rows, cols, time);
  
  // of all species together
  
  arma::cube P_sub = arma::zeros(rows, cols, deltat);
  arma::cube P = arma::zeros(rows, cols, time);
  
  arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  arma::cube Wu = arma::zeros(rows, cols, time);
  
  
  // actual model loop starts here
  
  for (t = 1; t< (time-1); t++){
    
    for (i=0; i< rows; i++) {
      
      for (j=0; j< cols; j++ ){
        
        //initialise cubes at t= 0
        h(i,j,0) = 10.0;
        // initial condition: homogeneous vegetation cover
        P(i,j,0) = 150.0;
        P_A(i,j,0) = 50.0;
        P_B(i,j,0) = 50.0;
        P_C(i,j,0) = 50.0;
        
        M(i,j,0) = 50.0;
        Svir(i,j,0) = 50.0;
        CM(i,j,0) = 0.0;
        SmI(i,j,0) = 0.0;
        SmM(i,j,0) = 0.0;
        Smh(i,j,0) = 0.0;
        Ch(i,j,0) = 0.0;
        seep(i,j,0) = 0.0;
        
        for (tt = 0; tt< (deltat-1); tt++){
          
          
          if(tt == 0) {
            h_sub(i,j,tt) = h(i,j,t-1);
            P_sub(i,j,tt) = P(i,j,t-1);
            
            P_subA(i,j,tt) = P_A(i,j,t-1);
            P_subB(i,j,tt) = P_B(i,j,t-1);
            P_subC(i,j,tt) = P_C(i,j,t-1);
            
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
            Rain_in =  Rain(t); // 10.0 *
          } else {
            Rain_in = 0.0;
          }
          
          
          
          // adjust infiltration rate
          if(h_sub(i,j,tt) < (timeincr * K_s_in)) {
            alpha_i = 1.0;
          } else {
            alpha_i = 1.0-((h_sub(i,j,tt) - (timeincr * K_s_in))/h_sub(i,j,tt));
          }
          
          
          // Runoff
          q_sub(i,j,tt+1) = timeincr * OF(h_sub(i,j,tt), cn, Mn, (slope(i,j)*0.01));
          
          
          // Runon
          mat runon_store(rows, cols, fill::ones);
          mat ones = mat(rows, cols, fill::ones);
          runon_store = Surface(rows, cols, flowdir_new, write_flowdirTable(), q_sub.slice(tt), ones);
          runon_sub(i,j,tt+1) = runon_store(i,j);
          
          
          // calculate water depth on soil
          h_sub(i,j,tt+1) =  h_sub(i,j,tt) + Rain_in
            - (timeincr * Infil(h_sub(i,j,tt),P_sub(i,j,tt), alpha_i, k_in, W0_in)) - q_sub(i,j,tt+1) + seep_sub(i,j,tt) + runon_sub(i,j,tt+1);
          
          
          //    Ponding check and rewrite of elevation mat
          if (((h_sub(i,j,tt+1) + elev(i,j)*1000.0) > (diff(i,j)*1000.0)) & (diff(i,j) > 0.0)){ // times 1000 for [mm]
            
            mat elev_substitute(rows,cols);
            // elev_substitute = elev;
            elev_substitute(i,j) =  elev(i,j) + diff(i,j);
            
            // calling TauDEM again
            int bla = call_Taudem(elev_substitute);
            // flow directions
            mat flowdir_adj = flowdir_load();
            //slope of each cell
            mat slope_adj = slope_load();
            
            q_sub(i,j,tt+1) = timeincr * OF(h_sub(i,j,tt+1), cn, Mn, slope_adj(i,j));
            
            runon_store = Surface(rows, cols, flowdir_adj, write_flowdirTable(), q_sub.slice(tt+1), ones);
            runon_sub(i,j,tt+1) = runon_store(i,j);
            
            h_sub(i,j,tt+1) =  h_sub(i,j,tt+1) + Rain_in
              - (timeincr * Infil(h_sub(i,j,tt+1),P_sub(i,j,tt+1), alpha_i, k_in, W0_in)) - q_sub(i,j,tt+1) + seep_sub(i,j,tt+1) + runon_sub(i,j,tt+1);
            
          }
          
          //Infiltration
          I_sub(i,j,tt+1) = timeincr * Infil(h_sub(i,j,tt+1), P_sub(i,j,tt) , alpha_i, k_in, W0_in); 
          
          
          // plant water uptake
          WU_subA(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subA(i,j,tt) , gmax_in, k1_inA); 
          WU_subB(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subB(i,j,tt) , gmax_in, k1_inB); 
          WU_subC(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subC(i,j,tt) , gmax_in, k1_inC); 
          
          WU_sub(i,j,tt+1) = WU_subA(i,j,tt+1) + WU_subB(i,j,tt+1) + WU_subC(i,j,tt+1);
          
          // adjustment soil moisture
          M_sub(i,j,tt+1) = M_sub(i,j,tt) + I_sub(i,j,tt+1) - WU_sub(i,j,tt+1);
          
          //           
          // plant growth
          Gr_subA(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subA(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_inA, P0); //, P0, sigmaP, CM_sub(i,j,tt));
          Gr_subB(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subB(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_inB, P0);
          Gr_subC(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subC(i,j,tt), P_sub(i,j,tt), c_in, gmax_in, k1_inC, P0); 
          //  plant Mortality
          Mo_subA(i,j,tt) = timeincr * Mo(P_subA(i,j,tt), d_inA, sigmaPA, CM_sub(i,j,tt));
          Mo_subB(i,j,tt) = timeincr * Mo(P_subB(i,j,tt), d_inB, sigmaPB, CM_sub(i,j,tt));
          Mo_subC(i,j,tt) = timeincr * Mo(P_subC(i,j,tt), d_inC, sigmaPC, CM_sub(i,j,tt));
          
          // seed transport with runoff/runon species A
          mat c2A = timeincr * c02*exp(-sigmaPA*CM_sub.slice(tt));
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c2A(i,j)))
          {
            qsd_subA(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subA(i,j,tt); //(Saco, 2007)
            //qsd_sub(i,j,tt) = c1 * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
          }
          if((q_sub(i,j,tt+1)*c1) > c2A(i,j)){ //(Saco, 2007)
            qsd_subA(i,j,tt) = c2A(i,j) * P_subA(i,j,tt);
          }
          // seed transport with runoff/runon species B
          mat c2B = timeincr * c02*exp(-sigmaPB*CM_sub.slice(tt));
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c2B(i,j)))
          {
            qsd_subB(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subB(i,j,tt); //(Saco, 2007)
            //qsd_sub(i,j,tt) = c1 * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
          }
          if((q_sub(i,j,tt+1)*c1) > c2B(i,j)){ //(Saco, 2007)
            qsd_subB(i,j,tt) = c2B(i,j) * P_subB(i,j,tt);
          }
          // seed transport with runoff/runon species C
          mat c2C = timeincr * c02*exp(-sigmaPC*CM_sub.slice(tt));
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c2C(i,j)))
          {
            qsd_subC(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subC(i,j,tt); //(Saco, 2007)
            //qsd_sub(i,j,tt) = c1 * q_sub(i,j,tt)*P_sub(i,j,tt); //(Saco, 2007)
          }
          if((q_sub(i,j,tt+1)*c1) > c2C(i,j)){ //(Saco, 2007)
            qsd_subC(i,j,tt) = c2C(i,j) * P_subC(i,j,tt);
          }
          
          
          //seed runon Species A
          mat runonsd_storeA(rows, cols, fill::ones);
          runonsd_storeA = Surface(rows, cols, flowdir_new, write_flowdirTable(),qsd_subA.slice(tt),ones);
          runonsd_subA(i,j,tt) = runonsd_storeA(i,j);
          
          //seed runon Species B
          mat runonsd_storeB(rows, cols, fill::ones);
          runonsd_storeB = Surface(rows, cols, flowdir_new, write_flowdirTable(),qsd_subB.slice(tt),ones);
          runonsd_subB(i,j,tt) = runonsd_storeB(i,j);
          
          //seed runon Species C
          mat runonsd_storeC(rows, cols, fill::ones);
          runonsd_storeC = Surface(rows, cols, flowdir_new, write_flowdirTable(),qsd_subC.slice(tt),ones);
          runonsd_subC(i,j,tt) = runonsd_storeC(i,j);
          
          //Germination reduction due to salt
          double germ = (1.0-(M_sub(i,j,tt+1)/Svir_sub(i,j,tt)));
          // 
          // seed diffusion (into all directions due to wind/animals)
          mat seed_diff_loss = seedDiffusionLoss(rows, cols, P_sub.slice(tt),  Dp,  timeincr);
          mat seed_diff_gain = seedDiffusionGain(rows, cols, write_DiffdirectionTable(), P_sub.slice(tt),  Dp,  timeincr);
          
          //  Plant biomass balance
          
          P_subA(i,j,tt+1) = P_subA(i,j,tt) + Gr_subA(i,j,tt) - Mo_subA(i,j,tt) - qsd_subA(i,j,tt) + germ * runonsd_subA(i,j,tt) - seed_diff_loss(i,j) + germ * seed_diff_gain(i,j) 
            + zeta * interference(rows,cols, i,j,P_sub.slice(tt), dx, b1, b2, q1A, q2A);
          P_subB(i,j,tt+1) = P_subB(i,j,tt) + Gr_subB(i,j,tt) - Mo_subB(i,j,tt) - qsd_subB(i,j,tt) + germ * runonsd_subB(i,j,tt) - seed_diff_loss(i,j) + germ * seed_diff_gain(i,j) 
            + zeta * interference(rows,cols, i,j,P_sub.slice(tt), dx, b1, b2, q1B, q2B);
          P_subC(i,j,tt+1) = P_subC(i,j,tt) + Gr_subC(i,j,tt) - Mo_subC(i,j,tt) - qsd_subC(i,j,tt) + germ * runonsd_subC(i,j,tt) - seed_diff_loss(i,j) + germ * seed_diff_gain(i,j) 
            + zeta * interference(rows,cols, i,j,P_sub.slice(tt), dx, b1, b2, q1C, q2C);
          
          P_sub(i,j,tt+1) = P_subA(i,j,tt+1) + P_subB(i,j,tt+1) + P_subC(i,j,tt+1);
          
          // Zr_in averaged from different root depths
          double Zr_in = (P_subA(i,j,tt+1)/P_sub(i,j,tt+1)*Zr_inA) + (P_subB(i,j,tt+1)/P_sub(i,j,tt+1)*Zr_inB) + (P_subC(i,j,tt+1)/P_sub(i,j,tt+1)*Zr_inC);
          
          //vertical water flux (capillary rise/drainage)
          flux_sub(i,j,tt) = L_n(M_sub(i,j,tt+1),Zras(i,j),n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);
          
          //           
          //adjustment for soil mositure balance
          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) +  (timeincr * flux_sub(i,j,tt));
          
          // Seepage
          if(M_sub(i,j,tt+1) > (s_fc_in * n_in * Zr_in)){
            
            seep_sub(i,j,tt) = M_sub(i,j,tt+1)-(s_fc_in * n_in * Zr_in);
            // soil moisture adjustment due to seepage (if condition is met)
            M_sub(i,j,tt+1) = M_sub(i,j,tt+1) - seep_sub(i,j,tt);
          }
          //// Subsurface lateral flow
          mat Subsurfacerunon_store = Subsurface(rows, cols, flowdir_new, write_flowdirTable(), M_sub.slice(tt+1),  ones, Dm, timeincr);
          Subsrunon_sub(i,j,tt) = Subsurfacerunon_store(i,j);
          
          //adjustment for soil mositure balance (complete now)
          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) + Subsrunon_sub(i,j,tt) - (Dm * timeincr * M_sub(i,j,tt+1));
          
          // salt leaching
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
          SmI_sub(i,j,tt+1) = SmI_sub(i,j,tt) + (I_sub(i,j,tt+1) * Ch_sub(i,j,tt));
          //           
          //           // salt mass coming in through lateral subsurface flow
          mat runonSubsSalt_store(rows, cols, fill::ones);
          runonSubsSalt_store = Subsurface(rows, cols, flowdir_new, write_flowdirTable(), M_sub.slice(tt+1), CM_sub.slice(tt), Dm, timeincr);
          runonSubsSalt_sub(i,j,tt) = runonSubsSalt_store(i,j);
          //           
          // #salt mass balance in soil
          SmM_sub(i,j,tt+1) = SmI_sub(i,j,tt+1) + U_salt(i,j,tt) - L_salt(i,j,tt) - (CM_sub(i,j,tt) * seep_sub(i,j,tt)) + runonSubsSalt_sub(i,j,tt) -
            (Dm * timeincr * CM_sub(i,j,tt) * M_sub(i,j,tt+1));
          
          //  salt concentration in soil
          CM_sub(i,j,tt+1) = SmM_sub(i,j,tt+1)/M_sub(i,j,tt+1);
          
          
          
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
          Svir_sub(i,j,tt+1) = n_in * Zr_in * (pow((h1bar_in * 0.1),(1.0/b_in))) *
            pow((h1bar_in * 0.1 * pow((M_sub(i,j,tt+1)/(n_in * Zr_in)),-b_in))+(3.6 * (CM_sub(i,j,tt+1)*(1.0/58.44))),(-1.0/b_in));
          
          // Stress gradient hypothesis (less competition with more salinity)
          double beta = Svir_sub(i,j,tt+1)/M_sub(i,j,tt+1);
          
          b2A = b2A * beta;
          b2B = b2B * beta;
          b2C = b2C * beta;
          
          b1A = 1.0 - b2A;
          b1B = 1.0 - b2B;
          b1C = 1.0 - b2C;
          
          //           
          // surface salt transport with runon
          mat salt_runon_store(rows, cols, fill::ones);
          salt_runon_store = Surface(rows, cols, flowdir_new, write_flowdirTable(),q_sub.slice(tt+1), CM_sub.slice(tt+1));
          salt_runon_sub(i,j,tt) = salt_runon_store(i,j);
          
          //     Balance of salt mass in overland flow depth h
          Smh_sub(i,j,tt+1) = Smh_sub(i,j,tt) + (ConcConst_in * Rain_in) + (CM_sub(i,j,tt+1) * seep_sub(i,j,tt)) - (q_sub(i,j,tt+1)*Ch_sub(i,j,tt))
            +  salt_runon_sub(i,j,tt);
          
          // salt concentration in overland flow depth h
          Ch_sub(i,j,tt+1) = (Smh_sub(i,j,tt+1)/h_sub(i,j,tt+1));
          
          //           
          //// salt balance
          //           
          mb_salt_sub(i,j,tt) = SmI_sub(i,j,tt) + U_salt(i,j,tt) + runonSubsSalt_sub(i,j,tt) - (CM_sub(i,j,tt)*seep_sub(i,j,tt)) -
            L_salt(i,j,tt) - (Dm * timeincr * CM_sub(i,j,tt) * M_sub(i,j,tt));
          // # checking the mass balance
          mb_sub(i,j,tt) = I_sub(i,j,tt) - WU_sub(i,j,tt) + (flux_sub(i,j,tt) * timeincr);
          
        }
        
        // Aggregating the subdaily to daily values
        
        h(i,j,t) = h_sub(i,j,(deltat-2));
        Wu(i,j,t) = WU_sub(i,j,(deltat-2));
        
        WU_A(i,j,t) = WU_subA(i,j,(deltat-2));
        WU_B(i,j,t) = WU_subB(i,j,(deltat-2));
        WU_C(i,j,t) = WU_subC(i,j,(deltat-2));
        
        P(i,j,t) = P_sub(i,j,(deltat-2));
        
        P_A(i,j,t) = P_subA(i,j,(deltat-2));
        P_B(i,j,t) = P_subB(i,j,(deltat-2));
        P_C(i,j,t) = P_subC(i,j,(deltat-2));
        
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
        double summb_salt = 0.0;
        
        double sumqsdA = 0.0;
        double sumqsdB = 0.0;
        double sumqsdC = 0.0;
        
        double sumrunonsdA = 0.0;
        double sumrunonsdB = 0.0;
        double sumrunonsdC = 0.0;
        
        double sumrunonSubsSalt = 0.0;
        double sumsalt_runon = 0.0;
        double sumSubsrunon = 0.0;
        
        for(int tt = 0; tt < deltat; tt++)
        {
          sumI += I_sub(i,j,tt);
          sumrunon += runon_sub(i,j,tt);
          sumq += q_sub(i,j,tt);
          sumflux += flux_sub(i,j,tt) * timeincr;
          summb += mb_sub(i,j,tt);
          summb_salt += mb_salt_sub(i,j,tt);
          
          sumqsdA += qsd_subA(i,j,tt);
          sumrunonsdA += runonsd_subA(i,j,tt);
          
          sumqsdB += qsd_subB(i,j,tt);
          sumrunonsdB += runonsd_subB(i,j,tt);
          
          sumqsdC += qsd_subC(i,j,tt);
          sumrunonsdC += runonsd_subC(i,j,tt);
          
          sumrunonSubsSalt += runonSubsSalt_sub(i,j,tt);
          sumsalt_runon += salt_runon_sub(i,j,tt);
          sumSubsrunon += Subsrunon_sub(i,j,tt);
          
        }
        
        q(i,j,t) = sumq;
        runon(i,j,t) = sumrunon;
        In(i,j,t) = sumI;
        flux(i,j,t) = sumflux;
        mb(i,j,t) = summb;
        mb_salt(i,j,t) = summb_salt;
        
        qsdA(i,j,t) = sumqsdA;
        runonsdA(i,j,t) = sumrunonsdA;
        
        qsdB(i,j,t) = sumqsdB;
        runonsdB(i,j,t) = sumrunonsdB;
        
        qsdC(i,j,t) = sumqsdC;
        runonsdC(i,j,t) = sumrunonsdC;
        
        runonSubsSalt(i,j,t) = sumrunonSubsSalt;
        salt_runon(i,j,t) = sumsalt_runon;
        Subsrunon(i,j,t) = sumSubsrunon;
        
      }
    }
  }
  
  // Creating output list
  
  List fields;
  arma::field<arma::cube> f1(32);
  f1( 0 ) = h;
  f1( 1 ) = q;
  f1( 2 ) = In;
  f1( 3 ) = runon;
  f1( 4 ) = Wu;
  f1( 5 ) = WU_A;
  f1( 6 ) = WU_B;
  f1( 7 ) = WU_C;
  f1( 8 ) = P;
  f1( 9 ) = P_A;
  f1( 10 ) = P_B;
  f1( 11 ) = P_C;
  f1( 12 ) = flux;
  f1( 13 ) = M;
  f1( 14 ) = SmM;
  f1( 15 ) = CM;
  f1( 16 ) = mb;
  f1( 17 ) = Svir;
  f1( 18 ) = SmI;
  f1( 19 ) = Smh;
  f1( 20 ) = Ch;
  f1( 21 ) = qsdA;
  f1( 22 ) = qsdB;
  f1( 23 ) = qsdC;
  f1( 24 ) = runonsdA;
  f1( 25 ) = runonsdB;
  f1( 26 ) = runonsdC;
  f1( 27 ) = seep;
  f1( 28 ) = runonSubsSalt;
  f1( 29 ) = salt_runon;
  f1( 30 ) = Subsrunon;
  f1( 31 ) = mb_salt;
  fields["field<cube>"] = f1;
  
  
  List output = List::create(_["fields : field<cube>"] = f1);
  
  return output;
  
}

/*** R
# results <- SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,saltpar = saltpar1,
#                                  dims = list(rows=rows,cols=cols,time=time, Z=Z),  Rain=Rain,
#                                  diverseInput = list (deltat = deltat, gslp = gslp, ext=ext, Zr=Zr, Dm = Dm, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, sigmaP= sigmaP,
#                                                       c1 = c1, c02 = c02, Dp = Dp, b1 = b1, b2 = b2, q1 = q1, q2 = q2, zeta = zeta))
#   library(rasterVis)
#   coul = brewer.pal(8, "YlGn")
#   coul = colorRampPalette(coul)(100)
#   
#   
#   qr<-brick(results$fields[[6]][2:19,2:19,1:100])
#   
#   levelplot(qr,main="P [g/m^2] ",sub="day 1 to day 50, salt from gw, randomly varied alpha and lambda", col.regions = coul) #col.regions = YlGn.colors(20))
# # 
# results$fields[[6]][2:19,2:19,90:91]
# # results$fields[[6]][2:19,2:19,94]
# # results$fields[[6]][2:19,2:19,1:10]
# # results$fields[[6]][2:19,2:39,100]
# # results$fields[[1]][2:19,2:19,280]
# # results$fields[[6]][1:20,1:20,1:20]
#  animate(qr, n=1)
# # # 
# # results$fields[[2]]
# results$fields[[6]]
# f1( 0 ) = h;
# f1( 1 ) = q;
# f1( 2 ) = In;
# f1( 3 ) = runon;
# f1( 4 ) = Wu;
# f1( 5 ) = P;
# f1( 6 ) = flux;
# f1( 7 ) = M;
# f1( 8 ) = SmM;
# f1( 9 ) = CM;
# f1( 10 ) = mb;
# f1( 11 ) = Svir;
# f1( 12 ) = SmI;
# f1( 13 ) = Smh;
# f1( 14 ) = Ch;
# f1( 15 ) = qsd;
# f1( 16 ) = runonsd;
# f1( 17 ) = seep;
# f1( 18 ) = runonSubsSalt;
# f1( 19 ) = salt_runon;
# f1( 20 ) = Subsrunon;
# f1( 21 ) = mb_salt;

*/