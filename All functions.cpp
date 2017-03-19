#include <Rcpp.h>
using namespace Rcpp;

//Sandy Clay Loam


// Rain_function<-function(time){ Rain <- rep(1, time)
//   return(Rain)}
// Rain<-Rain_function(time=time)
  
  
  
int const rows = 10;
int const cols = 10;
int const time = 10;
int const deltat =12;
  
List Constants_cpp(List soilpar,
                   List saltpar,
                   List vegpar) 
{
  
  
  // SOIL 
  
  
  const double n = soilpar["n"];
  const double K_s = soilpar["K_s"];
  const double b = soilpar["b"];
  const double s_fc = soilpar["s_fc"];
  const double psi_s_bar = soilpar["psi_s_bar"];
  const double h1bar = soilpar["h1bar"];
  const double hb = soilpar["hb"];
  const double Mn = soilpar["Mn"];
  const double cn = soilpar["cn"];
  double alpha_i = soilpar["alpha_i"]; 
  
  soilpar["n"] = 0.367; //porosity
  soilpar["K_s"] = 52.08*10; //Hydraulic conductivity mm/day 
  soilpar["b"] = 6.4069; // campbell's b
  soilpar["s_fc"] = 0.2677/n; // Field capacity
  soilpar["psi_s_bar"] = -1.2E-3; // bubbling pressure
  soilpar["h1bar"] =  -psi_s_bar;
  soilpar["hb"] = psi_s_bar*pow(-10,5);//  mm
  soilpar["Mn"] = 10; // Manning's n
  soilpar["cn"] = 1; // runoff conversion factor
  soilpar["alpha_i"] = 1;//#maximum infiltration rate per day, This needs to be a fraction of h (p117 Saco and Moreno-Las Heras) 
  
  
  
  
  // VEGETATION
  const double Zr = vegpar["Zr"];
  const double k = vegpar["k"];
  const double W0 = vegpar["W0"];
  const double gmax = vegpar["gmax"];
  const double c = vegpar["c"];
  const double k1 = vegpar["k1"];
  const double d = vegpar["d"]; //fraction of plant mortality
  
  vegpar["Zr"] = 400; //mm, Grass
  vegpar["d"] = 0.24;//Saco et al, 2013
  vegpar["c"] = 0.10;//Saco et al, 2013
  vegpar["k1"] = 5;//Saco et al, 2013
  vegpar["gmax"] = 0.05;//Saco et al, 2013
  vegpar["W0"] = 0.2;//Saco et al, 2013
  vegpar["k"] = 12;//Saco et al, 2013
  
  
  
  
  // SALT                            
  const double ConcConst = saltpar["ConcConst"];
  const double CMgw = saltpar["CMgw"];
  const double f = saltpar["f"];
  
  saltpar["ConcConst"] = 0;  //ConcConst is the concentration of the salt in the infiltrating water in g/l
  saltpar["CMgw"] = 0; //CMgw is the goundwater salt concentration  in g/l
  saltpar["f"] = 0.8; //f is the soil salt leaching efficiency (whether some salt is retained)
}



// Storage arrays for the daily time steps declared and initialized

// soil moisture [mm]
                      double M[rows][cols][time] = {{0}};
                      
                      //  plant biomass density
                      double P[rows][cols][time]= {{0}};
                      //h
                      double h[rows][cols][time]= {{0}};
                      //CM
                      double CM[rows][cols][time]= {{0}};
                      // SmI
                      double SmI[rows][cols][time]= {{0}};
                      // SmM
                      double SmM[rows][cols][time]= {{0}};
                      // In
                      double In[rows][cols][time]= {{0}};
                      //Svir
                      double Svir[rows][cols][time]= {{0}};
                      // flux
                      double flux[rows][cols][time]= {{0}};
                      //q
                      double q[rows][cols][time]= {{0}};
                      //runon
                      double runon[rows][cols][time]= {{0}};
                      //mb
                      double mb[rows][cols][time]= {{0}};

//
// Storage arrays for the SUBdaily time steps declared and initialized
                                    
                                    //h sub
                                    double h_sub[rows][cols][deltat]= {{0}};
                                    //P sub
                                    double P_sub[rows][cols][deltat]= {{0}};
                                    //M sub
                                    double M_sub[rows][cols][deltat]= {{0}};
                                    //CM sub
                                    double CM_sub[rows][cols][deltat]= {{0}};
                                    // SmIsub
                                    double SmI_sub[rows][cols][deltat]= {{0}};
                                    // SmMsub
                                    double SmM_sub[rows][cols][deltat]= {{0}};
                                    // Isub
                                    double I_sub[rows][cols][deltat]= {{0}};
                                    //Svirsub
                                    double Svir_sub[rows][cols][deltat]= {{0}};
                                    // fluxsub
                                    double flux_sub[rows][cols][deltat]= {{0}};
                                    //qsubsub
                                    double q_sub[rows][cols][deltat]= {{0}};
                                    //runonsub
                                    double runon_sub[rows][cols][deltat]= {{0}};
                                    //mbsub
                                    double mb_sub[rows][cols][deltat]= {{0}};
                                    //Gr_sub
                                    double Gr_sub[rows][cols][deltat]= {{0}};
                                    //Mo_sub
                                    double Mo_sub[rows][cols][deltat]= {{0}};
                                    //WU_sub
                                    double WU_sub[rows][cols][deltat]= {{0}};

                                    // Salt leaching 
                                    double L_salt[rows][cols][deltat]= {{0}};
                                    //salt rise
                                    double U_salt[rows][cols][deltat]= {{0}};
                                    
                                    
                                                              ///   OLDS    
                                                              double h_old[rows][cols]={{0}};
                                                              double P_old[rows][cols]={{0}};
                                                              double M_old[rows][cols]={{0}};
                                                              double SmI_old[rows][cols]={{0}};
                                                              double CM_old[rows][cols]={{0}};
                                                              double Svir_old[rows][cols]={{0}};

// FUNCTIONS


        // Overland flow RUNOFF, kinematic wave approach, as used in from Saco et al 2013
        
        double OF(double h, const double soilpar["cn"], double slope){
          
          double qq = (soilpar["cn"]/soilpar["Mn"])*(pow(h,1.666667))*sqrt(slope));
          return qq;
        }
        

//vertical water flux function (capillary rise and drainage), eq from Salvucci 1993

        double L_n(double M, double Z, double n, double Zr, double b, double K_s, double psi_s_bar){
          
          double hb = psi_s_bar*pow(10,5);
          double s=M/(n*Zr); // extract n and Zr from list and define them
          double psi = hb*pow(s,-b);
          double s_fc = pow((Z/hb),(-1/b));// define hb and b
          double m = 2 + 3/b;
          double qf = pow((Z/hb),(-m))-(pow((psi/hb),-m)/(1+pow((psi/hb),-m)))+(m-1)*pow((Z/hb),-m);
          double flux = K_s * qf;
          
          return flux;
          
        }
        
// VEGETATION FUNCTIONS
        
        //Plant water uptake function WU
        
        double WU(double M, double P, double gmax, double k1 ) {  /// not quite happy with the list item calling...list par
          
          double Wu = gmax*(M/(M+k1))*P;   // function is called WU, output is called Wu, cannot be the same
          return Wu;
          
        }
        
        
        //Plant Growth function Gr
        
        double Gr(double M, double P, double c, double gmax, double k1){
          
          double Gro = c*WU(M,P,gmax,k1);
          return Gro;
        }
        
        
        // Plant mortality function Mo
        
        double Mo(double P, double M, double Svir, double d ){
          
          double Mort=P*(d*(M/Svir));
          return Mort;
        }  
        
        
        // Infiltration function Infil
        
        double Infil(double h, double P, double alpha_i, double k, double W0){
          
          double I=alpha_i*h*(P+k*W0)/(P+k);
          return I;
        }
        
        


        
       //   source("Rasters.R")
        //  source("Rainfall.R")
     

     //// BIG OLD BALANCES FUNCTION 
     
#include <array>
#include <iostream>
#include <tuple>
using namespace std;
     
     // std::array <double, 3> myarray; // declare an integer array with length 3

double balances2D(double Rain, List par, List soilpar, List vegpar){
  
  int i;
  int j;
  int t;
  int tt;
  int deltat;
  int rain;
  int rows;
  int cols;
  //  int arr[rows * cols];
  int time;
  double timeincr = 1/deltat;
  
  
  for (i=1; i< rows; i++) {
    
    for (j=1; j< cols; j++ ){
      
      for (t=2; t< rain; t++){  // it is actually the length of rain, has to be defined, Rain has to be rain ?
        
        for (tt=1; tt< (deltat-1); tt++){  // deltat has to be defined
          
          // the OLDs
          // the following DOES NOT work this way, array value assignment is different in Cpp
// 
//           void revalue(double r, double bucky[], int n){
//             for(int i=0;i<n;i++){
//                 bucky[i] *=r;
//             }
//           })
          

          
          // cin >> h_old[i];
          
          if(tt == 1) {
            h_old[i][j] = h[i][j][t-1];
             } else {
            h_old[i][j] =  h_sub[i][j][tt];
              }
              if(tt == 1) {
                P_old[i][j] = P[i][j][t-1];
              } else {
                P_old[i][j] =  P_sub[i][j][tt];
              }
                if(tt == 1) {
                  M_old[i][j]= M[i][j][t-1];
                } else {
                  M_old[i][j] =  M_sub[i][j][tt];
                }
                    if(tt == 1) {
                      SmI_old[i][j]= SmI[i][j][t-1];
                    } else {
                      SmI_old[i][j] =  SmI_sub[i][j][tt];
                    }
                        if(tt == 1) {
                          CM_old[i][j]= CM[i][j][t-1];
                        } else {
                          CM_old[i][j] =  CM_sub[i][j][tt];
                        }
                            if(tt == 1) {
                              Svir_old[i][j]= Svir[i][j][t-1];
                            } else {
                              Svir_old[i][j] = Svir_sub[i][j][tt];
                            }
          
          double slope; // think about this again
          
          q_sub[i][j][tt+1] = OF(h_old[i][j], soilpar["cn"], soilpar["Mn"], slope)*timeincr; // define slope!!
          
          
      
          
          //      runon_sub[i,j,tt+1] <-rn[i,j]*q_sub[i,j,tt]
          
          //      h_sub[i,j,tt+1] <- h.old[i,j] + ifelse(tt==1,(10*Rain[t]),0) - Infil(h.old[i,j], P.old[i,j],par)*timeincr - q_sub[i,j,tt] + runon_sub[i,j,tt]
          
          //    par$alpha_i <- ifelse(h_sub[i,j,tt+1]<soilpar$K_s*timeincr, 1,(1-(h_sub[i,j,tt+1]-soilpar$K_s*timeincr)/h_sub[i,j,tt+1]))
          //       I_sub[i,j,tt+1] <- Infil(h.old[i,j], P.old[i,j],par)*timeincr
          // Water uptake           
          WU_sub[i][j][tt] = WU(M_sub[i][j][tt+1],P_old[i][j],soilpar["gmax"], soilpar["k1"])*timeincr;
          // Growth            
          Gr_sub[i][j][tt] = Gr(Svir_old[i][j], P_old[i][j], soilpar["c"], soilpar["gmax"], soilpar["k1"])*timeincr; 
          //Mortality
          Mo_sub[i][j][tt] = Mo(P_old[i][j], M_old[i][j], Svir_old[i][j],d)*timeincr;
          // Plant biomass balance             
          P_sub[i][j][tt+1] = P_old[i][j] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; /// not sure if this all is ok this way or too close to R
          
          // Water balance before drainage
          M_sub[i][j][tt+1] = M_old[i][j] + I_sub[i][j][tt] - WU_sub[i][j][tt];
          
          // Drainage/Capillary rise (vertical water flux)          
          
          flux_sub[i][j][tt+1] = L_n(M_sub[i][j][tt+1],Zras[i][j],soilpar=soilpar,vegpar=vegpar));  /// how yo read in Zras,,, change the soilpar and vegpar stuff
          
          // Adjustment for M including flux
          //
          M_sub[i][j][tt+1] = M_sub[i][j][tt+1] +  flux_sub[i][j][tt+1]*timeincr; 

          
          // Salt balance
          
          
          // salt leaching
          // 
          if(flux_sub[i][j][tt+1] <0 ) {
            L_salt[i][j][tt+1] = soilpar["f"] *CM_sub[i][j][tt+1] * flux_sub[i][j][tt+1]*timeincr;
          } else {
            L_salt[i][j][tt+1] = 0;
          }
            
            // salt upflow
            if(flux_sub[i][j][tt+1] >0 ) {
              U_salt[i][j][tt+1] = soilpar["f"] * CM_sub[i][j][tt+1] * flux_sub[i][j][tt+1]*timeincr;
            } else {
              L_salt[i][j][tt+1] = 0;
            }
            
            
// # salt mass coming in with infiltration
                  SmI_sub[i][j][tt+1] = SmI_old[i][j] + I_sub[i][j][tt] * soilpar["ConcConst"];
                    
// #salt mass in soil
                    SmM_sub[i][j][tt+1] = SmI_sub[i][j][tt+1] + U_salt[i][j][tt+1] - L_salt[i][j][tt+1];
                    
// # salt concentration in soil
                    CM_sub[i][j][tt+1] = (SmM_sub[i][j][tt+1]/M_sub[i][j][tt+1])*(1/58.44);         
                      
// # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
                      Svir_sub[i][j][tt+1] = soilpar['n']* vegpar['Zr']*(pow((soilpar['h1bar']* pow(10,-1)),(1/soilpar$b))) *
                        
                        
                        ((soilpar$h1bar*10^-1)*(M_sub[i][j][tt+1]/(soilpar$n*vegpar$Zr))^(-soilpar$b)+(3.6*CM_sub[i][j][tt+1]))^(-1/soilpar$b)
                                                                          
//  Svir_sub[i,j,tt+1]<-soilpar$n*vegpar$Zr*((soilpar$h1bar*10^-1)^(1/soilpar$b))*
// ((soilpar$h1bar*10^-1)*(M_sub[i,j,tt+1]/(soilpar$n*vegpar$Zr))^(-soilpar$b)+(3.6*CM_sub[i,j,tt+1]))^(-1/soilpar$b)
//                         
// # checking the mass balance!
                        mb_sub[i][j][tt] = I_sub[i][j][tt] - WU_sub[i][j][tt] + flux_sub[i][j][tt]*timeincr;  
                          
                          
                          
        }
// # Aggregating the substep results to daily values.
        
        P[i][j][t] = P_sub[i][j][deltat];
        M[i][j][t] = M_sub[i][j][deltat];
        h[i][j][t] = h_sub[i][j][deltat];
        CM[i][j][t] = CM_sub[i][j][deltat];
        SmI[i][j][t] = SmI_sub[i][j][deltat];
        SmM[i][j][t] = SmM_sub[i][j][deltat];
        Svir[i][j][t] = Svir_sub[i][j][deltat];
        // int main()
        // {
        // 
          
          double sumI = std::accumulate(I_sub[i][j][1], I_sub[i][j][deltat], 0);
          In[i][j][t] = sumI;
          
          double sumflux = std::accumulate(flux_sub[i][j][1], flux_sub[i][j][deltat], 0);
          flux[i][j][t]= sumflux; 
          
          double sumq = std::accumulate(q_sub[i][j][1], q_sub[i][j][deltat], 0);
          q[i][j][t] = sumq;
          
          double sumrunon = std::accumulate(runon_sub[i][j][1], runon_sub[i][j][deltat], 0);
          runon[i][j][t] = sumq;
          
          double summb = std::accumulate(mb_sub[i][j][1], mb_sub[i][j][deltat], 0);
          mb[i][j][t] = summb;

        
   
            
            
        }
      }
    }
  }
  // too many brackets?! 
  
}

Out <- list(P=P[,,],M=M[,,],h=h[,,], CM=CM[,,], SmM=SmM[,,], In=In[,,], flux=flux[,,], Svir=Svir[,,],h=h[,,], q=q[,,],mb=mb[,,], runon=runon[,,])
return(Out)

}





        
#// bits and pieces
#// int arr[rows * cols];
 
                      
                      
                   






