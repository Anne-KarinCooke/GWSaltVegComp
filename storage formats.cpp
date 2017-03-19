// Storage formats, arrays etc

int const rows = 10;
int const cols = 10;
int const time = 10;
int const deltat =12;

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
          //M sub
          double M_sub[rows][cols][deltat]= {{0}};
          //P sub
          double P_sub[rows][cols][deltat]= {{0}};
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

            
///       
                  double h_old[rows][cols]={{0}};
                  double P_old[rows][cols]={{0}};
                  double M_old[rows][cols]={{0}};
                  double SmI_old[rows][cols]={{0}};
                  double CM_old[rows][cols]={{0}};
                  double Svir_old[rows][cols]={{0}};

     
    
                      
                      
                
                      
                      
                      