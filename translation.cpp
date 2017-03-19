#include <Rcpp.h>
using namespace Rcpp;

// Balances function 

#include <iostream>  /// why do I need this, read up on it
using namespace std;// what does this mean ?

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
  double timeinc = 1/deltat;

  
  for (i=1; i< rows; ) {
    
    for (j=1; j< cols; ){
      
      for (t=2; t< rain;){  // it is actually the length of rain, has to be defined, Rain has to be rain ?
        
        for (tt=1; tt< (deltat-1);){  // deltat has to be defined
          
          // h.old[i,j] <- ifelse(tt==1,h[i,j,t-1],h_sub[i,j,tt]) 
          // P.old[i,j] <- ifelse(tt==1,P[i,j,t-1],P_sub[i,j,tt])  
          // M.old[i,j] <- ifelse(tt==1,M[i,j,t-1],M_sub[i,j,tt])
          // SmI.old[i,j] <-ifelse(tt==1,SmI[i,j,t-1],SmI_sub[i,j,tt])
          // CM.old[i,j] <-ifelse(tt==1,CM[i,j,t-1],CM_sub[i,j,tt])
          // Svir.old[i,j] <-ifelse(tt==1,Svir[i,j,t-1],Svir_sub[i,j,tt])

          double q_sub[rows][cols][time]; // how to define array and where
          
      
         
          q_sub[i][j][tt+1] = OF(h, cn, Mn, slope)*timeincr;
            
          
          
          An element in 2-dimensional array is accessed by using the subscripts, i.e., row index and column index of the array. For example:
            
            int val = a[2][3];
          
    //      runon_sub[i,j,tt+1] <-rn[i,j]*q_sub[i,j,tt]

    //      h_sub[i,j,tt+1] <- h.old[i,j] + ifelse(tt==1,(10*Rain[t]),0) - Infil(h.old[i,j], P.old[i,j],par)*timeincr - q_sub[i,j,tt] + runon_sub[i,j,tt]

      //    par$alpha_i <- ifelse(h_sub[i,j,tt+1]<soilpar$K_s*timeincr, 1,(1-(h_sub[i,j,tt+1]-soilpar$K_s*timeincr)/h_sub[i,j,tt+1]))
     //       I_sub[i,j,tt+1] <- Infil(h.old[i,j], P.old[i,j],par)*timeincr
// Water uptake           
            WU_sub[i][j][tt] = WU(M_sub[i][j][tt+1],P.old[i][j],gmax, k1)*timeincr;
// Growth            
            Gr_sub[i][j][tt] = Gr(Svir.old[i][j], P.old[i][j], c, gmax, k1)*timeincr; 
//Mortality
            Mo_sub[i][j][tt] = Mo(P.old[i][j], M.old[i][j], Svir.old[i][j],d)*timeincr;
// Plant biomass balance             
            P_sub[i][j][tt+1] = P.old[i][j] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt]; /// not sure if this all is ok this way or too close to R
          
// Water balance before drainage
            M_sub[i][j][tt+1] = M.old[i][j] + I_sub[i][j][tt] - WU_sub[i][j][tt];
          
// Drainage/Capillary rise (vertical water flux)          

         //  flux_sub[i][j][tt+1] = do.call(L_n,list(M=M_sub[i,j,tt+1],Z=Zras[i,j],soilpar=soilpar,vegpar=vegpar))
            
// Adjustment for M including flux
//            M_sub[i,j,tt+1] <- M_sub[i,j,tt+1] + flux_sub[i,j,tt+1]*timeincr #- Diff[i,j,tt] 
            

            
            
// Salt balance
            
            
// salt leaching
            L_salt[i,j,tt+1] <- ifelse(flux_sub[i,j,tt+1]<0, par$f*CM_sub[i,j,tt+1]*flux_sub[i,j,tt+1]*timeincr,0) # leaching of salt
              
// salt upflow
              U_salt[i,j,tt+1] <- ifelse(flux_sub[i,j,tt+1]>0, par$CM.gw*flux_sub[i,j,tt+1]*timeincr,0) # rise of salt
          
          
            
                        }
                    }
                }
            }
  
         }


  
  time <- 10

deltat<-12
Rain_function<-function(time){ Rain <- rep(1, time)
  return(Rain)}
Rain<-Rain_function(time=time)
  
## Source functions
  setwd("H:/Thesis project model/R project/GWSaltVegComp")
    source("Rasters.R")
    source("Rainfall.R")
    source("Infiltration.R")
    source("Flux.R")
    source("Vegetation functions.R")
    source("Runoff.R")
    source("Constants.R")
    source("storage formats.R")
    
    
    
#################################################################################################################
### DISCRETIZATION
    
    
    timeincr= 1/deltat
      
      M[,,1] <- 5
    h[,,1] <- 10
    P[,,1] <- 10
    CM[,,1] <- 0
    Svir[,,1] <- s_fc
      
#################################################################################################################
      
      
##################################################################################################################
##################################################################################################################
#  BALANCES FUNTION STARTS
      
      cppFunction(' 
                    
                    ')    
        
#// bits and pieces
#// int arr[rows * cols];
        
   
          
    
# re-calculate water balance
# 2. before leaching
            M_sub[i,j,tt+1] <- M.old[i,j] + I_sub[i,j,tt] - WU_sub[i,j,tt] 
            
            
            
# 3. calculate leaching and capillary rise amount
            flux_sub[i,j,tt+1]<-do.call(L_n,list(M=M_sub[i,j,tt+1],Z=Zras[i,j],soilpar=soilpar,vegpar=vegpar))
              
              
# Divergence, soil moisture diffusion, still to be fixed
##Diff[i,j,tt]<- M_sub[i,j,tt+1]*(Dm*timeincr)*((Z[g+1]-Z[i,j])/dist)  
              
              
# 4. final adjust soil moisture for leaching/rise AND DIVERGENCE
              M_sub[i,j,tt+1] <- M_sub[i,j,tt+1] + flux_sub[i,j,tt+1]*timeincr #- Diff[i,j,tt] 
              
              
# calculate saltbalance
              
              
# Salt leaching
              L_salt[i,j,tt+1] <- ifelse(flux_sub[i,j,tt+1]<0, par$f*CM_sub[i,j,tt+1]*flux_sub[i,j,tt+1]*timeincr,0) # leaching of salt
                
# salt uplfow
                U_salt[i,j,tt+1] <- ifelse(flux_sub[i,j,tt+1]>0, par$CM.gw*flux_sub[i,j,tt+1]*timeincr,0) # rise of salt
                  
# salt mass coming in with infiltration
                  SmI_sub[i,j,tt+1]<- SmI.old[i,j] + I_sub[i,j,tt]*par$ConcConst 
                    
#salt mass in soil
                    SmM_sub[i,j,tt+1] <- SmI_sub[i,j,tt+1] + U_salt[i,j,tt+1] - L_salt[i,j,tt+1]
                    
# salt concentration in soil
                    CM_sub[i,j,tt+1]<- (SmM_sub[i,j,tt+1]/M_sub[i,j,tt+1])*(1/58.44)         
                      
# Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
                      Svir_sub[i,j,tt+1]<-soilpar$n*vegpar$Zr*((soilpar$h1bar*10^-1)^(1/soilpar$b))*
                        ((soilpar$h1bar*10^-1)*(M_sub[i,j,tt+1]/(soilpar$n*vegpar$Zr))^(-soilpar$b)+(3.6*CM_sub[i,j,tt+1]))^(-1/soilpar$b)
                        
# checking the mass balance!
                        mb_sub[i,j,tt] <- I_sub[i,j,tt] - WU_sub[i,j,tt] + flux_sub[i,j,tt]*timeincr  
                          
                          
                          
                          }
# Aggregating the substep results to daily values.
                          
                          P[i,j,t] = P_sub[i,j,deltat]
                          M[i,j,t] = M_sub[i,j,deltat]
                          h[i,j,t] = h_sub[i,j,deltat] #### modified
                          CM[i,j,t] = CM_sub[i,j,deltat]
#      SmM[t,g] = SmI[t,g] = SmM_sub[deltat,g]  ###
                          SmI[i,j,t] = SmI_sub[i,j,deltat]
                          SmM[i,j,t] = SmM_sub[i,j,deltat]
                          In[i,j,t]= sum(I_sub[i,j,])
                          Svir[i,j,t] = Svir_sub[i,j,deltat]
                          flux[i,j,t]= sum(flux_sub[i,j,])
                          q[i,j,t] = sum(q_sub[i,j,]) ####modified
                          
                          runon[i,j,t] = sum(runon_sub[i,j,])
                          mb[i,j,t] = sum(mb_sub[i,j,])
                          
                          }
                          }  
                          
                          
                          }
                          
                          
                          
                          Out <- list(P=P[,,],M=M[,,],h=h[,,], CM=CM[,,], SmM=SmM[,,], In=In[,,], flux=flux[,,], Svir=Svir[,,],h=h[,,], q=q[,,],mb=mb[,,], runon=runon[,,])
                          return(Out)
                          
                          }
                          
                          
                          
 