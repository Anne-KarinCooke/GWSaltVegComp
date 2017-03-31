#world lonlat raster
r <- raster(ncol=10,nrow=10)
r[] <- 1
r[48] <- 2
r[66:68] <- 3
d <- gridDistance(r,origin=2,omit=3) 
plot(d)

library(gdistance)


r   <- raster(nrows=100,ncols=100,xmn=0,ymn=0,xmx=100,ymx=100)
r[] <- rep(1, ncell(r))
plot(r)
h16  <- transition(r, transitionFunction=function(x){1},16,symm=FALSE)
h16   <- geoCorrection(h16, scl=FALSE)
A <- c(20,50)
h16.acc <- accCost(h16,A)
plot(h16.acc)

Pras <- raster(P[,,tt])

xx <-accCost(geoCorrect(transition(r, transitionFunction=function(x){1},16,symm=FALSE), scl=FALSE),Pras[i,j])
  
  
# values(h16.acc)



## 
for(every cell with species A){
  


P_sub[i,j,tt+1] <- P.old[i,j] + Gr_sub[i,j,tt]- Mo_sub[i,j,tt]
+ integrate(kernel(gdistance))*P_sub[i+1,j+1,tt]

b1*exp(-gdistance*gdistance) - b2*(-gdistance*gdistance)

density(P_sub[i,j,tt],kernel=c("gaussian"))

}






# Plant water uptake
WU <- function(M,P,par){ 
  # using Svir in here means scaling Svir back to M, easier to do at Svir in balances  
  #  WU=par$gmax*((M*(1+Svir))/(((M*(1+Svir))+par$k1)))*P 
  WU=par$gmax*(M/((M+par$k1)))*P 
  
  return(WU)
}

#Plant Growth

Gr <- function(M,P,par) { 
  
  Gr = par$c*WU(M,P,par)
  
  return(Gr)
}


## Plant mortality function WITH SALT INFLUENCE, WITH VIRTUAL SATURATION

Mo <- function(P,M,Svir,par) {
  # needs to be M/Svir because both are "large" numbers
  # you want a number ~1 for multiplication, or <0.1 for addition
  Mo = P*(par$d*(M/Svir))
  
  return(Mo)
  
}

x <- c(0, 1, 1.1, 1.5, 1.9, 2.8, 2.9, 3.5)
 n <- length(x)

xgrid <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01)

h <- 0.4
 bumps <- sapply(x, function(a) gauss((xgrid - a)/h)/(n * h))
 
 
 
 
 
 # Plant water uptake
 WU <- function(M,P,par){ 
   
   WU=par$gmax*(M/((M+par$k1)))*P 
   
   return(WU)
 }
 
 #Plant Growth
 
 Gr <- function(M,P,par){ 
   
   Gr = par$c*WU(M,P,par)
   
   return(Gr)
 }
 
 
 ## Plant mortality function WITH SALT INFLUENCE, WITH VIRTUAL SATURATION
 
 Mo <- function(P,M,Svir,par) {
   # needs to be M/Svir because both are "large" numbers
   # you want a number ~1 for multiplication, or <0.1 for addition
   Mo = P*(par$d*(M/Svir))
   
   return(Mo)
   
 }
 
 # Biomass density species A
 P_subA <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) 
 # Biomass density species B
 P_subB <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  
 # Biomass density species C
 P_subC <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) 
 
 # Growth of biomass species A
 Gr_subA <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  
 # Growth of biomass species B
 Gr_subB <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 # Growth of biomass species C
 Gr_subC <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 
 # Mortality of biomass species A
 Mo_subA<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  
 # Mortality of biomass species B
 Mo_subB<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 # Mortality of biomass species C
 Mo_subC<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 
 # Water uptake in mm species A
 WU_subA<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  
 # Water uptake in mm species B
 WU_subB<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 # Water uptake in mm species C
 WU_subC<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
 
 
 # Water uptake, sum of WU_subA,B,C
 
 WU_sub[i,j,tt] <- WU_subA[i,j,tt] + WU_subB[i,j,tt] + WU_subC[i,j,tt]
 
 WU_sub[i,j,tt] <- WU(M_sub[i,j,tt+1],P.old[i,j],vegpar)*timeincr 
 
 # growth rate
 Gr_sub[i,j,tt] <- Gr(M=Svir.old[i,j], P.old[i,j],par)*timeincr 
 # Mortality
 Mo_sub[i,j,tt]<- Mo(P.old[i,j], M=M.old[i,j], Svir=Svir.old[i,j], par)*timeincr
 
 # calculate plant biomass balance
 P_sub[i,j,tt+1] <- P.old[i,j] + Gr_sub[i,j,tt]- Mo_sub[i,j,tt]
 
 
 
 ##############################################################################
 ## Response of 3 diff functional plant types (Maas and Hoffman, 1977)
 ## (neither temporal threshold nor water logging considered yet )
 ##############################################################################
 
 ### Species A Halophyte
 
   # lower threshold
   th_lowA =
   # low
   conc_lowA =
   # high
   conc_highA =
   # upper threshold
   th_upA =
     
 ### Species B salt-TOLERANT NON-Halophyte
   
   # lower threshold -> irrelevant for this functional type
   # th_low =
   # # low
   # conc_low =
   # high
   conc_highB = 
   # upper threshold
   th_upB =

### Species C salt-SENSITIVE NON-Halophyte
   # lower threshold
   th_lowC =
   # low
   conc_lowC =
   # high
   conc_highC =
   # upper threshold
   th_upC =

   
   
   
### HALOPHYTE
   if (CM_sub[i,j,tt] >= conc_lowA){
     Gr_subA = ## INCREASED by some factor
       }
 
         if (CM_sub[i,j,tt] >= conc_highA){
         Gr_subC = ##  reduced by some factor
         }
 
           if (CM_sub[i,j,tt] >= th_upA){
             Gr_subC = 0
           }
           
 ### Salt TOLERANT NON HALOPHYTE     
       if (CM_sub[i,j,tt] >= conc_highB){
         Gr_subB = ## reduced by some factor
          }
           if (CM_sub[i,j,tt] >= th_upB){
             Gr_subB[i,j,tt] = 0
           }
     
 ### Salt SENSITIVE NON HALOPHYTE  
     if (CM_sub[i,j,tt] >= conc_lowC){
       Gr_subC = ## reduced by some factor
        }
 
           if (CM_sub[i,j,tt] >= conc_highC){
             Gr_subC = ## even further reduced by some factor
             }
   
                 if (CM_sub[i,j,tt] >= th_upC){
                   Gr_subC = 0
                 }
 

# Waterlogging conditions
# time scale of days!
# Check waterlogging conditions, arbirtrary set to 90& of field capacity for 3 days in a row (after 3 days oxygen is assumed to be consumed and nitrogen reduced)
   
   if((M[i,j,t]/(n*Zr)) >= 0.9*soilpar$s_fc && (M[i,j,t-1]/(n*Zr))>= 0.9*soilpar$s_fc && (M[i,j,t-2]/(n*Zr))>= 0.9*soilpar$s_fc)
   {
     ## Growth factor c decreases
   }
 
 
 ##***********GERMINATION*********************************************************************************************************************
 ##
 ## Since I am NOT modelling the path and fate of every single seed and just implement seed dispersal anisotropic with the laplacian operator,
 ## to find out whether it is germination, I just ask the condition whether P of the previous timestep was zero (for that species)
 ## and if yes, the salinity levels are checked against a given threshold for germination
 
 # to dos: define germination salinity threshold germ_th for every individual species
 
 if(P[i,j,t-1] == 0 && CM[i,j,t] >= germ_th){  
   P[i,j,t] <- 0
 }
   
   