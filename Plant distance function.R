## Plant interaction and salt response code***************************************************
#*********************************************************************************************

# 1. Old plant functions as a reminder
# 2. Definition of storage arrays for plant variables for every single species
# 3. Water uptake in that cell, as a sum of water uptake of all present species
# 4. Setting of specific threshold concentration constants
# 5. Plant salt responses (influence on growth) 
# 6. Impact of Waterlogging
# 7. Salt impact on Germiantion
# 8. Salt impact on long-term seed production
# 9. Interaction functions (distance betweeen plants)


# 1. *****************************************************************************************
# # Plant water uptake
# WU <- function(M,P,par){ 
#   # using Svir in here means scaling Svir back to M, easier to do at Svir in balances  
#   #  WU=par$gmax*((M*(1+Svir))/(((M*(1+Svir))+par$k1)))*P 
#   WU=par$gmax*(M/((M+par$k1)))*P 
#   
#   return(WU)
# }
# 
# #Plant Growth
# 
# Gr <- function(M,P,par) { 
#   
#   Gr = par$c*WU(M,P,par)
#   
#   return(Gr)
# }
# 
# 
# ## Plant mortality function WITH SALT INFLUENCE, WITH VIRTUAL SATURATION
# 
# Mo <- function(P,M,Svir,par) {
#   # needs to be M/Svir because both are "large" numbers
#   # you want a number ~1 for multiplication, or <0.1 for addition
#   Mo = P*(par$d*(M/Svir))
#   
#   return(Mo)
#   
# }


# 2. ***************************************************************************************** 
#  Definition of storage arrays for plant variables for every single species
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
 
 # 3. ******************************************************************************************************
 # Water uptake in that cell, as a sum of water uptake of all present species
 # Water uptake, sum of WU_subA,B,C
 
 WU_sub[i,j,tt] <- WU_subA[i,j,tt] + WU_subB[i,j,tt] + WU_subC[i,j,tt]
 

 
 # 4. ******************************************************************************************************
 #Setting of specific threshold concentration constants
 ## Response of 3 diff functional plant types (Maas and Hoffman, 1977)
 ## (neither temporal threshold nor water logging considered yet )

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
   # th_low 
   # # low
   # conc_low -> irrelevant for this functional type
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

 # 5. *************************************************************************************************  
 # Plant salt responses (influence on growth)    
   
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
 
# 6. ***************************************************************************************************
# Waterlogging conditions
# time scale of days!
# Check waterlogging conditions, arbirtrary set to 90& of field capacity for 3 days in a row (after 3 days oxygen is assumed to be consumed and nitrogen reduced)
   
   if((M[i,j,t]/(n*Zr)) >= 0.9*soilpar$s_fc && (M[i,j,t-1]/(n*Zr))>= 0.9*soilpar$s_fc && (M[i,j,t-2]/(n*Zr))>= 0.9*soilpar$s_fc)
   {
     ## Growth factor c decreases
   }

 # 7. ********************************************************************************************************************************
 ## Germination
 ## Since I am NOT modelling the path and fate of every single seed and just implement seed dispersal anisotropic with the laplacian operator,
 ## to find out whether it is germination, I just ask the condition whether P of the previous timestep was zero (for that species)
 ## and if yes, the salinity levels are checked against a given threshold for germination
 
 # to dos: define germination salinity threshold germ_th for every individual species
 
 if(P[i,j,t-1] == 0 && CM[i,j,t] >= germ_th){  
   P[i,j,t] <- 0
 }
   
# 8. ****************************************************************************************************************
# Long-term fertility reduction (reduced seed production due to salinity)   

   # if (CM[i,j,t] > 3 months... then seed dispersal is reduced...
   
   CM[i,j,t+x] - CM[i,j,t]
   
   filter(CM[i,j,], CM[i,j,t]  )
   filterWindow <- seq(1,90)
   filter(x, filterWindow, method = "convolution", sides=1)
   
   # filter(x, filter, method = c("convolution", "recursive"),
   #        sides = 2, circular = FALSE, init)

         
# 9. *************************************************************************************************************   
#  Interaction functions (distance betweeen plants)   
   
library(raster)  
library(igraph)
library(gdistance)

 # Lefever, R., Lejeune, O., & Couteron, P. (2001). Generic modelling of vegetation patterns. A case study of tiger bush in sub-saharian sahel. In Mathematical models for biological pattern formation (pp. 83-112). Springer New York.
 # Eq. 7 on page 89
 # Weight function 
 # L = parameter that controls steepness with which weighting function varies
 # r is the distance 
 
 wfunc <- function(r, L){
   w <- (1/(2*pi*(L*L)))*exp(-(abs(r)*abs(r))/(2*(L*L)))
   return(w)
 }
 # This is basically a kernel function! The R package "spatialfil" provides such a function, I actually didnt need to write it myself
 #convKernel(sigma = 1.4, k = c("gaussian", "LoG", "sharpen", "laplacian",
 #                             "emboss", "sobel"))
 #install.packages("spatialfil")
 # library(spatialfil)
 # convKernel(sigma=1.4,k="gaussian")

 # Distance of a cell to all other cells, considering 16 directions 
 # http://personal.colby.edu/personal/m/mgimond/Spatial/Distance_rook_vs_queen_vs_knight.html
 #Distance <-accCost(geoCorrection(transition(Pras, transitionFunction=function(x){1},16,symm=FALSE), scl=FALSE),A)
 # insert some kind of raster and a point A from which these distances shall be calculated
 # plot(Distance)
 
 if (P_subA[i,j,tt] > 0) 
   # Raster with cells that 
   
   r   <- raster(nrows=100,ncols=100,xmn=0,ymn=0,xmx=100,ymx=100)
 Pras <- r # raster(P[,,tt])
 values(Pras) <- 1

 # Distance <- brick(r, nl=ncell(Pras))
 # Distance <- brick(r, nl=i*j)
   for (i in nrow(Pras)){
     for (j in ncol(Pras)){

       
      a<- as.array(values(Distance <-accCost(geoCorrection(transition(Pras, transitionFunction=function(x){1},16,symm=FALSE),scl=FALSE),c(i,j))), dim=c(i,j))
       
       # interference <- wfunc(Distance[i,j],0.3)*Pras
     #   P_sub[i,j,tt+1] <- P.old[i,j] + Gr_sub[i,j,tt]- Mo_sub[i,j,tt] + interference
      }
   }
 a
 
 r   <- raster(nrows=100,ncols=100,xmn=0,ymn=0,xmx=100,ymx=100)
 Pras <- r # raster(P[,,tt])
 values(Pras) <- 1
 
 # Distance <- array(dim=c(nrow(Pras),ncol(Pras),(nrow(Pras)*ncol(Pras))))
 # 
 # DistanceBrick <- brick(nrow=nrow(Pras),ncol=ncol(Pras), ncell=(nrow(Pras)*ncol(Pras)), nlayers=(nrow(Pras)*ncol(Pras)))
 # values(DistanceBrick) <-1
 # DistanceBrick[,,,,1]
 # 
 # for (i in nrow(Pras)){
 #   for (j in ncol(Pras)){
 #     
 #     a<-brick(accCost(geoCorrection(transition(Pras, transitionFunction=function(x){1},16,symm=FALSE),scl=FALSE),c(i,j)),nlayers=(i*j))
 #     
 #   }
 # }
 # a
 

 r   <- raster(nrows=5,ncols=5,xmn=0,ymn=0,xmx=100,ymx=100)
 Pras <- r # raster(P[,,tt])
 values(Pras) <- 1
 Prasmatr<-as.matrix(Pras)
 
 DistArray <- array(dim=c(nrow(Prasmatr),ncol(Prasmatr),(nrow(Prasmatr)*ncol(Prasmatr))))
 

 for (i in 1:nrow(Prasmatr)){
    for (j in 1:ncol(Prasmatr)){
 # This creates (first raster then matrix) of distances from the cell c(i,j) to all other cells in the grid
  A<-c(i,j)
 Distance<-accCost(geoCorrection(transition(Pras, transitionFunction=function(x){1},16,symm=FALSE),scl=FALSE),A)
 DistMatr<- as.matrix(Distance)
 DistMatrRot <-t(DistMatr)[,ncol(DistMatr):1] 

  DistArray[,,i*j]<- DistMatrRot
 
    }
   }
 
 DistArray[,,]
 
 apply(DistArray,)
 
 # Distance[3,6,3*6]
 # Sys.time(d)
 #   d<-function(){
 # Distance <-accCost(geoCorrection(transition(Pras, transitionFunction=function(x){1},16,symm=FALSE),scl=TRUE),c(Pras[,]))
 # }
 # acc
 # plot(Distance)
 # plot(interference)
 #   

 
 ##
 # Plant dispersal
 # wind, animals...
 Dp <- 0.3 #m^2d^-1
 
 

 
 #Recycling bin
 
 # g <- make_ring(5)
 # laplacian_matrix(g)
 # laplacian_matrix(g, norm=TRUE)
 # laplacian_matrix(g, norm=TRUE, sparse=FALSE)
 
 # r   <- raster(nrows=100,ncols=100,xmn=0,ymn=0,xmx=100,ymx=100)
 # r[] <- rep(1, ncell(r))
 # plot(r)
 # h16  <- transition(r, transitionFunction=function(x){1},16,symm=FALSE)
 # h16   <- geoCorrection(h16, scl=FALSE)
 # A <- c(20,50)
 # h16.acc <- accCost(h16,A)
 # plot(h16.acc)
 
 # r <- raster(ncol=10,nrow=10)
 # r[] <- 1
 # r[48] <- 2
 # r[66:68] <- 3
 # d <- gridDistance(r,origin=2,omit=3) 
 # plot(d)
 # x <- c(0, 1, 1.1, 1.5, 1.9, 2.8, 2.9, 3.5)
 # n <- length(x)
 # 
 # xgrid <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01)
 # 
 # h <- 0.4
 # bumps <- sapply(x, function(a) gauss((xgrid - a)/h)/(n * h))
 # 
 # 
 # 
 # 
 