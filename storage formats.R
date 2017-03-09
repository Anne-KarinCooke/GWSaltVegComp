# 
# grids <- ncell(raster)  ### dont forget to define raster


# Storage vectors for the daily steps 

M <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain))  # soil moisture [mm]

####################################
###Difference! h is now defined as a raster object; # flow depth [mm]
h <- brick(raster,nl=length(Rain))
values(h)<-0
###################################  
P <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) #biomass density []
CM<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # Salt concentration in soil water in g/L or g/mm
SmI<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # Salt mass in infiltrating water [g]
SmM <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # Salt mass in soil water [g]
In <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # infiltration [mm]
Svir <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # virtual saturation
flux<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain))  # vertical water flux, capillary rise and drainage


q<-brick(raster,nl=length(Rain))
values(q)<-0
###RUNON
runon<-brick(raster,nl=length(Rain))
values(runon)<-0


Diff<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain)) # divergence of soil moisture from one grid cells to the next


# We decided to split the numerical calculations for the daily into 12 substeps.
deltat <- 12 # split in 12 increments
#mass balance
mb <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(Rain))


# Storage vectors for the substeps are initialized.

M_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat)) # soil moisture

#########################################################

###Difference! h_sub is now defined as a raster object; # flow depth [mm] 
h_sub<- brick(raster,nl=length(deltat)) # flow depth in [mm]
values(h_sub)<-0

########################################################



I_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  #Infiltration
WU_sub <-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # Water uptake in mm
P_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # Biomass density
Gr_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # Growth of biomass
Mo_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # Mortality of biomass
SmI_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  #salt mass in infiltration water
SmM_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # salt mass in soil
CM_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat)) # salt concentration in soil 
Svir_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # virtual saturation
flux_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))  # drainage/capillary rise


##q_sub<-matrix(0,nrow= deltat, ncol =grids)  # overland flow  
q_sub<-brick(raster,nl=length(deltat))
values(q_sub)<-0

### RUNON
runon_sub<-brick(raster,nl=length(deltat))
values(runon_sub)<-0


Diff_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat)) # divergence of soil moisture from one grid cells to the next

###h.old<-rep(0,grids)
h.old<- raster(raster) # flow depth in [mm]
values(h.old)<-0

######################
P.old<- raster(raster) 
values(P.old)<-0
M.old<- raster(raster) 
values(M.old)<-0
SmI.old<- raster(raster) 
values(SmI.old)<-0
CM.old<- raster(raster) 
values(CM.old)<-0
Svir.old<- raster(raster) 
values(Svir.old)<-0



L_salt<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))
U_salt<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))

mb_sub<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=length(deltat))

# #initital values to start the simulation.
# M[1,] <- 5
# h[1,] <- 10 
# P[1,] <- 10
# CM[1,]<- 0
# Svir[1,] <- s_fc



