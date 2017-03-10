# 
# grids <- ncell(raster)  ### dont forget to define raster


# Storage vectors for the daily steps 

M <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(Rain)))  # soil moisture [mm]

####################################
###
h<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain))))  # vertical water flux, capillary rise and dlength(length(Rain))age



###################################  
P <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) #biomass density []
CM<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # Salt concentration in soil water in g/L or g/mm
SmI<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # Salt mass in infiltrating water [g]
SmM <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # Salt mass in soil water [g]
In <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # infiltration [mm]
Svir <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # virtual saturation
flux<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain))))  # vertical water flux, capillary rise and dlength(length(Rain))age

q<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain))))  # vertical water flux, capillary rise and dlength(length(Rain))age

# q<-brick(raster,nl=length(Rain))
# 
# values(q)<-0
###RUNON

runon<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain))))


Diff<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),length(length(Rain)))) # divergence of soil moisture from one grid cells to the next


# We decided to split the numerical calculations for the daily into 12 substeps.
deltat <- 12 # split in 12 increments
#mass balance
mb <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
            


# Storage vectors for the substeps are initialized.

M_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) # soil moisture


#########################################################

###Difference! h_sub is now defined as a raster object; # flow depth [mm] 
# h_sub<- brick(raster,nl=deltat) # flow depth in [mm]
# values(h_sub)<-0
h_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) # h


########################################################



I_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  #Infiltration
WU_sub <-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # Water uptake in mm
P_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # Biomass density
Gr_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # Growth of biomass
Mo_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # Mortality of biomass
SmI_sub <- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  #salt mass in infiltration water
SmM_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # salt mass in soil
CM_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) # salt concentration in soil 
Svir_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # virtual saturation
flux_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))  # dlength(Rain)age/capillary rise


##q_sub<-matrix(0,nrow= deltat, ncol =grids)  # overland flow  
# q_sub<-brick(raster,nl=deltat)
# values(q_sub)<-0
q_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))


### RUNON

runon_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))



Diff_sub<- array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat)) # divergence of soil moisture from one grid cells to the next

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



L_salt<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))
U_salt<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))

mb_sub<-array(matrix(0,nrow= nrow(raster), ncol =ncol(raster)),dim=c(nrow(raster),ncol(raster),deltat))




