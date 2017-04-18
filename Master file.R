## MASTER FILE

# This file sets 


### Sourcing the grid, taudem etc, generates flowdir (ang) and slope raster
source("Rasters.R")

### RUN TIME (how many days)
time <- 100
### Raster size
ext <- 40 ## EXTEND of PLOT in [m]
Z <- 3000 ##Groundwater depth in mm from 0 elevation
rows <- 10## rows of cells
cols <- 10 ## columns of cells

# alpha <- c(0.6,1.5) # 
alpha <- seq(0.6,1.5,by=0.1)
lambda <- seq(0.1,1,by=0.1)
# lambda <- c(0.1,1) #  
delta <- 0

# ## RAINFALL GENERATION
for (k in 1:length(alpha)) {
  
  for (l in 1:length(lambda)) {
    # generate the rainfall
    Rain <- Precip(time,alpha[k],lambda[l],delta)
    Rainlist <- data.frame(Precip(time,alpha[k],lambda[l],delta))
  }}

Rain <- c(Rain)


### Sourcing the grid, taudem etc, generates flowdir (ang) and slope raster
source("Rasters.R")
### Source the runon raster (flowdir)
source("flowdir.R")
#Runon raster generated from flowdir.R
rn_matrix<- as.matrix(rn,nrow= nrow(rn), ncol=ncol(rn))
# slope raster generated from Rasters.R
slp_matrix<- as.matrix(slp,nrow= nrow(slp), ncol=ncol(slp))
Zras_matrix <- as.matrix(Zras,nrow= nrow(Zras), ncol=ncol(Zras))


  
sourceCpp("Model_heteroRain_heteroTopo.cpp")

soilpar_in <- soil_simple()
vegpar_in <- veg_simple()
saltpar_in <- salt_simple()
result<- SurfaceSoilSaltWBGRID(soilpar=soil_simple(), vegpar=veg_simple(),
                               saltpar = salt_simple(), dims = list(rows=rows,cols=cols,time=time),
                               alpha_i =1.0, cn=0.01, Mn=0.04, Rain=Rain, slope=slp_matrix,Zras=Zras_matrix, rn=rn_matrix)


result$fields[[1]][1:10,1:10,20]
#str(result$fields)
result$fields[[6]]





# #   
#   result<- SurfaceSoilSaltWBGRID(soilpar=soil_simple(), vegpar=veg_simple(),
#                                  saltpar = salt_simple(), dims = list(rows=10,cols=10,time=100),
#                                  alpha_i =1.0, cn=0.01, Mn=0.04, Rain=Rain, slope=0.001,Zras=1000.0)
#   result$fields[[1]][1:10,1:10,2]
# #str(result$fields)
#   result$fields[[6]]


# SurfaceSoilSaltWBGRID(alpha_i =1.0, cn=0.01, Mn=0.04, Rain=1.0, slope=0.001,Zras=1000.0, soilpar=soilpar_in, vegpar=vegpar_in,saltpar=saltpar_in) 

# hdata <- as.data.frame(result$fields[1])
# write.table(hdata, "C:/Users/acoo7451/Desktop/hdata.txt", sep="\t")

# Grid_run()
