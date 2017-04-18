## MASTER FILE

# This file sets 


### Sourcing the grid, taudem etc, generates flowdir (ang) and slope raster
source("Rainfall.R")

### RUN TIME (how many days)
time <- 100
### Raster size
ext <- 200 ## EXTEND of PLOT in [m]
Z <- 3000 ##Groundwater depth in mm from 0 elevation
rows <- 100## rows of cells
cols <- 100 ## columns of cells

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



# Sourcing the cpp functions that define the constants for soil, veg and salt
sourceCpp("soilfun.cpp")
sourceCpp("vegfun.cpp")
sourceCpp("saltfun.cpp")
# 

# # creating parameter lists
#soilpar1 <- Soil_cpp("S Clay Loam")
soilpar1 <- Soil_cpp("C Sand")
vegpar1 <-Veg_cpp("Fantasy Tree")
saltpar1 <- Salt_cpp("Both")  ## other options: "Rain", "Both"


result<- SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,
                               saltpar = saltpar1, dims = list(rows=rows,cols=cols,time=time),
                               alpha_i =1.0, cn=0.01, Mn=0.04, Rain=Rain, slope=slp_matrix,Zras=Zras_matrix, rn=rn_matrix)


result$fields[[1]][1:10,1:10,2]
#str(result$fields)

result$fields[[11]]

library(rasterVis)
qr<-brick(result$fields[[6]][1:10,1:10,20:40])
levelplot(qr,main="P [g/m^2] ",sub="day 20 to day 40") 

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
