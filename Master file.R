## MASTER FILE

# This file sets 


### Sourcing the grid, taudem etc, generates flowdir (ang) and slope raster
source("Rainfall.R")

### RUN TIME (how many days)
time <- 100
### Raster size
ext <- 200 ## EXTEND of PLOT in [m]
Z <- 3000.0 ##Groundwater depth in mm from 0 elevation
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
### preparations, put them somewhere else later
flowdir <- ang # ang is angle flowdir[i,j] from DInf TauDEM as raster
flowdir[is.na(flowdir)] <- 8.0   ###********************The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...
flowdir <- as.matrix(flowdir,nrow= nrow(flowdir), ncol=ncol(flowdir))
# slope raster generated from Rasters.R transformed into matrix
slp_matrix<- as.matrix(slp,nrow= nrow(slp), ncol=ncol(slp))
Zras_matrix <- as.matrix(Zras,nrow= nrow(Zras), ncol=ncol(Zras))

write.table(flowdir, "flowdir.txt")
write.table(slp_matrix, "slp_matrix.txt")
write.table(Zras_matrix, "Zras_matrix.txt")


# Sourcing the cpp functions that define the constants for soil, veg and salt
sourceCpp("soilfun.cpp")
sourceCpp("vegfun.cpp")
sourceCpp("saltfun.cpp")
# 

# # creating parameter lists
#soilpar1 <- Soil_cpp("S Clay Loam")
soilpar1 <- Soil_cpp("S Clay Loam")
vegpar1 <-Veg_cpp("Fantasy Tree")
saltpar1 <- Salt_cpp("Groundwater")  ## other options: "Rain", "Both", "None", "Groundwater"

sourceCpp("Model_largechanges.cpp")
system.time(
result <-SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,
                                saltpar = saltpar1, dims = list(rows=rows,cols=cols,time=time),
                                alpha_i =1.0, cn=0.01, Mn=0.04, Rain=Rain, slope=slp_matrix,Zras=Zras_matrix, flowdir = flowdir))
qr<-brick(result$fields[[8]][2:9,2:9,50:100])
levelplot(qr,main="P [g/m^2] ",sub="day 20 to day 40") 



result$fields[[4]][1:10,1:10,80:90]
          
#system.time()
result$fields[[1]][1:10,1:10,90]
#str(result$fields)
flowdir
result$fields[[8]][2:9,2:9,150:200]
result$fields[[18]]
result$fields[[2]][2:9,2:9,170:200]

library(rasterVis)
qr<-brick(result$fields[[6]][2:9,2:9,1:300])
#levelplot(qr,main="P [g/m^2] ",sub="day 20 to day 40") 
animate(qr, n=1)






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
