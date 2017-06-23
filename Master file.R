## MASTER FILE

# This file sets 
library("devtools")
find_rtools()
library("Rcpp")
library("RcppArmadillo")

source("Rainfall.R")

time <- 300

# Rainlist <- list()
# # ## RAINFALL GENERATION
# for (k in 1:length(alpha)) {
#   
#   for (l in 1:length(lambda)) {
#     # generate the rainfall
#     # Rain <- Precip(time=1:5,alpha[k],lambda[l],delta)
#     Rainlist[[(k-1)+l]] <- data.frame(Precip(time,alpha[k],lambda[l],delta))
#   }}
# #Rain <- c(Rain)
# Rain <- Rainlist[[1]][,1]
# alpha <- seq(0.6,1.5,by=0.1)
# lambda <- seq(0.1,1,by=0.1)
# # alpha <- c(0.6,1.5) # 
# # lambda <- c(0.1,1) #
delta <- 0
alpha <- 0.91
lambda <- 0.16
Rain <- Precip(time,alpha,lambda,delta)
Rain <- c(Rain)

sourceCpp("soilfun.cpp")
sourceCpp("vegfun_static.cpp")
sourceCpp("saltfun.cpp") 
##sourceCpp("Model_withPlantInterference.cpp")  

# # creating parameter lists
#soilpar1 <- Soil_cpp("S Clay Loam")
soilpar1 <- Soil_cpp("S Clay Loam")
vegpar1 <-Veg_cpp()
saltpar1 <- Salt_cpp("Groundwater")  ## other options: "Rain", "Both", "None", "Groundwater"

### Raster size

Z <- 2000.0 ##Groundwater depth in mm from 0 elevation
rows <- 30## rows of cells
cols <- 30## columns of cells

# "diverseInput"
deltat <- 12 # temporal discretization, subdaily timesteps
gslp <- 0.05 # hillslope [%]
# infiltration
alpha_i <- 1.0
#Kinematic wave (runoff) paramters
cn <-0.01
Mn <- 0.04 # Manning's n
# EXTEND of PLOT in [m]
ext <- 200.0 
# soil moisture diffusivity
Dm <- 0.27 #Saco and Moreno-de las Heras,2013
# coefficient impact of plant interference
zeta <- 1.0
# ecosystem carrying capacity
P0 <-  500.0
# sensitivity to salinity
sigmaP <- 2.0
# interference parameters, competition and facilitation
b1 <- 0.9
b2 <- 0.1
q1 <-  1.0
q2 <-  2.0
# seed dispersal and diffusion
c1 <- 2.25 # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002   #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.0003 #Saco and Moreno-de las Heras,2013
Zr <- 400.0 # mm, Grass

#sourcing and running the model
sourceCpp("Model_revised.cpp")  
# sourceCpp("ModelTaudem.cpp")  

# visualization
library(rasterVis)
coul = brewer.pal(8, "YlGn")
coul = colorRampPalette(coul)(100)
qr<-brick(results$fields[[6]][4:26,4:26,150:200])
# results$fields[[10]][5:16,5:16,173]
levelplot(qr,main="P [g/m^2] ",sub="day 150 to day 200, salt from gw, alpha and lambda like in Tennant Creek, NT", col.regions = coul) #col.regions = YlGn.colors(20))
# animate(qr, n=1)
# results$fields[[6]][2:19,2:19,1:100]