## MASTER FILE

# This file sets 
library("devtools")
find_rtools()
library("Rcpp")
library("RcppArmadillo")

source("Rainfall.R")
# alpha <- c(0.6,1.5) # 
alpha <- seq(0.6,1.5,by=0.1)
lambda <- seq(0.1,1,by=0.1)
# lambda <- c(0.1,1) #  
delta <- 0
time <- 8

# Rainlist <- list()
# # ## RAINFALL GENERATION
# for (k in 1:length(alpha)) {
#   
#   for (l in 1:length(lambda)) {
#     # generate the rainfall
#     # Rain <- Precip(time=1:5,alpha[k],lambda[l],delta)
#     Rainlist[[(k-1)+l]] <- data.frame(Precip(time,alpha[k],lambda[l],delta))
#   }}
# 
# #Rain <- c(Rain)
# Rain <- Rainlist[[1]][,1]
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

Rain



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
rows <- 20## rows of cells
cols <- 20 ## columns of cells

# "diverseInput"
deltat <- 5 # temporal discretization, subdaily timesteps
gslp <- 0.02 # hillslope [%]
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
zeta <- 0.2
# ecosystem carrying capacity
P0 <-  500.0
# sensitivity to salinity
sigmaP <- 0.05
# interference parameters, competition and facilitation
b1 <- 0.9
b2 <- 0.1
q1 <-  0.1
q2 <-  0.3
# seed dispersal and diffusion
c1 <- 2.25; # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002 ;  #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.3 #Saco and Moreno-de las Heras,2013
Zr <- 400.0 # mm, Grass

sourceCpp("Model_revised.cpp")  
results <- SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,saltpar = saltpar1,
                                 dims = list(rows=rows,cols=cols,time=time, Z=Z),  Rain=Rain,
                                 diverseInput = list (deltat = deltat, gslp = gslp, ext=ext, Zr=Zr, Dm = Dm, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, sigmaP= sigmaP,
                                                      c1 = c1, c02 = c02, Dp = Dp, b1 = b1, b2 = b2, q1 = q1, q2 = q2, zeta = zeta))
library(rasterVis)
coul = brewer.pal(8, "YlGn")
coul = colorRampPalette(coul)(100)


qr<-brick(results$fields[[6]][2:19,2:19,10:40])

levelplot(qr,main="P [g/m^2] ",sub="day 10 to day 40, salt from gw, randomly varied alpha and lambda", col.regions = coul) #col.regions = YlGn.colors(20))
# animate(qr, n=1)
# results$fields[[6]][2:19,2:19,1:100]