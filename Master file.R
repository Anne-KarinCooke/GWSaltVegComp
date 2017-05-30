## MASTER FILE

# This file sets 


source("Rainfall.R")
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


sourceCpp("soilfun.cpp")
sourceCpp("vegfun.cpp")
sourceCpp("saltfun.cpp")
# 

# # creating parameter lists
#soilpar1 <- Soil_cpp("S Clay Loam")
soilpar1 <- Soil_cpp("S Clay Loam")
vegpar1 <-Veg_cpp("Fantasy Tree")
saltpar1 <- Salt_cpp("Groundwater")  ## other options: "Rain", "Both", "None", "Groundwater"

### RUN TIME (how many days)
time <- 5.0
### Raster size

Z <- 3000.0 ##Groundwater depth in mm from 0 elevation
rows <- 5.0## rows of cells
cols <- 5.0## columns of cells

# "diverseInput"
deltat <- 12 # temporal discretization, subdaily timesteps
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
P0 <-  200.0
# sensitivity to salinity
sigmaP <- 0.5
# interference parameters, competition and facilitation
 b1 <- 0.1
 b2 <- 0.9
 q1 <-  2.0
 q2 <-  4.0
# seed dispersal and diffusion
c1 <- 2.25; # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002 ;  #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.3 #Saco and Moreno-de las Heras,2013

sourceCpp("Model_withPlantInterference.cpp")  

SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,saltpar = saltpar1, 
                      dims = list(rows=rows,cols=cols,time=time, Z=Z),
                      diverseInput = list (deltat = deltat, gslp = gslp, ext=ext, Dm = Dm, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, sigmaP= sigmaP,
                      c1 = c1, c02 = c02, Dp = Dp, b1 = b1, b2 = b2, q1 = q1, q2 = q2, zeta = zeta), 
                      Rain=Rain) 


### results Visualization


# library(rasterVis)
# qr<-brick(result$fields[[17]][2:9,2:9,10:110])
# levelplot(qr,main="P [g/m^2] ",sub="day 1 to day 100, salt from gw, randomly varied alpha and lambda",col.regions = terrain.colors(20)) 
# animate(qr, n=1)
# qr<-brick(result$fields[[8]][2:9,2:9,50:100])
# levelplot(qr,main="P [g/m^2] ",sub="day 20 to day 40") 
# result$fields[[11]][1:10,1:10,90:100]
# result$fields[[6]][1:10,1:10,300:400]


