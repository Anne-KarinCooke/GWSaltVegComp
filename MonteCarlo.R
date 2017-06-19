### Script for simulation


# necessities:

# packages for Rtools, Rcpp, RcppArma...
library("devtools")
find_rtools() # check wheter it finds Rtools, otherwise Rcpp won't run
library("Rcpp")
library("RcppArmadillo")

# packages for the patch and diversity analysis 
library("SDMTools")
library("vegan") # 

# package for tibbles, suppossably better than data.frames
library("tibble")
# for random field DEM generation (surface elevation variance and corr.length)
library("geoR")


# ****************************************************************************************************************************************

# duration of simulation
time <- 17520 # days, 48 years

## Generate rainfall data
delta <- 0
alpha <- 0.91 #cm/event
lambda <- 0.16 #d/event

Rain <- Precip(time,alpha,lambda,delta)
Rain1 <- c(Rain*10) # to be in mm

# Read in actual rainfall data

TennantAirport <- read.csv("M:/Master thesis/tennant creek airport/IDCJAC0009_015135_1800_Data.csv", header= T, sep = ",")
Rain2 <-TennantAirport$Rainfall.amount..millimetres.
#length(Rain)

# source the parameter list functions
sourceCpp("soilfun.cpp")
sourceCpp("vegfun_static.cpp")
# # creating parameter lists
soilpar1 <- Soil_cpp("S Clay Loam") # soil type -> choose one and run simulation with that
vegpar1 <-Veg_cpp() # only contains fixed veg parameters. the ones varied are in simInput

# sOURCE THE MODEL
#source("FINAL_NAME.cpp")

# PARAMETERS THAT STAY THE SAME
# size of the grid
rows <- 30## rows of cells
cols <- 30 ## columns of cells
# EXTEND of PLOT in [m]
ext <- 200.0 
# "diverseInput"
deltat <- 12 # temporal discretization, subdaily timesteps
# infiltration
alpha_i <- 1.0
#Kinematic wave (runoff) paramters
cn <-0.01
Mn <- 0.04 # Manning's n
# soil moisture diffusivity
Dm <- 0.27 #Saco and Moreno-de las Heras,2013
# coefficient impact of plant interference
zeta <- 0.2
# ecosystem carrying capacity
P0 <-  500.0
# seed dispersal and diffusion
c1 <- 2.25; # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002 ;  #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.3 #Saco and Moreno-de las Heras,2013



# ****************************************************************************************************************************************
## SIMULATION
## Monte Carlo simulation

# number of runs
runs = 10000

# read in the parameter ranges
MC_par <- read.csv("MC_Parameters.txt", header=F)
attach(MC_par)
par_MC <- as.data.frame(matrix(0,nrow=runs,ncol=nrow(MC_par), byrow=F))
colnames(par_MC) <- MC_par[,1]
set.seed(runs)
# generate random combination of the varied input 
for (i in 1:nrow(MC_par)) {
  par_MC[,i] <- runif(runs,MC_par[i,2],MC_par[i,3])
}

# declaring the input lists
dims = list(rows=rows,cols=cols,time=time, Z=Z, ext = ext)

fixedInput = list(deltat = deltat, Dm = Dm, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, 
                  c1 = c1, c02 = c02, Dp = Dp, zeta = zeta, f=f)

simInput = list(Z =Z,ConcConst = ConcConst, CMgw = CMgw,gslp = gslp, 
                dA = dA, k1A = k1A,  b1A = b1A, b2A = b2A, q1A = q1A, q2A = q2A, sigmaPA= sigmaPA,
                dB = dB, k1B = k1B,  b1B = b1B, b2B = b2B, q1B = q1B, q2B = q2B, sigmaPB= sigmaPB,
                dC = dC, k1C = k1C,  b1C = b1C, b2C = b2C, q1C = q1C, q2C = q2C, sigmaPC= sigmaPC,
                ZrA =ZrA, ZrB =ZrB, ZrC =ZrC,
                sigma2 = sigma2, range = range)


# Setting up the Store for the simulation results
Store <- data_frame(simInput,  
                    P= array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
                    P_A = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
                    P_B = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
                    P_C = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
                    SmM = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
                    Zmatrix, elev_data_new) 

for (j in 1:runs) {
  
  simInput$Z <- Store$Z[j] ##Groundwater depth in mm from 0 elevation
  simInput$Zr <- Store$Zr[j] ## root depth
  simInput$ConcConst <- Store$ConcConst[j] ## salt concentration
  simInput$CMgw <- Store$CMgw[j]## salt concentration groundwater 
  simInput$gslp <- Store$gslp[j] # hillslope [%]
  
  ## species A
  simInput$k1A <- Store$k1A[j]## half saturation constant plant water uptake
  simInput$b1A <- Store$b1A[j] ## plant interference param. facilitation
  simInput$b2A <- Store$b2A[j]## plant interference param. competition
  simInput$q1A <- Store$q1A[j]## plant interference param. range faiclitation
  simInput$q2A <- Store$q2A[j]## plant interference param. range competition
  simInput$sigmaPA <- Store$sigmaPA[j] # sensitivity to salinity
  simInput$dA <- Store$dA[j]## plant mortality
  
  ## species B
  simInput$k1B <- Store$k1B[j]## half saturation constant plant water uptake
  simInput$b1B <- Store$b1B[j] ## plant interference param. facilitation
  simInput$b2B <- Store$b2B[j]## plant interference param. competition
  simInput$q1B <- Store$q1B[j]## plant interference param. range faiclitation
  simInput$q2B <- Store$q2B[j]## plant interference param. range competition
  simInput$sigmaPB <- Store$sigmaPB[j] # sensitivity to salinity
  simInput$dB <- Store$dB[j]## plant mortality
  
  ## species C
  simInput$k1C <- Store$k1C[j]## half saturation constant plant water uptake
  simInput$b1C <- Store$b1C[j] ## plant interference param. facilitation
  simInput$b2C <- Store$b2C[j]## plant interference param. competition
  simInput$q1C <- Store$q1C[j]## plant interference param. range faiclitation
  simInput$q2C <- Store$q2C[j]## plant interference param. range competition
  simInput$sigmaPC <- Store$sigmaPC[j] # sensitivity to salinity
  simInput$dC <- Store$dC[j]## plant mortality
  
  ## soil surface elevation heterogeneity parameters
  simInput$sigma2 <- Store$sigma2[j] # sigma squared, variance of surface elevation heterogeneity, normally distr
  simInput$range <- Store$range[j] # correlation length of surface elevation heterogeneity, normally distr
  
  # DEM generation
  A <- c(simInput$sigma2, simInput$range)
  # generating a Gaussian random field with geoR-package function grf()
  # grid 5 times finer than needed
  sim <-grf(((rows*cols)*5), grid = "reg", cov.pars = A)
  # Taking every 5th element
  elev_data <- as.matrix(sim$data)
  elev_data[,seq(nrows(elev_data), ncol(elev_data),5)]
  # adding the global slope to the elevation (function add_Slope_to_Elev() defined in model script)
  Store$elev_data_new <- add_Slope_to_Elev(elev_data, dims$rows, simInput$gslp, dims$ext)
  # recalculation of surface elevation variance
  sigma2_new <- var(elev_data_new) 
  # new variance is used in store (important for later lin.reg.)
  simInput$sigma2 <- sigma2_new
  # distance to groundwater table for every cell (function Z_matrix() defined in model script)
  Store$Zmatrix <-  Z_matrix(elev_data_new, rows, cols, Z)
  
  # calling the model 
  results <- SurfaceSoilSaltWBGRID(soilpar=soilpar1, vegpar=vegpar1,
                                   dims = dims, Rain=Rain, fixedInput = fixedInput, 
                                   simInput = simInput, elev = elev_data_new, Zras = Zmatrix)
  
  # storing the interesting results
  warmUp <- 100 ## days; start up phase of the model (maybe needs to be set higher)
  
  Store$P[j] <- results$fields[[9]][warmUp:time]
  Store$P_A[j] <- results$fields[[10]][warmUp:time] 
  Store$P_B[j] <- results$fields[[11]][warmUp:time]
  Store$P_C[j] <- results$fields[[12]][warmUp:time]
  Store$SmM[j] <- results$fields[[15]][warmUp:time] 
  
  
}

# ****************************************************************************************************************************************
## saving the results
write.table(Store, "Store.txt")

# End *************************************************************************************************************************************
# ****************************************************************************************************************************************

