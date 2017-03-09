# R= 0.0821 # gas constant [ L atm mol-1 K-1 ]
# Temp= 298 # temperature in Kelvin,  25 degrees Celsius
# Mm= 58.44 # molar mass of NaCl in g/mol

#SOIL :
# Sandy Clay Loam
n<-0.367 # porosity
# more soil variables for evaporation & losses
# Hydraulic conductivity
K_s<-52.08*10 # mm/day
# campbell's b
b<-6.4069 # neurotheta sandy clay loam
# van Genuchten parameters
#     avg <- 0.0521
#     nvg <- 1.237
s_fc<-0.2677/n # Field capacity
# This is the bubbling pressure
psi_s_bar<--1.2E-3 #

h1bar =  -psi_s_bar 
hb = psi_s_bar*-10^5 # mm


## OVERLAND FLOW parameter
# cn= 1 #conversion factor cn
# Mn= 10 #Manning's n, Mn
# Sl = 0.09 # Slope Sl



# #................................................
# Vegetation 1 (Grass)
# paspalum secateum F-I and R-I, 2004

Zr = 400 # soil depth (mm)   Check Also Table 2...Fernandez-Illescas and Rodriguez-Iturbe...2001


# parameters describing the root zone   
vegpar <- list(Zr = Zr)


# parameters describing the soil
soilpar <- list(b = b, n = n, s_fc = s_fc, K_s = K_s, 
                psi_s_bar = psi_s_bar, h1bar = h1bar, hb = hb, cn=cn, Mn=Mn, Sl=Sl)
# parameters describing plant dynamics and salt features

alpha_i=1 #maximum infiltration rate per day, This needs to be a fraction of h (p117 Saco and Moreno-Las Heras) 

k=12 # Saco et al, 2013
W0=0.2 # Saco et al, 2013
gmax=0.05 # Saco et al, 2013
k1=5 # Saco et al, 2013
c=10  # Saco et al, 2013
f= 0.8  # f is the soil salt leaching efficiency (whether some salt is retained)
ConcConst = 0.1 # ConcConst is the concentration of the salt in the infiltrating water in g/l
CM.gw = 0.1 # salt concentration in groundwater
d=0.24 # fraction of plant mortality

par <- list(alpha_i=alpha_i,k=k, W0=W0, gmax=gmax, k1=k1, c=c, f=f, ConcConst=ConcConst, CM.gw= CM.gw, d=d)




# library(raster)
# ## Generate a RasterLayer object
# 
# rasterObj <- raster(ncol=4, nrow=4, xmn=0, xmx=40, ymn=0, ymx=40)
# # set.seed(0)
# values(rasterObj) <- runif(ncell(rasterObj))
# 
# #values have to be a list or matrix of data....
# 
# # plot(rasterObj)
# 
# # res(rasterObj)
# 
# 
# 
# plot(rasterObj)
# 
# ###### GRID CELLS


# grids <- ncell(rasterObj)
# grids <- 4
# grid <-c(1:grids)
# 
# ### DIVERGENCE, diffusion
# Dm is the soil moisture diffusivity parameter
# Dm = 0.27 # from Saco et al, 2013, in m^2/d
# #distance between the "buckets" in mm
# dist= 5000

#mm GROUNDWATER DEPTH from the surface
# Z<-c(grid)
# Z[1]<- 3000 # mm 
# for (i in 1:length(grid)){
#   Z[i+1] <- Z[i]-dist*sin(soilpar$Sl)
# }
#  Z[2]<- Z[1]-dist*sin(soilpar$Sl)
# Z[3]<- Z[2]-dist*sin(soilpar$Sl) # downhill the distance to the groundwater becomes smaller
# Z
