##################################################################################################################
##################################################################################################################
##################################################################################################################
###***************RASTER SCRIPT
## this file defines the raster the model is based on, defines its size and resolution,defines elevation 

library(raster)
library(shapefiles)
library(rgdal)

setwd("H:/real salt/Logan/TauDEM/TauDEM5Exe/")  ############### need to tidy this up later
######################### TauDEM/Dinf functions need to be in same directory as data...


####################################################################################################################
####################################################################################################################
##************BASE RASTER**************************************************************************************
## Generate a Raster object as a base raster, defines its size and resolution, number of cells
raster<- raster(ncol=10, nrow=10, xmn=0, xmx=40, ymn=0, ymx=40)

####################################################################################################################
####################################################################################################################
##************ELEVATION RASTER**************************************************************************************
elev <- raster ## elevation 
set.seed(100)
values(elev) <- runif(ncell(elev),0,1000)

####################################################################################################################
####################################################################################################################
##************DISTANCE TO GROUNDWATER RASTER************************************************************************


Z=3000 #### Groundwater depth in mm from 0 elevation

Zras<-raster(elev)
values(Zras)<-values(elev)+Z
plot(elev)
plot(Zras)

### Zras is the Raster with distance to groundwater table surface to the soil surface in mm


##************Preparation for TauDEM functions**********************************************************************
# needs to be .tif format for some reason
writeRaster(elev, filename="elev.tif", format="GTiff", overwrite=TRUE)
z=raster("elev.tif")
plot(z)

##********** Pitremove*********************************************************************
system("mpiexec -n 8 pitremove -z elev.tif -fel elevfel.tif")
fel=raster("elevfel.tif")
plot(fel)


####################################################################################################################
####################################################################################################################
##************ TAUDEM DInfINITY FLOW DIRECTIONS Raster**************************************************************
#

# DInf flow directions
system("mpiexec -n 8 DinfFlowdir -ang elevang.tif -slp elevslp.tif -fel elevfel.tif",show.output.on.console=F,invisible=F)
### Angles
ang=raster("elevang.tif")
plot(ang)
### Slope
slp=raster("elevslp.tif")
plot(slp)

###############################################################################################################
# the following code is from the TauDEM script from Tarboton's website, but wasnt used here (and there is heaps more)

# # D8 flow directions
# system("mpiexec -n 8 D8Flowdir -p rasterObjp.tif -sd8 rasterObjd8.tif -fel rasterObjfel.tif",show.output.on.console=F,invisible=F)
# p=raster("rasterObjp.tif")
# plot(p)
# sd8=raster("rasterObjd8.tif")
# plot(sd8)
# 
# # Contributing area
# system("mpiexec -n 8 AreaD8 -p rasterObjp.tif -ad8 rasterObjad8.tif")
# ad8=raster("rasterObjad8.tif")
# plot(log(ad8))
# #zoom(log(ad8))
# 
# # Grid Network 
# system("mpiexec -n 8 Gridnet -p rasterObjp.tif -gord rasterObjgord.tif -plen rasterObjplen.tif -tlen rasterObjtlen.tif")
# gord=raster("rasterObjgord.tif")
# plot(gord)
# #zoom(gord)
