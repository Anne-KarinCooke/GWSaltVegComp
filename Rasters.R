##################################################################################################################
##################################################################################################################
##################################################################################################################
###***************RASTER SCRIPT
## this file defines the raster the model is based on, defines its size and resolution,defines elevation 
#setwd("H:/real salt/GWSaltVegComp/thesis project/GWSaltVegComp")
library(raster)
library(shapefiles)
library(rgdal)

#  setwd("H:/real salt/Logan/TauDEM/TauDEM5Exe/")  ############### need to tidy this up later
######################### TauDEM/Dinf functions need to be in same directory as data...


####################################################################################################################
####################################################################################################################
##************BASE RASTER**************************************************************************************
## Generate a Raster object as a base raster, defines its size and resolution, number of cells

## EXTEND of PLOT in [m]
raster<- raster(ncol=cols, nrow=rows, xmn=0, xmx=ext, ymn=0, ymx=ext)


####################################################################################################################
####################################################################################################################
##************ELEVATION RASTER**************************************************************************************
elev <- raster ## elevation 
set.seed(100)
microdepth<-0.5 #[m]
values(elev) <- runif(ncell(elev),0,microdepth)


####################################################################################################################
####################################################################################################################
##************GLOBAL SLOPE************************************************************************
###  
gslp=0.02 
elev[1,]<-elev[1,]+(gslp*ext)
for (i in 2:nrow(elev)){
  elev[i,]<-elev[i,]+(gslp*(ext-((ext/nrow(elev))*(i-1))))
}
#(elev)

# terrslp <- terrain(fel, oprt="slope")
# plot(terrslp)
# values(terrslp)
# values(slp)
# plot(slp)
# terrain


## recalculation of 

####################################################################################################################
####################################################################################################################
##************DISTANCE TO GROUNDWATER RASTER************************************************************************


#Z=3000 #### Groundwater depth in mm from 0 elevation

Zras<-raster(elev)
values(Zras)<-(values(elev)*1000)+Z

# plot(elev)
# plot(Zras)

### Zras is the Raster with distance to groundwater table surface to the soil surface in mm


##************Preparation for TauDEM functions**********************************************************************
# needs to be .tif format for some reason
writeRaster(elev, filename="elev.tif", format="GTiff", overwrite=TRUE)
z=raster("elev.tif")
 # plot(z)

##********** Pitremove*********************************************************************
system("mpiexec -n 8 pitremove -z elev.tif -fel elevfel.tif")
fel=raster("elevfel.tif")
#plot(fel, main="pits removed")



####################################################################################################################
####################################################################################################################
##************ TAUDEM DInfINITY FLOW DIRECTIONS Raster**************************************************************
#

# DInf flow directions
system("mpiexec -n 8 DinfFlowdir -ang elevang.tif -slp elevslp.tif -fel elevfel.tif",show.output.on.console=F,invisible=F)
### Angles
ang=raster("elevang.tif")
# plot(ang)

### Slope
slp=raster("elevslp.tif")
slp[is.na(slp)] <- 0


# 
# plot(elev)
# path<-flowPath(elev,elev[1,1])
# plot(path)

# fd<-terrain(elev, opt="flowdir")
# path <- flowPath(fd, 20)
# xy <- xyFromCell(fd, path)
# plot(elev)
# lines(xy)

# plot(slp)

# Dinf contributing area
# system("mpiexec -n 8 AreaDinf -ang elevang.tif -sca elevsca.tif")
# sca=raster("elevsca.tif")
# plot(sca)
# plot(log(sca))
# zoom(log(sca))
