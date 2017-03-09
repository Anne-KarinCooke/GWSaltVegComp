


library(raster)
## Generate a RasterLayer object

rasterObj <- raster(ncol=10, nrow=10, xmn=0, xmx=40, ymn=0, ymx=40)
# set.seed(0)
values(rasterObj) <- runif(ncell(rasterObj),0,10)

writeRaster(rasterObj, filename="rasterObj.tif", format="GTiff", overwrite=TRUE)


#  R script to run TauDEM

library(raster)
library(shapefiles)
library(rgdal)

# Set working directory to your location

z=raster("rasterObj.tif")

plot(z)

setwd("H:/real salt/Logan/TauDEM/TauDEM5Exe/")

# Pitremove
system("mpiexec -n 8 pitremove -z rasterObj.tif -fel rasterObjfel.tif")
fel=raster("rasterObjfel.tif")
plot(fel)


# # D8 flow directions
# system("mpiexec -n 8 D8Flowdir -p rasterObjp.tif -sd8 rasterObjd8.tif -fel rasterObjfel.tif",show.output.on.console=F,invisible=F)
# p=raster("rasterObjp.tif")
# plot(p)
# sd8=raster("rasterObjd8.tif")
# plot(sd8)
# 
# 
# 
# # Contributing area
# system("mpiexec -n 8 AreaD8 -p rasterObjp.tif -ad8 rasterObjad8.tif")
# ad8=raster("rasterObjad8.tif")
# plot(log(ad8))
# #zoom(log(ad8))
# 
# 
# # Grid Network 
# system("mpiexec -n 8 Gridnet -p rasterObjp.tif -gord rasterObjgord.tif -plen rasterObjplen.tif -tlen rasterObjtlen.tif")
# gord=raster("rasterObjgord.tif")
# plot(gord)
# #zoom(gord)

# DInf flow directions
system("mpiexec -n 8 DinfFlowdir -ang rasterObjang.tif -slp rasterObjslp.tif -fel rasterObjfel.tif",show.output.on.console=F,invisible=F)
ang=raster("rasterObjang.tif")
plot(ang)
values(ang)
slp=raster("rasterObjslp.tif")
plot(slp)


# 
# # Dinf contributing area
# system("mpiexec -n 8 AreaDinf -ang rasterObjang.tif -sca rasterObjsca.tif")
# sca=raster("rasterObjsca.tif")
# plot(log(sca))
# # zoom(log(sca))

# # Threshold
# system("mpiexec -n 8 Threshold -ssa rasterObjad8.tif -src rasterObjsrc.tif -thresh 100")
# src=raster("rasterObjsrc.tif")
# plot(src)
# # zoom(src)
# 
# # a quick R function to write a shapefile
# makeshape.r=function(sname="shape",n=1)
# {
#   xy=locator(n=n)
#   points(xy)
#   
#   #Point
#   dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
#   ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
#   ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
#   write.shapefile(ddShapefile, sname, arcgis=T)
# }
# 
# makeshape.r("ApproxOutlets")
# 
# # Move Outlets
# system("mpiexec -n 8 moveoutletstostreams -p rasterObjp.tif -src rasterObjsrc.tif -o approxoutlets.shp -om Outlet.shp")
# outpt=read.shp("outlet.shp")
# approxpt=read.shp("ApproxOutlets.shp")
# 
# plot(src)
# points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
# points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)
# 
# zoom(src)
# # 
# 
# # Contributing area upstream of outlet
# system("mpiexec -n 8 Aread8 -p rasterObjp.tif -o Outlet.shp -ad8 rasterObjssa.tif")
# ssa=raster("rasterObjssa.tif")
# plot(ssa) 
# 
# 
# # Threshold
# system("mpiexec -n 8 threshold -ssa rasterObjssa.tif -src rasterObjsrc1.tif -thresh 2000")
# src1=raster("rasterObjsrc1.tif")
# plot(src1)
# zoom(src1)
# 
# # Stream Reach and Watershed
# system("mpiexec -n 8 Streamnet -fel rasterObjfel.tif -p rasterObjp.tif -ad8 rasterObjad8.tif -src rasterObjsrc1.tif -o outlet.shp -ord rasterObjord.tif -tree rasterObjtree.txt -coord rasterObjcoord.txt -net rasterObjnet.shp -w rasterObjw.tif")
# plot(raster("rasterObjord.tif"))
# zoom(raster("rasterObjord.tif"))
# plot(raster("rasterObjw.tif"))
# 
# # Plot streams using stream order as width
# snet=read.shapefile("rasterObjnet")
# ns=length(snet$shp$shp)
# for(i in 1:ns)
# {
#   lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
# }
# 
# # Peuker Douglas stream definition
# system("mpiexec -n 8 PeukerDouglas -fel rasterObjfel.tif -ss rasterObjss.tif")
# ss=raster("rasterObjss.tif")
# plot(ss)
# zoom(ss)
# 
# #  Accumulating candidate stream source cells
# system("mpiexec -n 8 Aread8 -p rasterObjp.tif -o outlet.shp -ad8 rasterObjssa.tif -wg rasterObjss.tif")
# ssa=raster("rasterObjssa.tif")
# plot(ssa)
# 
# #  Drop Analysis
# system("mpiexec -n 8 Dropanalysis -p rasterObjp.tif -fel rasterObjfel.tif -ad8 rasterObjad8.tif -ssa rasterObjssa.tif -drp rasterObjdrp.txt -o outlet.shp -par 5 500 10 0")
# 
# # Deduce that the optimal threshold is 300 
# # Stream raster by threshold
# system("mpiexec -n 8 Threshold -ssa rasterObjssa.tif -src rasterObjsrc2.tif -thresh 300")
# plot(raster("rasterObjsrc2.tif"))
# 
# # Stream network
# system("mpiexec -n 8 Streamnet -fel rasterObjfel.tif -p rasterObjp.tif -ad8 rasterObjad8.tif -src rasterObjsrc2.tif -ord rasterObjord2.tif -tree rasterObjtree2.dat -coord rasterObjcoord2.dat -net rasterObjnet2.shp -w rasterObjw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)
# 
# plot(raster("rasterObjw2.tif"))
# snet=read.shapefile("rasterObjnet2")
# ns=length(snet$shp$shp)
# for(i in 1:ns)
# {
#   lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
# }
# 
# # Wetness Index
# system("mpiexec -n 8 SlopeAreaRatio -slp rasterObjslp.tif -sca rasterObjsca.tif -sar rasterObjsar.tif", show.output.on.console=F, invisible=F)
# sar=raster("rasterObjsar.tif")
# wi=sar
# wi[,]=-log(sar[,])
# plot(wi)
# 
# # Distance Down
# system("mpiexec -n 8 DinfDistDown -ang rasterObjang.tif -fel rasterObjfel.tif -src rasterObjsrc2.tif -m ave v -dd rasterObjdd.tif",show.output.on.console=F,invisible=F)
# plot(raster("rasterObjdd.tif"))
# 
