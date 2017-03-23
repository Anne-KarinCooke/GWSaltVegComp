###
######### Climate loop for grid model#####################################################################
##########################################################################################################
##########################################################################################################

# SOURCE
setwd("H:/Thesis project model/R project/GWSaltVegComp")


### RUN TIME (how many days)
time <- 20
### Raster size
ext <- 40 ## EXTEND of PLOT in [m]
rows <- 5## rows of cells
cols <- 5 ## columns of cells
# Discretization, subdaily timesteps
deltat<-12

### Sourcing the grid, taudem etc, generates flowdir (ang) and slope raster
source("Rasters.R")
### Source the runon raster (flowdir)
source("flowdir.R")
### Sourcing the rainfall
source("Rainfall.R")


#Runon raster generated from flowdir.R
rn_matrix<- as.matrix(rn,nrow= nrow(rn), ncol=ncol(rn))

# slope raster generated from Rasters.R
slp_matrix<- as.matrix(slp,nrow= nrow(slp), ncol=ncol(slp))

#### RAINFALL GENERATION

##UNIFORM
Rain_function<-function(time){ Rain <- rep(1, time)
return(Rain)}
Rain <-list(Rain_function(time=time)) # Rain as a list

## HETEROGENEOUS
# set.seed(100)
# alpha <- c(0.6,1.5) #  alpha <- seq(0.6,1.5,by=0.1) 
# lambda <- c(0.1,1) #  lambda <- seq(0.1,1,by=0.1)
# delta <- 0
# Rain <- Precip(time,alpha[k],lambda[l],delta)

# ## RAINFALL GENERATION
# for (k in 1:length(alpha)) {
#   
#   for (l in 1:length(lambda)) {
#     # generate the rainfall
#     Rain <- Precip(time,alpha[k],lambda[l],delta)
#     Rainlist <- data.frame(Precip(time,alpha[k],lambda[l],delta))
#   }}


# soilpar["Mn"] = 10.0; // Manning's n
#         soilpar["cn"] = 1.0; // runoff conversion factor
# soilpar["alpha_i"] = 1.0;//#maximum infiltration rate per day, This needs to be a fraction of h (p117 Saco and Moreno-Las Heras)

install.packages("Rcpp")
library(Rcpp)
# Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
# Sys.setenv("PKG_LIBS"="-fopenmp")
install.packages("devtools")
library(devtools)

# Sourcing the cpp functions that define the constants for soil, veg and salt
sourceCpp("soilfun.cpp")
sourceCpp("vegfun.cpp")
sourceCpp("saltfun.cpp")


# creating parameter lists
soilpar <- Soil_cpp("S Clay Loam")
vegpar <-Veg_cpp("Fantasy Tree")
saltpar<- Salt_cpp("Groundwater")


Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("Allfunctions.cpp")
evalCpp("Allfunctions.cpp")
balances2D

sourceCpp("Runoff.cpp")

sourceCpp("Flux.cpp")
evalCpp("Runoff.cpp")
#include "storageformats.hpp"

sourceCpp("Infiltration.cpp")



#include "Flux.hpp"
#include "Runoff.hpp"
#include "Infiltration.hpp"









### cREATING DATA STORAGE OF MODEL OUTPUTS
Store <- list()
                                   
    system.time(                                
    Store<-list(balances2D(Rain=Rain, par=par, soilpar=soilpar, vegpar=vegpar))
    )
    
 Store[[1]]$mb
  mb[20,20,]
  
  
 datamb<-data.frame(Store[[1]]$mb)
    
mb<-brick(Store[[1]]$mb)

levelplot(mb, main="mb")
    
#write.table(Store,"10x10Raster20DaysDeltat12.txt")
    
# Plants
p<-brick(Store[[1]]$P)
# plot(p)

#SOil moisture
sm<-brick(Store[[1]]$M)
# plot(sm)

#h
hh<-brick(Store[[1]]$h)
# plot(hh)

# CM
cm<-brick(Store[[1]]$CM)
# plot(cm)

#flux
fl<-brick(Store[[1]]$flux)

#runoff
qr<-brick(Store[[1]]$q)
levelplot(qr,main="Runoff", sub="20 days, deltat 12, 10*10 raster")

#runon
rr<-brick(Store[[1]]$runon)
levelplot(rr,main="Runon", sub="20 days, deltat 12, 10*10 raster")


library(rasterVis)
levelplot(cm, main="Salt concentration in soil", sub="20 days, deltat 12, 10*10 raster")
levelplot(p,main="Plant biomass density", sub="20 days, deltat 12, 10*10 raster")
levelplot(sm,main="Soil moisture", sub="20 days, deltat 12, 10*10 raster")
levelplot(fl,main="vertical water flux (rise/drainage)", sub="20 days, deltat 12, 10*10 raster")
                


#       sub_store[[j]] <-data.frame(alpha_o=rep(alpha[k],time),
#                                   lambda_o=rep(lambda[l],time),
#                                   
#                                   
#                                   
#                                   balances2D(Rain, par=par, soilpar=soilpar, #plotit=T,
#                                              vegpar=vegpar)) 
      
      
  #     
  #   }
  #   Store[[i]] <- sub_store
  # }
#   gstore[[g]] <-Store 
# }




## Results 


#Plotting M and P for different lambdas
# 
# 
# require(ggplot2)
# lambda_sum <- do.call(rbind,Store[][[2]])
# lambda_sum$time <- rep(1:time,length(lambda))
# 
# # Gridcell 1
# pl <- ggplot(lambda_sum,aes(x=time,y=P.1, colour="P (plant biomass density [gm-2])")) + geom_line()
# pl  <- pl + geom_line(aes(x=time,y=M.1, colour="Moisture"))  
# pl <- pl +  facet_wrap(~lambda_o)   
# pl  + ggtitle("Plant biomass P and soilmoisture M for varying lambdas") +  geom_line(aes(x=time, y=SmM.1, colour= "S (soil salt mg/L")) + theme(plot.title = element_text(lineheight=.4))
# 
# #Gridcell 2
# pl <- ggplot(lambda_sum,aes(x=time,y=P.2, colour="P (plant biomass density [gm-2])")) + geom_line()
# pl  <- pl + geom_line(aes(x=time,y=M.2, colour="Moisture"))  
# pl <- pl +  facet_wrap(~lambda_o)   
# pl  + ggtitle("Plant biomass P and soilmoisture M for varying lambdas") +  geom_line(aes(x=time, y=SmM.2, colour= "S (soil salt mg/L")) + theme(plot.title = element_text(lineheight=.4))
# 
# #Gricell 3
# pl <- ggplot(lambda_sum,aes(x=time,y=P.3, colour="P (plant biomass density [gm-2])")) + geom_line()
# pl  <- pl + geom_line(aes(x=time,y=M.3, colour="Moisture"))  
# pl <- pl +  facet_wrap(~lambda_o)   
# pl  + ggtitle("Plant biomass P and soilmoisture M for varying lambdas") +  geom_line(aes(x=time, y=SmM.3, colour= "S (soil salt mg/L")) + theme(plot.title = element_text(lineheight=.4))
# 
# #Gricell 4
# pl <- ggplot(lambda_sum,aes(x=time,y=P.4, colour="P (plant biomass density [gm-2])")) + geom_line()
# pl  <- pl + geom_line(aes(x=time,y=M.4, colour="Moisture"))  
# pl <- pl +  facet_wrap(~lambda_o)   
# pl  + ggtitle("Plant biomass P and soilmoisture M for varying lambdas") +  geom_line(aes(x=time, y=SmM.4, colour= "S (soil salt mg/L")) + theme(plot.title = element_text(lineheight=.4))
# 
# 
# 
# 
# 
# 
# #Plotting Soil water salt concentration for different alphas and lambdas
# 
# 
# lambda_sum_all <- do.call(rbind,do.call(rbind,Store))
# lambda_sum_all$time <- rep(rep(1:time,length(lambda)),length(alpha))
# 
# 
# pa <- ggplot(lambda_sum_all,aes(x=time,y=SmM,col=as.factor(lambda_o))) + geom_line(linetype=1) 
# # pa <- pa + + scale_color_gradient(low="blue", high="red")
# pa <- pa  + facet_wrap(~alpha_o) + ggtitle("Soil water salt mass [g] for varying alphas and lambdas") + theme(plot.title = element_text(lineheight=.4))
# pa  
# 
# 
# 
# lambda_sum$mb
# 
# 
# 
# 
# # Plotting M and P for different lambdas
# # 
# # ```{r}
# # require(ggplot2)
# # lambda_sum <- do.call(rbind,Store[][[2]])
# # lambda_sum$time <- rep(1:time,length(lambda))
# # 
# # pl <- ggplot(lambda_sum,aes(x=time,y=P, colour="P (plant biomass density [gm-2])")) + geom_line()
# # pl  <- pl + geom_line(aes(x=time,y=M, colour="Moisture"))  
# # pl <- pl +  facet_wrap(~lambda_o, ncol=2)   #, colour=lambda_o (put this in aes-brackets) 
# # pl  + ggtitle("Plant biomass P and soilmoisture M for varying lambdas") +  geom_line(aes(x=time, y=SmM, colour= "S (soil salt mg/L")) + theme(plot.title = element_text(lineheight=.8, face="bold"))
# # ```
# # 
# # Plotting Soil water salt concentration for different alphas and lambdas
# # 
# # ```{r}
# # lambda_sum_all <- do.call(rbind,do.call(rbind,Store))
# # lambda_sum_all$time <- rep(rep(1:time,length(lambda)),length(alpha))
# # 
# # 
# # pa <- ggplot(lambda_sum_all,aes(x=time,y=SmM,col=lambda_o)) + geom_line(linetype=1) 
# # pa <- pa + scale_color_gradient(low="blue", high="red") + facet_wrap(~alpha_o, ncol=2) + ggtitle("Soil water salt mass [g] for varying alphas and lambdas") + theme(plot.title = element_text(lineheight=.8, face="bold"))
# # pa  
# # ```

# Comparing results

# Biomass without salt
# summary(lambda_sum_woSalt$P)
# # Biomass with salt
# summary(lambda_sum$P)
# #Soil moisture without salt
# # summary(lambda_sum_woSalt$M)
# #Soil moisture with salt
# summary(lambda_sum$M)
# summary(lambda_sum$Svir)
# # summary(1+lambda_sum$Svir)
# # summary(-lambda_sum$Svir*n*Zr)
# 
# 


