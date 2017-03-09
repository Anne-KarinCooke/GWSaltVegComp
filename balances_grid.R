## Source functions
setwd("H:/real salt/Salinization-model")
source("Rainfall.R")
source("Infiltration.R")
source("Flux.R")
source("Vegetation functions.R")
source("Runoff.R")
source("Constants.R")
source("storage formats.R")

#################################################################################################################

timeincr= 1/deltat

##################################################################################################################
##################################################################################################################
#  BALANCES FUNTION STARTS

balances2D <- function(Rain, par,
                       soilpar,
                       vegpar){ 
  
 #plotit=F,
  
  for (i in nrow(raster)) { 
    
    for (j in ncol(raster)){
    
    for (t in 2:length(Rain)){
      
      for (tt in 1:(deltat-1)) {
        
        
        
        h.old[i,j] <- ifelse(tt==1,h[i,j][t-1],h_sub[i,j][tt]) #raster
        
        P.old[i][j] <- ifelse(tt==1,P[i][j][t-1],P_sub[i][j][tt])  # arrays
        M.old[i][j] <- ifelse(tt==1,M[i][j][t-1],M_sub[i][j][tt])
        SmI.old[i][j] <-ifelse(tt==1,SmI[i][j][t-1],SmI_sub[i][j][tt])
        CM.old[i][j] <-ifelse(tt==1,CM[i][j][t-1],CM_sub[i][j][tt])
        Svir.old[i][j] <-ifelse(tt==1,Svir[i][j][t-1],Svir_sub[i][j][tt])
        
        q_sub[i,j][tt]<- OF(h=h.old[i,j], soilpar=soilpar,slope=slp)
        
        runon_sub[i,j][tt]<-runon_fun(flowdir=flowdir,q=q_sub[i,j][tt])
        
        
        #### how to define runon correctly
        
        h_sub[i][j][tt+1] <- h.old[i,j] + ifelse(tt==1,(10*Rain[t]),0) - Infil(h.old[i][j], P.old[i][j],par)*timeincr - q_sub[i,j][tt] + runon[i,j][tt]
        

        # Infiltration
        par$alpha_i <- ifelse(h_sub[i,j][tt+1]<soilpar$K_s*timeincr, 1,(1-(h_sub[i,j][tt+1]-soilpar$K_s*timeincr)/h_sub[i,j][tt+1]))
        I_sub[i][j][tt] <- Infil(h.old[i,j], P.old[i][j],par)*timeincr
       
        
        #  1. Update soil moisture with infiltration
        
        M_sub[i][j][tt+1] <- M.old[i,j] + I_sub[i][j][tt] #+ Diff[tt,g-1]     # plus soil moisture diffusing from grid cell higher up  
        
        # Now do all plant uptake and growth
        # water uptake by plants: include infiltration in available water
        
        
        WU_sub[i][j][tt] <- WU(M_sub[i][j][tt+1],P.old[i,j],par)*timeincr 
        
        # growth rate
        Gr_sub[i][j][tt] <- Gr(M=Svir.old[i,j], P.old[i,j],par)*timeincr 
        # Mortality
        Mo_sub[i][j][tt]<- Mo(P.old[i,j], M=M.old[i,j], Svir=Svir.old[i,j], par)*timeincr
        
        
        # calculate plant biomass balance
        P_sub[i][j][tt+1] <- P.old[i,j] + Gr_sub[i][j][tt]- Mo_sub[i][j][tt] 
        
        
        # re-calculate water balance
        # 2. before leaching
        M_sub[i][j][tt+1] <- M.old[i,j] + I_sub[i][j][tt] - WU_sub[i][j][tt] #- L_sub[tt] 
        
        
        
        # 3. calculate leaching and capillary rise amount
        flux_sub[i][j][tt+1]<-do.call(L_n,list(M=M_sub[i][j][tt+1],Z=Zras[i,j],soilpar=soilpar,vegpar=vegpar))
        
        
        # Divergence, soil moisture diffusion, still to be fixed
        ##Diff[i][j][tt]<- M_sub[i][j][tt+1]*(Dm*timeincr)*((Z[g+1]-Z[i,j])/dist)  
        
        
        # 4. final adjust soil moisture for leaching/rise AND DIVERGENCE
        M_sub[i][j][tt+1] <- M_sub[i][j][tt+1] + flux_sub[i][j][tt+1]*timeincr #- Diff[i][j][tt] 
        
        
        # calculate saltbalance
        
        
        # Salt leaching
        L_salt[i][j][tt+1] <- ifelse(flux_sub[i][j][tt+1]<0, par$f*CM_sub[i][j][tt+1]*flux_sub[i][j][tt+1]*timeincr,0) # leaching of salt
        
        # salt uplfow
        U_salt[i][j][tt+1] <- ifelse(flux_sub[i][j][tt+1]>0, par$CM.gw*flux_sub[i][j][tt+1]*timeincr,0) # rise of salt
        
        # salt mass coming in with infiltration
        SmI_sub[i][j][tt+1]<- SmI.old[i,j] + I_sub[i][j][tt]*par$ConcConst 
        
        #salt mass in soil
        SmM_sub[i][j][tt+1] <- SmI_sub[i][j][tt+1] + U_salt[i][j][tt+1] - L_salt[i][j][tt+1]
        
        # salt concentration in soil
        CM_sub[i][j][tt+1]<- (SmM_sub[i][j][tt+1]/M_sub[i][j][tt+1])*(1/58.44)         
        
        # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
        Svir_sub[i][j][tt+1]<-soilpar$n*vegpar$Zr*((soilpar$h1bar*10^-1)^(1/soilpar$b))*
          ((soilpar$h1bar*10^-1)*(M_sub[i][j][tt+1]/
                                    (soilpar$n*vegpar$Zr))^(-soilpar$b)
           +(3.6*CM_sub[i][j][tt+1]))^(-1/soilpar$b)
        
        # checking the mass balance!
        mb_sub[i][j][tt] <- I_sub[i][j][tt] - WU_sub[i][j][tt] + flux_sub[i][j][tt] - q_sub[i,j][tt]
        
      } 
      
      # Aggregating the substep results to daily values.
      
      P[i][j][t] = P_sub[i][j][deltat]
      M[i][j][t] = M_sub[i][j][deltat]
      h[i,j][t] = h_sub[i,j][deltat] #### modified
      CM[i][j][t] = CM_sub[i][j][deltat]
      #      SmM[t,g] = SmI[t,g] = SmM_sub[deltat,g]  ###
      SmI[i][j][t] = SmI_sub[i][j][deltat]
      SmM[i][j][t] = SmM_sub[i][j][deltat]
      In[i][j][t]= sum(I_sub[i][j])
      Svir[i][j][t] = Svir_sub[i][j][deltat]
      flux[i][j][t]= sum(flux_sub[i][j])
      q[i,j][t] = sum(q_sub[i,j]) ####modified
      runon[i,j][t] = sum(runon_sub[i,j])
      mb[i][j][t] = sum(mb_sub[i][j])
      
      
    }
    
    
    # Plotting
    
#     if (plotit==T) {  
#       plot(M[,g], type="l",ylim=c(0,100),xlim=c(0,time),xlab=("time [d]"), main=paste("lambda=", lambda[j],"alpha=", alpha[i], "gridcell=", grid[g]))
#       points(Rain*10, type="h", col="skyblue")
#       
#       abline(h=0, col="Gray50",lwd=2,lty=2)
#       
#       lines(SmM[,g],type="l", col="red")
#       lines(CM[,g],type="l", col="purple")
#       lines(P[,g],type="l", col="green")
#       
      
      #  legend("topright",cex=1, pt.cex=0.4, c("Moisture [mm]","Rainfall [mm]*10","overland flow depth[mm] ","salt mass in soil water [g]", "salt concentration in soil water [g/l]", "Plant biomass density [g/m^2]"),
      #           col=c("black","skyblue","blue","red","purple","green"),lty=1)
      #  
      
    }
  }
  Out <- data.frame(P=P,M=M,h=h, CM=CM, SmM=SmM, In=In, flux=flux, Svir=Svir,h=h, q=q, mb=mb)  ### *** CHECK AGAIN IF DATAFRAME IS RIGHT FORMAT
  return(Out)
}
}
