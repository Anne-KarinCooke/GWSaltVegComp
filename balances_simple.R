balances <- function(Rain, par, plotit=F,
                     soilpar,
                     vegpar){
  
  
  
  # Storage vectors for the daily steps are initialized.
  
  M <- rep(0,length(Rain)) # soil moisture [mm]
  h <- rep(0,length(Rain)) # infiltration depth [mm]
  P <- rep(0,length(Rain)) #biomass density []
  CM<- rep(0,length(Rain)) # Salt concentration in soil water in g/L or g/mm
  SmI<- rep(0,length(Rain)) # Salt mass in infiltrating water [g]
  SmM <- rep(0,length(Rain)) # Salt mass in soil water [g]
  In <- rep(0,length(Rain)) # infiltration [mm]
  
  Svir <- rep(0,length(Rain)) # virtual saturation
  flux<- rep(0,length(Rain)) # drainage and capillary rise flux, according to sign
  
  #mass balance
  mb <- rep(0,length(Rain))
  
  
  
  # Initial values to start the simulation.
  
  
  M[1] <- 10
  h[1] <- 10 
  P[1] <- 30
  CM[1]<- 0
  Svir[1] <- M[1]
  
  
  
  # We decided to split the numerical calculations for the daily into 12 substeps.
  
  deltat <- 12 # split in 12 increments
  
  
  # Storage vectors for the substeps are initialized.
  
  M_sub <- rep(0,deltat)
  h_sub <- rep(0,deltat)
  I_sub <- rep(0,deltat)
  WU_sub <-rep(0,deltat) # Water uptake in mm
  P_sub <- rep(0,deltat) 
  Gr_sub <- rep(0,deltat) # Growth of biomass
  Mo_sub<- rep(0,deltat) # Mortality of biomass
  SmI_sub <- rep(0,deltat) 
  SmM_sub<- rep(0,deltat) 
  CM_sub<- rep(0,deltat) 
  Svir_sub <- rep(0,deltat) # virtual saturation
  
  flux_sub<-rep(0,deltat) # calculates leakage loss without evaporation loss
  
  U_salt <-rep(0,deltat) # Salt mass rising
  L_salt <-rep(0,deltat) # salt mass drained
  
  #mass balance
  mb_sub <-rep(0,deltat)
  
  
  timeincr= 1/deltat
  
  for (t in 2:length(Rain)){
    
    for (tt in 1:(deltat-1)) {
      
      h.old <- ifelse(tt==1,h[t-1],h_sub[tt])
      P.old <- ifelse(tt==1,P[t-1],P_sub[tt])
      M.old <- ifelse(tt==1,M[t-1],M_sub[tt])
      SmI.old <-ifelse(tt==1,SmI[t-1],SmI_sub[tt])
      CM.old <-ifelse(tt==1,CM[t-1],CM_sub[tt])
      Svir.old <-ifelse(tt==1,Svir[t-1],Svir_sub[tt])
      
      
      # Balance for water depth on soil
      h_sub[tt+1] <- h.old + ifelse(tt==1,(10*Rain[t]),0) 
      #- Infil(h.old, P.old,par)*timeincr
      
      # Infiltration
      par$alpha_i <- ifelse(h_sub[tt+1]<soilpar$K_s*timeincr, 1,
                            (1-(h_sub[tt+1]-soilpar$K_s*timeincr)/h_sub[tt+1]))
      # Calculate infiltration and recalculate h_sub   
      I_sub[tt] <- Infil(h.old, P.old,par)*timeincr
      h_sub[tt+1] <- h_sub[tt+1] - I_sub[tt] 
      
      # Now do all plant uptake and growth
      # water uptake by plants: include infiltration in available water
      
      WU_sub[tt] <- WU(M=Svir.old,P.old,par)*timeincr 
      
      # growth rate
      Gr_sub[tt] <- Gr(M=Svir.old, P.old,par)*timeincr 
      # Mortality
      Mo_sub[tt]<- Mo(P.old,M=M.old, Svir=Svir.old, par)*timeincr
      # calculate plant biomass balance
      P_sub[tt + 1] <- P.old + Gr_sub[tt]- Mo_sub[tt] 
      
      # 1. Update soil moisture with infiltration
      M_sub[tt + 1] <- M.old + I_sub[tt] - WU_sub[tt]      
      
      # 2. calculate leaching and capillary rise amount
      flux_sub[tt+1]<-do.call(L_n,list(M=M_sub[tt+1],Z=Z,soilpar=soilpar,vegpar=vegpar))
      
      # 3. final adjust soil moisture for leaching or capillary rise
      
      M_sub[tt + 1] <-  M_sub[tt + 1] + flux_sub[tt+1]*timeincr
      
      
      # calculate saltbalance
      
      # Salt leaching
      L_salt[tt+1] <- ifelse(flux_sub[tt+1]<0,
                             par$f*CM_sub[tt+1]*flux_sub[tt+1]*timeincr,0)
      # leaching of salt
      
      # salt upflow
      U_salt[tt+1] <- ifelse(flux_sub[tt+1]>0, 
                             par$CM.gw*flux_sub[tt+1]*timeincr,0) 
      # rise of salt
      
      # salt mass coming in with infiltration
      SmI_sub[tt+1]<- SmI.old + I_sub[tt]*par$ConcConst 
      
      #salt mass in soil
      SmM_sub[tt+1] <- SmI_sub[tt+1] + U_salt[tt+1] - L_salt[tt+1]
      # Calculate salt concentration in the soil
      CM_sub[tt+1]<- (SmM_sub[tt+1]/M_sub[tt+1])*(1/58.44) # 
      
      # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
      Svir_sub[tt + 1]<-soilpar$n*vegpar$Zr*((soilpar$h1bar*10^-1)^(1/soilpar$b))*
        ((soilpar$h1bar*10^-1)*(M_sub[tt + 1]/
                                  (soilpar$n*vegpar$Zr))^(-soilpar$b)
         +(3.6*CM_sub[tt + 1]))^(-1/soilpar$b)
      
      mb_sub[tt] <- I_sub[tt] - WU_sub[tt] + flux_sub[tt]
      
    } 
    
    # Aggregating the substep results to daily values.
    
    P[t] = P_sub[deltat]
    M[t] = M_sub[deltat]
    h[t] = h_sub[deltat]
    CM[t] = CM_sub[deltat]
    SmM[t] = SmM_sub[deltat] 
    SmI[t]=SmI_sub[deltat]
    In[t]= sum(I_sub)
    flux[t] = sum(flux_sub)
    Svir[t] = Svir_sub[deltat]
    mb[t] = sum(mb_sub[deltat])
    
  }
  
  
  # Plotting
  
  if (plotit==T) {  
    plot(M, type="l",ylim=c(-10,100),xlim=c(0,time),xlab=("time [d]"), main=paste(" lambda=", lambda[j],"alpha=", alpha[i]))
    points(Rain*10, type="h", col="skyblue")
    
    lines(h,type="l", col="blue")
    abline(h=0, col="Gray50",lwd=2,lty=2)
    
    lines(mb*100,type="l", col="red")
    lines(CM,type="l", col="purple")
    lines(P/10,type="l", col="green")
    lines(flux,type="l", col="orange")
    
    
    #  legend("topright", title="Variables",cex=1, pt.cex=0.2, c("Moisture [mm]","Rainfall [mm]*10","overland flow depth[mm] ","salt mass in soil water [g]", "salt concentration in soil water [g/l]", "Plant biomass density [g/m^2]/10"),
    #          col=c("black","skyblue","blue","red","purple","green"),lty=1)
    #  
  }
  
  Out <- data.frame(P=P,M=M,h=h, CM=CM, SmM=SmM, In=In, flux=flux, Svir=Svir, mb=mb)
  return(Out)
}