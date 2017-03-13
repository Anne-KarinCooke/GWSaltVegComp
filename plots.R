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
