### Statistical analysis

data <- read.table("Store.txt", header=T)

### PATCH AND DIVERSITY STATISTICS
library("SDMTools")

# data preparation
years <- time/365
Patch_store <- list()
P_yearly <- list()

for (i in 1:runs){
      #aggregating yearly values
  P_yearly[i] <- aggregate(Store$P, nfrequency =365 , FUN = "mean")
}
for (i in 1:runs){
# calculating patch statistics
P_all_stat <- ClassStat(P_yearly[i] , cellsize = dx, bkgd = NA, latlon = FALSE)
#save results
Patch_store <- as.data_frame(i,P_all_stat[i])
}

# for (j in 1:years){ 
attach(Patch_store)

for (i in 1:runs){

dens <- P_all_stat[i]$patch.density 
shape <- P_all_stat[i]$mean.shape.index 

trend_dens = ma(as.ts(dens[i]), order = 4, centre = F) ## specify parameters...
trend_shape= ma(as.ts(shape[i]), order = 4, centre = F)

if (trend_dens[i] <= 0.01 & trend_shape <= 0.01)
  
  }

  ## find the year where the trend does not change much anymore
  ## and average time - that year and
  

# # P_A$n.patches
# # P_A$mean.patch.area
# # P_A$sd.patch.area
# # P_A$patch.density
# # P_A$mean.shape.index
# # P_A$sd.shape.index
# P_A$total.area 


  for (i in 1:runs){
  
  library("vegan")
  div_mat <- matrix(0, nrow= runs, ncol=4)
  
  for (i in 1:runs){
    
    div_mat[i,1] <- i
    div_mat[i,2]<-sum(P_A_mat)
    div_mat[i,3]<-sum(P_B_mat)
    div_mat[i,4]<-sum(P_C_mat)
  }
  
  shannon <- diversity(div_mat, index = "shannon", MARGIN = 1, base = exp(1))
  simpson <- diversity(div_mat, index = "simpson", MARGIN = 1)
  
  ## how is the div output summed up?


  ## from one 
  ## might need to calculate the patch areas seperately and then check ratio
  
  Store_Analysis[i] <- data_frame(data$simInput[i], P_all_stat[i], P_A_stat[i], P_B_stat[i], P_C_stat[i], SmM_mat[i], shannon[i], simpson[i])
}

write.table(Store_Analysis, "PatchStatDiversity.txt")

# Store <- data_frame(simInput,  
#                     P= array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
#                     P_A = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
#                     P_B = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
#                     P_C = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
#                     SmM = array(matrix(0, nrow=rows, ncol=cols), dim = c(rows, cols, runs)),
#                     Zmatrix, eleve_data_new) 





### Linear regression

linReg <- read.table("PatchStatDiversity.txt")
attach(linReg)

# P_A$n.patches
# P_A$mean.patch.area
# P_A$sd.patch.area
# P_A$patch.density
# P_A$mean.perim.area.ratio
# P_A$sd.perim.area.ratio
# P_A$mean.shape.index
# P_A$sd.shape.index

# Store_Analysis[i] <- data_frame(data$simInput[i], P_all_stat[i], P_A_stat[i], P_B_stat[i], P_C_stat[i], SmM_mat[i], shannon[i], simpson[i])
# }
# Z =Z,ConcConst = ConcConst, CMgw = CMgw,gslp = gslp, 
# dA = dA, k1A = k1A,  b1A = b1A, b2A = b2A, q1A = q1A, q2A = q2A, sigmaPA= sigmaPA,
# dB = dB, k1B = k1B,  b1B = b1B, b2B = b2B, q1B = q1B, q2B = q2B, sigmaPB= sigmaPB,
# dC = dC, k1C = k1C,  b1C = b1C, b2C = b2C, q1C = q1C, q2C = q2C, sigmaPC= sigmaPC,
# ZrA =ZrA, ZrB =ZrB, ZrC =ZrC,
# sigma2 = sigma2, range = range)

# A
Stat_A <- lm(P_A$n.patches ~ )
Plants_Stat <- lm( ~ Z + Zr + d + ConcConst + CM.gw + c + alpha + lambda, data = mydata)
summary(Plants_Stat)



# from summary, read R2, p and correlation and write into table
# 1.	Geometrical analysis: The same patch geometry analysis as described earlier will be applied. Consequently, the resulting patch geometry dataset will be compared with the mono-species version. 
# 
# 2.	The Mann–Whitney–Wilcoxon test shall be conducted to test whether the patch statistics of the multi-species simulations are significantly different from the mono-species runs.

