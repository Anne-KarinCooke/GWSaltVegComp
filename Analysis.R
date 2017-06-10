### Statistical analysis

data <- read.table("Store.txt")

for 

## Check stable state condition!!!
## Check whether standard dev 
startEq ## start of equilibirum condition

# average values from that pint in time onwards

# make new table with P_A, P_B, P_C, SmM averaged values (just one matrix each)

### PATCH AND DIVERSITY STATISTICS


 library("SDMTools")

# get a matrix for each P abc and salt: P_A_mat

#P_A  matrix with P values of species A (average of stable state)

P_A <- ClassStat(P_A_mat, cellsize = dx, bkgd = NA, latlon = FALSE)
P_B <- ClassStat(P_B_mat, cellsize = dx, bkgd = NA, latlon = FALSE)
P_C <- ClassStat(P_C_mat, cellsize = dx, bkgd = NA, latlon = FALSE)

# P_A$n.patches
# P_A$mean.patch.area
# P_A$sd.patch.area
# P_A$patch.density
# P_A$mean.perim.area.ratio
# P_A$sd.perim.area.ratio
# P_A$mean.shape.index
# P_A$sd.shape.index

Store_Analysis <- data_frame(P_A, P_B, P_C, SmM, )


# do geostat with SmM

library("vegan")
data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")

shannon <- diversity(x, index = "shannon", MARGIN = 1, base = exp(1)) ## whole matrix not just rows? 
simpson <- diversity(x, index = "simpson", MARGIN = 1)

# how to define x, so that it includes all the 3 species classes?




# f1( 0 ) = h;
# f1( 1 ) = q;
# f1( 2 ) = In;
# f1( 3 ) = runon;
# f1( 4 ) = Wu;
# f1( 5 ) = P;
# f1( 6 ) = flux;
# f1( 7 ) = M;
# f1( 8 ) = SmM;
# f1( 9 ) = CM;
# f1( 10 ) = mb;
# f1( 11 ) = Svir;
# f1( 12 ) = SmI;
# f1( 13 ) = Smh;
# f1( 14 ) = Ch;
# f1( 15 ) = qsd;
# f1( 16 ) = runonsd;
# f1( 17 ) = seep;
# f1( 18 ) = runonSubsSalt;
# f1( 19 ) = salt_runon;
# f1( 20 ) = Subsrunon;
# f1( 21 ) = mb_salt;





### Linear regression

for (i in 1:runs) {
  
  Plants_Stat <- lm(Plants ~ Z + Zr + d + ConcConst + CM.gw + c + alpha + lambda, data = mydata)
  summary(Plants_Stat)
  
}

# from summary, read R2, p and correlation and write into table

