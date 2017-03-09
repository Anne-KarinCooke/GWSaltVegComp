# for (k in ncell(raster)){

source("Runoff.R")

#   OF<- function(h, soilpar=soilpar, slope=slp) #still need to changed soilpar to slope
#     
#     qq=(soilpar$cn/soilpar$Mn)*h^(5/3)*slope^(1/2)
# 
# return(qq)
# }


### Rewrite q
##############q_sub[tt,g]<-OF(h=h.old[g],soilpar=soilpar) 

### q stuff is in balances code
# q<-brick(raster,nl=length(Rain))
# values(q)<-0
# q<-brick(raster,nl=length(deltat))
# values(q)<-0
# 
# q_sub[i,j][tt]<- OF(h=h.old[i,j], soilpar=soilpar,slope=slp)


### define runon, from the right cell to the following...######RUNON

# runon<-brick(raster,nl=length(Rain))
# 
# runon_sub<-brick(raster,nl=length(deltat))
# values(q)<-0


# for (k in nrow(raster)){
#   
#   for (l in ncol(raster)){
#     
runon<-raster(raster)
    q<-raster(raster)
values(q) <- runif(ncell(q),0,5)


### preparations, put them somewhere else later
flowdir <- ang # ang is angle flowdir[i,j] from DInf TauDEM as raster
flowdir[is.na(flowdir)] <- 8   ###********************The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...


###
runon_fun<-function(flowdir,q){


##### For Loop to fill runon with runoff from the right cell:
for (i in nrow(flowdir)){
  
      for (j in ncol(flowdir)) {
        
  if(flowdir[i,j]==0)
  {
    runon[i,j+1]<-q[i,j] #right

  }
  if(flowdir[i,j]==(pi/4))
    
  {
    runon[i-1,j+1]<-q[i,j] #top right
  }
  if(flowdir[i,j]==(pi/2))
    
  {
    runon[i-1,j]<-q[i,j] #top
  }
  if(flowdir[i,j]==(3*pi/4))
    
  {
    runon[i-1,j-1]<-q[i,j] #top left
  }
  if(flowdir[i,j]==pi)
    
  {
    runon[i,j-1]<-q[i,j] #left
  }
  if(flowdir[i,j]==(5*pi/4))
    
  {
    runon[i+1,j-1]<-q[i,j] #bottom left
  }
  if(flowdir[i,j]==(3*pi/2))
    
  {
    runon[i+1,j]<-q[i,j] #bottom 
  }
  #***
      }}
  
  
  if(flowdir[i,j]==(7*pi/4))
    
  {
    runon[i+1,j+1]<-q[i,j] #bottom tight
  }
  
  if(flowdir[i,j]>0 && flowdir[i,j]<(pi/4))
  {
    if(flowdir[i,j]/((pi/2)-(pi/4))<0.5) {
      runon[i,j+1]<-q[i,j]*(1-flowdir[i,j]/((pi/4)-0)) ## more to the right
      runon[i-1,j+1] <- q[i,j]*((flowdir[i,j]/((pi/4)-0)))
    }
    else{
     runon[i-1,j+1]<-q[i,j]*(1-(flowdir[i,j]/((pi/4)-0))) ##more to the top right
     runon[i,j+1]<- q[i,j]*((flowdir[i,j]/((pi/4)-0)))
    }}

    
     if(flowdir[i,j]>(pi/4) && flowdir[i,j]<(pi/2))
       {
         if((flowdir[i,j]-(pi/4))/((pi/4)-0)<0.5) {
        runon[i-1,j+1]<-q[i,j]*((1-flowdir[i,j]-(pi/4))/((pi/4)-0)) ## more to the top right
        runon[i-1,j]<-q[i,j]*((flowdir[i,j]-(pi/4))/((pi/4)-0))
         }
         else{
           runon[i-1,j]<-q[i,j]*((1-(flowdir[i,j]-(pi/4))/((pi/4)-0))) ##more to the top 
           runon[i-1,j+1]<-q[i,j]*((flowdir[i,j]-(pi/4))/((pi/4)-0))
            }}
      
        
          if(flowdir[i,j]>(pi/2) && flowdir[i,j]<(3*pi/4)){
         
           if((flowdir[i,j]-(pi/2))/((pi/4)-0)<0.5) {
              runon[i-1,j]<-q[i,j]*((1-(flowdir[i,j]-(pi/2)))/((pi/4)-0)) ## more to the top 
              runon[i-1,j-1]<-q[i,j]*((flowdir[i,j]-(pi/2))/((pi/4)-0))
           }
           else{
              runon[i-1,j-1]<-q[i,j]*((1-(flowdir[i,j]-(pi/2))/((pi/4)-0))) ##more to the top left
              runon[i-1,j]<-q[i,j]*((flowdir[i,j]-(pi/2))/((pi/4)-0))
              }}
          
             if(flowdir[i,j]>(3*pi/4) && flowdir[i,j]<pi)
               
              {
                 if((flowdir[i,j]-(3*pi/4))/((pi/4)-0)<0.5) {
                   runon[i-1,j-1]<-((1-(flowdir[i,j]-(3*pi/4)))/((pi/4)-0)) ## more to the top left
                   runon[i,j-1]<-q[i,j]*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                 }
                else{
                    runon[i,j-1]<-q[i,j]*((1-(flowdir[i,j]-(3*pi/4))/((pi/4)-0))) ##more to the left
                    runon[i-1,j-1]<-q[i,j]*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                  }}
             
                if(flowdir[i,j]>pi && flowdir[i,j]<(5*pi/4))
                {
                  if((flowdir[i,j]-(pi))/((pi/4)-0)<0.5) {
                    runon[i,j-1]<-q[i,j]*((1-(flowdir[i,j]-(pi)))/((pi/4)-0)) ## more to the left
                    runon[i+1,j-1]<-q[i,j]*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }
                  else{
                    runon[i+1,j-1]<-q[i,j]*((1-(flowdir[i,j]-(pi))/((pi/4)-0))) ##more to the bottom left
                    runon[i,j-1]<-q[i,j]*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }}
                  
                  if(flowdir[i,j]>(5*pi/4) && flowdir[i,j]<(3*pi/2))
                  {
                    if((flowdir[i,j]-(5*pi/4))/((pi/4)-0)<0.5) {
                      runon[i+1,j-1]<-q[i,j]*((1-(flowdir[i,j]-(5*pi/4)))/((pi/4)-0)) ## more to the bottom left
                      runon[i+1,j]<-q[i,j]*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                    else{
                        runon[i+1,j]<-q[i,j]*((1-(flowdir[i,j]-(5*pi/4))/((pi/4)-0))) ##more to the bottom
                        runon[i+1,j-1]<-q[i,j]*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                  }
    
                        
                    if(flowdir[i,j]>(3*pi/2) && flowdir[i,j]<(7*pi/4))
                    {
                      if((flowdir[i,j]-(3*pi/2))/((pi/4)-0)<0.5) {
                         runon[i+1,j]<-q[i,j]*((1-(flowdir[i,j]-(3*pi/2)))/((pi/4)-0)) ## more to the bottom 
                         runon[i+1,j+1]<-q[i,j]*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                      else{
                         runon[i+1,j+1]<-q[i,j]*((1-(flowdir[i,j]-(3*pi/2))/((pi/4)-0))) ##more to the bottom right
                         runon[i+1,j]<-q[i,j]*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                    }
                      
                      
                      if(flowdir[i,j]>(7*pi/4) && flowdir[i,j]<(2*pi)){
                      
                        if((flowdir[i,j]-(7*pi/4))/((pi/4)-0)<0.5) {
                          runon[i+1,j+1]<-q[i,j]*((1-(flowdir[i,j]-(7*pi/4)))/((pi/4)-0)) ## more to the bottom right
                          runon[i,j+1]<-q[i,j]*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        else{
                           runon[i,j+1]<-q[i,j]*((1-(flowdir[i,j]-(7*pi/4))/((pi/4)-0))) ##more to the right
                           runon[i+1,j+1]<-q[i,j]*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        if(flowdir[i,j]==8){
                          runon[i,j]<-1
                          
                        }
                     
                      }
return(runon)
  
        
                   }
                

# plot(runon_fun(flowdir=flowdir,q=q))

