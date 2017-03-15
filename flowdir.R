
### preparations, put them somewhere else later
flowdir <- ang # ang is angle flowdir[i,j] from DInf TauDEM as raster
flowdir[is.na(flowdir)] <- 8   ###********************The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...

runon_base <- raster(raster)
values(runon_base)<-0
###
runon_fun<-function(flowdir){
##### For Loop to fill runon with runoff from the right cell:
  
for (i in 1:nrow(flowdir)){
  
      for (j in 1:ncol(flowdir)) {
        
      if((i+1)<=10 && (j+1)<=10 && (i-1)>=1 && (j-1)>=1){
        
        
  if(flowdir[i,j]==0)
  {
    runon_base[i,j+1]<-1 #right
    
  }
  if(flowdir[i,j]==(pi/4))
    
  {
    runon_base[i-1,j+1]<-1 #top right
  }
  if(flowdir[i,j]==(pi/2))
    
  {
    runon_base[i-1,j]<-1 #top
  }
  if(flowdir[i,j]==(3*pi/4))
    
  {
    runon_base[i-1,j-1]<-1 #top left
  }
  if(flowdir[i,j]==pi)
    
  {
    runon_base[i,j-1]<-1 #left
  }
  if(flowdir[i,j]==(5*pi/4))
    
  {
    runon_base[i+1,j-1]<-1 #bottom left
  }
  if(flowdir[i,j]==(3*pi/2))
    
  {
    runon_base[i+1,j]<-1 #bottom 
  }

  
  
  if(flowdir[i,j]==(7*pi/4))
    
  {
    runon_base[i+1,j+1]<-1 #bottom tight
  }
  
  if(flowdir[i,j]>0 && flowdir[i,j]<(pi/4))
  {
    if(flowdir[i,j]/((pi/2)-(pi/4))<0.5) {
      runon_base[i,j+1]<- 1*(1-flowdir[i,j]/((pi/4)-0)) ## more to the right
      runon_base[i-1,j+1] <- 1*((flowdir[i,j]/((pi/4)-0)))
    }
    else{
     runon_base[i-1,j+1]<- 1*(1-(flowdir[i,j]/((pi/4)-0))) ##more to the top right
     runon_base[i,j+1]<- 1*((flowdir[i,j]/((pi/4)-0)))
    }}

    
     if(flowdir[i,j]>(pi/4) && flowdir[i,j]<(pi/2))
       {
         if((flowdir[i,j]-(pi/4))/((pi/4)-0)<0.5) {
        runon_base[i-1,j+1]<- 1*((1-flowdir[i,j]-(pi/4))/((pi/4)-0)) ## more to the top right
        runon_base[i-1,j]<- 1*((flowdir[i,j]-(pi/4))/((pi/4)-0))
         }
         else{
           runon_base[i-1,j]<- 1*((1-(flowdir[i,j]-(pi/4))/((pi/4)-0))) ##more to the top 
           runon_base[i-1,j+1]<- 1*((flowdir[i,j]-(pi/4))/((pi/4)-0))
            }}
      
        
          if(flowdir[i,j]>(pi/2) && flowdir[i,j]<(3*pi/4)){
         
           if((flowdir[i,j]-(pi/2))/((pi/4)-0)<0.5) {
              runon_base[i-1,j]<- 1*((1-(flowdir[i,j]-(pi/2)))/((pi/4)-0)) ## more to the top 
              runon_base[i-1,j-1]<- 1*((flowdir[i,j]-(pi/2))/((pi/4)-0))
           }
           else{
              runon_base[i-1,j-1]<- 1*((1-(flowdir[i,j]-(pi/2))/((pi/4)-0))) ##more to the top left
              runon_base[i-1,j]<- 1*((flowdir[i,j]-(pi/2))/((pi/4)-0))
              }}
          
             if(flowdir[i,j]>(3*pi/4) && flowdir[i,j]<pi)
               
              {
                 if((flowdir[i,j]-(3*pi/4))/((pi/4)-0)<0.5) {
                   runon_base[i-1,j-1]<-((1-(flowdir[i,j]-(3*pi/4)))/((pi/4)-0)) ## more to the top left
                   runon_base[i,j-1]<- 1*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                 }
                else{
                    runon_base[i,j-1]<- 1*((1-(flowdir[i,j]-(3*pi/4))/((pi/4)-0))) ##more to the left
                    runon_base[i-1,j-1]<- 1*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                  }}
             
                if(flowdir[i,j]>pi && flowdir[i,j]<(5*pi/4))
                {
                  if((flowdir[i,j]-(pi))/((pi/4)-0)<0.5) {
                    runon_base[i,j-1]<- 1*((1-(flowdir[i,j]-(pi)))/((pi/4)-0)) ## more to the left
                    runon_base[i+1,j-1]<- 1*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }
                  else{
                    runon_base[i+1,j-1]<-1*((1-(flowdir[i,j]-(pi))/((pi/4)-0))) ##more to the bottom left
                    runon_base[i,j-1]<-1*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }}
                  
                  if(flowdir[i,j]>(5*pi/4) && flowdir[i,j]<(3*pi/2))
                  {
                    if((flowdir[i,j]-(5*pi/4))/((pi/4)-0)<0.5) {
                      runon_base[i+1,j-1]<- 1*((1-(flowdir[i,j]-(5*pi/4)))/((pi/4)-0)) ## more to the bottom left
                      runon_base[i+1,j]<- 1*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                    else{
                        runon_base[i+1,j]<- 1*((1-(flowdir[i,j]-(5*pi/4))/((pi/4)-0))) ##more to the bottom
                        runon_base[i+1,j-1]<- 1*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                  }
    
                        
                    if(flowdir[i,j]>(3*pi/2) && flowdir[i,j]<(7*pi/4))
                    {
                      if((flowdir[i,j]-(3*pi/2))/((pi/4)-0)<0.5) {
                         runon_base[i+1,j]<- 1*((1-(flowdir[i,j]-(3*pi/2)))/((pi/4)-0)) ## more to the bottom 
                         runon_base[i+1,j+1]<- 1*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                      else{
                         runon_base[i+1,j+1]<- 1*((1-(flowdir[i,j]-(3*pi/2))/((pi/4)-0))) ##more to the bottom right
                         runon_base[i+1,j]<- 1*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                    }
                      
                      
                      if(flowdir[i,j]>(7*pi/4) && flowdir[i,j]<(2*pi)){
                      
                        if((flowdir[i,j]-(7*pi/4))/((pi/4)-0)<0.5) {
                          runon_base[i+1,j+1]<- 1*((1-(flowdir[i,j]-(7*pi/4)))/((pi/4)-0)) ## more to the bottom right
                          runon_base[i,j+1]<- 1*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        else{
                           runon_base[i,j+1]<- 1*((1-(flowdir[i,j]-(7*pi/4))/((pi/4)-0))) ##more to the right
                           runon_base[i+1,j+1]<- 1*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        
                     
                      }
      }
    }
}

return(runon_base)

        
                   }

rn <-runon_fun(flowdir=flowdir)
rn[is.na(rn)] <- 0
plot(rn)
