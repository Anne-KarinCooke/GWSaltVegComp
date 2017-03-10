
### preparations, put them somewhere else later
flowdir <- ang # ang is angle flowdir[i,j] from DInf TauDEM as raster
flowdir[is.na(flowdir)] <- 8   ###********************The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...

blabla <-raster(raster)
values(blabla) <- runif(ncell(raster),0,1000)
blabla
runon<-raster(raster)
###
runon_fun<-function(flowdir,qf){
##### For Loop to fill runon with runoff from the right cell:
  
for (i in 1:nrow(flowdir)){
  
      for (j in 1:ncol(flowdir)) {
        
        
  if(flowdir[i,j]==0)
  {
    runon[i,j+1]<-qf #right
    
  }
  if(flowdir[i,j]==(pi/4))
    
  {
    runon[i-1,j+1]<-qf #top right
  }
  if(flowdir[i,j]==(pi/2))
    
  {
    runon[i-1,j]<-qf #top
  }
  if(flowdir[i,j]==(3*pi/4))
    
  {
    runon[i-1,j-1]<-qf #top left
  }
  if(flowdir[i,j]==pi)
    
  {
    runon[i,j-1]<-qf #left
  }
  if(flowdir[i,j]==(5*pi/4))
    
  {
    runon[i+1,j-1]<-qf #bottom left
  }
  if(flowdir[i,j]==(3*pi/2))
    
  {
    runon[i+1,j]<-qf #bottom 
  }

  
  
  if(flowdir[i,j]==(7*pi/4))
    
  {
    runon[i+1,j+1]<-qf #bottom tight
  }
  
  if(flowdir[i,j]>0 && flowdir[i,j]<(pi/4))
  {
    if(flowdir[i,j]/((pi/2)-(pi/4))<0.5) {
      runon[i,j+1]<-qf*(1-flowdir[i,j]/((pi/4)-0)) ## more to the right
      runon[i-1,j+1] <- qf*((flowdir[i,j]/((pi/4)-0)))
    }
    else{
     runon[i-1,j+1]<-qf*(1-(flowdir[i,j]/((pi/4)-0))) ##more to the top right
     runon[i,j+1]<- qf*((flowdir[i,j]/((pi/4)-0)))
    }}

    
     if(flowdir[i,j]>(pi/4) && flowdir[i,j]<(pi/2))
       {
         if((flowdir[i,j]-(pi/4))/((pi/4)-0)<0.5) {
        runon[i-1,j+1]<-qf*((1-flowdir[i,j]-(pi/4))/((pi/4)-0)) ## more to the top right
        runon[i-1,j]<-qf*((flowdir[i,j]-(pi/4))/((pi/4)-0))
         }
         else{
           runon[i-1,j]<-qf*((1-(flowdir[i,j]-(pi/4))/((pi/4)-0))) ##more to the top 
           runon[i-1,j+1]<-qf*((flowdir[i,j]-(pi/4))/((pi/4)-0))
            }}
      
        
          if(flowdir[i,j]>(pi/2) && flowdir[i,j]<(3*pi/4)){
         
           if((flowdir[i,j]-(pi/2))/((pi/4)-0)<0.5) {
              runon[i-1,j]<-qf*((1-(flowdir[i,j]-(pi/2)))/((pi/4)-0)) ## more to the top 
              runon[i-1,j-1]<-qf*((flowdir[i,j]-(pi/2))/((pi/4)-0))
           }
           else{
              runon[i-1,j-1]<-qf*((1-(flowdir[i,j]-(pi/2))/((pi/4)-0))) ##more to the top left
              runon[i-1,j]<-qf*((flowdir[i,j]-(pi/2))/((pi/4)-0))
              }}
          
             if(flowdir[i,j]>(3*pi/4) && flowdir[i,j]<pi)
               
              {
                 if((flowdir[i,j]-(3*pi/4))/((pi/4)-0)<0.5) {
                   runon[i-1,j-1]<-((1-(flowdir[i,j]-(3*pi/4)))/((pi/4)-0)) ## more to the top left
                   runon[i,j-1]<-qf*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                 }
                else{
                    runon[i,j-1]<-qf*((1-(flowdir[i,j]-(3*pi/4))/((pi/4)-0))) ##more to the left
                    runon[i-1,j-1]<-qf*((flowdir[i,j]-(3*pi/4))/((pi/4)-0))
                  }}
             
                if(flowdir[i,j]>pi && flowdir[i,j]<(5*pi/4))
                {
                  if((flowdir[i,j]-(pi))/((pi/4)-0)<0.5) {
                    runon[i,j-1]<-qf*((1-(flowdir[i,j]-(pi)))/((pi/4)-0)) ## more to the left
                    runon[i+1,j-1]<-qf*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }
                  else{
                    runon[i+1,j-1]<-qf*((1-(flowdir[i,j]-(pi))/((pi/4)-0))) ##more to the bottom left
                    runon[i,j-1]<-qf*((flowdir[i,j]-(pi))/((pi/4)-0))
                  }}
                  
                  if(flowdir[i,j]>(5*pi/4) && flowdir[i,j]<(3*pi/2))
                  {
                    if((flowdir[i,j]-(5*pi/4))/((pi/4)-0)<0.5) {
                      runon[i+1,j-1]<-qf*((1-(flowdir[i,j]-(5*pi/4)))/((pi/4)-0)) ## more to the bottom left
                      runon[i+1,j]<-qf*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                    else{
                        runon[i+1,j]<-qf*((1-(flowdir[i,j]-(5*pi/4))/((pi/4)-0))) ##more to the bottom
                        runon[i+1,j-1]<-qf*((flowdir[i,j]-(5*pi/4))/((pi/4)-0))
                    }
                  }
    
                        
                    if(flowdir[i,j]>(3*pi/2) && flowdir[i,j]<(7*pi/4))
                    {
                      if((flowdir[i,j]-(3*pi/2))/((pi/4)-0)<0.5) {
                         runon[i+1,j]<-qf*((1-(flowdir[i,j]-(3*pi/2)))/((pi/4)-0)) ## more to the bottom 
                         runon[i+1,j+1]<-qf*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                      else{
                         runon[i+1,j+1]<-qf*((1-(flowdir[i,j]-(3*pi/2))/((pi/4)-0))) ##more to the bottom right
                         runon[i+1,j]<-qf*((flowdir[i,j]-(3*pi/2))/((pi/4)-0))
                      }
                    }
                      
                      
                      if(flowdir[i,j]>(7*pi/4) && flowdir[i,j]<(2*pi)){
                      
                        if((flowdir[i,j]-(7*pi/4))/((pi/4)-0)<0.5) {
                          runon[i+1,j+1]<-qf*((1-(flowdir[i,j]-(7*pi/4)))/((pi/4)-0)) ## more to the bottom right
                          runon[i,j+1]<-qf*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        else{
                           runon[i,j+1]<-qf*((1-(flowdir[i,j]-(7*pi/4))/((pi/4)-0))) ##more to the right
                           runon[i+1,j+1]<-qf*((flowdir[i,j]-(7*pi/4))/((pi/4)-0))
                        }
                        
                     
                      }
      }
    }

return(runon)


        
                   }

runon_fun(flowdir=flowdir,qf=blabla)
