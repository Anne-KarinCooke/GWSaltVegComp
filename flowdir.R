
### preparations, put them somewhere else later
flowdir <- ang # ang is angle flowdir[i,j] from DInf TauDEM as raster
flowdir[is.na(flowdir)] <- 8   ###********************The loop had problems with NA, so I changed NA from the blundaries to be 8. 8 is outside of 2pi...

###
runon_fun<-function(flowdir,q){
##### For Loop to fill runon with runoff from the right cell:
  
for (i in 1:5){
  
      for (j in 1:5) {
        
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
                        
                     
                      }
      }
    }

return(runon)


        
                   }



