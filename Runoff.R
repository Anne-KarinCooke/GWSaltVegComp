OF<- function(h, soilpar,slope){
  
  qq=(soilpar$cn/soilpar$Mn)*h^(5/3)*slope^(1/2)
  
  return(qq)
}

