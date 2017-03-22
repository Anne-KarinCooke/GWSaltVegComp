
// VEGETATION FUNCTIONS

//Plant water uptake function WU

double WU(double M, double P, double gmax, double k1 ) {  /// not quite happy with the list item calling...list par
  
  double Wu = gmax*(M/(M+k1))*P;   // function is called WU, output is called Wu, cannot be the same
  return Wu;
  
}


//Plant Growth function Gr

double Gr(double M, double P, double c, double gmax, double k1){
  
  double Gro = c*WU(M,P,gmax,k1);
  return Gro;
}


// Plant mortality function Mo

double Mo(double P, double M, double Svir, double d ){
  
  double Mort=P*(d*(M/Svir));
  return Mort;
}  
