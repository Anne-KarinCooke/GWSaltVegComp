#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Salt_cpp(std::string stype) {

     // default, NO SALT from nowhere
     double ConcConst = 0.0; //ConcConst is the concentration of the salt in the infiltrating water in g/l
     double CMgw = 0.0; //CMgw is the goundwater salt concentration  in g/l
     double f = 1; //f is the soil salt leaching efficiency (whether some salt is retained)
  
  
  if (stype == "Groundwater") {
    
    ConcConst = 0.0;
    CMgw = 0.1;
    f = 1;
    
  }
  
  if (stype == "Rain") {
    
    ConcConst = 0.1;
    CMgw = 0.0;
    f = 1;
    
  }
  
  if (stype == "Both") {
    
    ConcConst = 0.1;
    CMgw = 0.1;
    f = 1;
    
  }
  
  if (stype == "None") {
    
    ConcConst = 0.0;
    CMgw = 0.0;
    f = 1;
    
  }
  
  
  return(Rcpp::List::create(Rcpp::Named("ConcCOnst") = ConcConst,
                            Rcpp::Named("CMgw") = CMgw,
                            Rcpp::Named("f") = f));
                            
}
