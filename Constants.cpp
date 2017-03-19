#include <Rcpp.h>
using namespace Rcpp;

//Sandy Clay Loam
      
      
      List Constants_cpp(List soilpar,
                         List saltpar,
                         List vegpar) 
      {
        
        
        // SOIL 
        
        
        const double n = soilpar["n"];
        const double K_s = soilpar["K_s"];
        const double b = soilpar["b"];
        const double s_fc = soilpar["s_fc"];
        const double psi_s_bar = soilpar["psi_s_bar"];
        const double h1bar = soilpar["h1bar"];
        const double hb = soilpar["hb"];
        const double Mn = soilpar["Mn"];
        const double cn = soilpar["cn"];
        double alpha_i = soilpar["alpha_i"]; 
        
        soilpar["n"] = 0.367; //porosity
        soilpar["K_s"] = 52.08*10; //Hydraulic conductivity mm/day 
        soilpar["b"] = 6.4069; // campbell's b
        soilpar["s_fc"] = 0.2677/n; // Field capacity
        soilpar["psi_s_bar"] = -1.2E-3; // bubbling pressure
        soilpar["h1bar"] =  -psi_s_bar;
        soilpar["hb"] = psi_s_bar*pow(-10,5);//  mm
        soilpar["Mn"] = 10; // Manning's n
        soilpar["cn"] = 1; // runoff conversion factor
        soilpar["alpha_i"] = 1;//#maximum infiltration rate per day, This needs to be a fraction of h (p117 Saco and Moreno-Las Heras) 

        
        
        
        // VEGETATION
                    const double Zr = vegpar["Zr"];
                    const double k = vegpar["k"];
                    const double W0 = vegpar["W0"];
                    const double gmax = vegpar["gmax"];
                    const double c = vegpar["c"];
                    const double k1 = vegpar["k1"];
                    const double d = vegpar["d"]; //fraction of plant mortality
                    
                    vegpar["Zr"] = 400; //mm, Grass
                    vegpar["d"] = 0.24;//Saco et al, 2013
                    vegpar["c"] = 0.10;//Saco et al, 2013
                    vegpar["k1"] = 5;//Saco et al, 2013
                    vegpar["gmax"] = 0.05;//Saco et al, 2013
                    vegpar["W0"] = 0.2;//Saco et al, 2013
                    vegpar["k"] = 12;//Saco et al, 2013
                    
                    
                    
        
        // SALT
                                    const double ConcConst = saltpar["ConcConst"];
                                    const double CMgw = saltpar["CMgw"];
                                    const double f = saltpar["f"];
                                  
                                    saltpar["ConcConst"] = 0;  //ConcConst is the concentration of the salt in the infiltrating water in g/l
                                    saltpar["CMgw"] = 0; //CMgw is the goundwater salt concentration  in g/l
                                    saltpar["f"] = 0.8; //f is the soil salt leaching efficiency (whether some salt is retained)
      }
      

          
          