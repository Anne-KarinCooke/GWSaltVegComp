
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

int rows =10;
int cols=10;

  // [[Rcpp::export]]
  NumericMatrix call(NumericMatrix xfel, CharacterVector format = "GTiff", CharacterVector filename = "xfel.tif",
                     CharacterVector command = "mpiexec -n 8 DinfFlowdir -ang elevang.tif -slp elevslp.tif -fel elevfel.tif", int rows =10, int cols=10){
     
     Environment rgdal("package:rgdal");
     Environment Matrix("package:Matrix");
     Environment raster("package:raster");
     Function f1 = raster["raster"];
     Function f2 = rgdal["writeRaster"];
     Function f3 = Matrix["as.matrix"];

    // Environment env = Environment::global_env();
     f2(f1(xfel));
   
    
    system("mpiexec -n 8 DinfFlowdir -ang xang.tif -slp xslp.tif -fel xfel.tif");
    
    NumericMatrix slp(rows,cols);
    slp = f3(f1("xslp.tif"));
    NumericMatrix ang(rows,cols);
    ang = f3(f1("xang.tif"));
    
    return slp;
    return ang;

  }

  NumericMatrix xfel(rows, cols, runif(100));
  NumericMatrix result = call(xfel);
  /*** R
   xfel <- matrix(runif(100),nrow= 10, ncol =10)
  rows =10
  cols=10
   call(xfel) #as.matrix()
  */
