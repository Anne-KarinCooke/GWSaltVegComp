#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

mat Surface(int rows, int cols, mat flowdir, mat flowdirTable, mat q){
  //flowdirTable

//const double pi = 3.141593; 

mat runon(rows, cols,fill::zeros); 

//double number = (flowdir(i,j)/(pi/4));
      
int a;
//int x;
//int y;
int i;
int j;

  for (i=0; i< rows; i++) {
 // 
    for (j=0; j< cols; j++ ){    
     
        for (a=0; a < 9; a++) {
          
      if (flowdir(i,j) == flowdirTable(0,a)){   ///(pi/4)

     // x = flowdirTable(a,1);
     // y = flowdirTable(a,2);
      
      runon(i+flowdirTable(1,a),j+flowdirTable(2,a)) += q(i,j);
        
        }
       }
    }
 }
 return runon;
 
}
 
  
  
//   
// mat SurfRedistr( NumericMatrix flowdir, mat filler, mat origin, int i, int j,int rows, int cols, const double pi = 3.141593){
//   
//   mat destination(rows, cols,fill::zeros); 
//   
//   
//   int ii;
//   int jj;
//   
//   for (ii=1; ii< (rows-1); ii++) {
//     
//     for (jj=1; jj< (cols-1); jj++ ){
//       
//       // if(((ii+1) <= cols) &  ((jj+1) <= rows) & ((ii-1) >= 1) & ((jj-1) >= 1)){
//       
//       
//       
//       if(flowdir(ii,jj) == 0)
//       {
//         destination(ii,jj+1) += filler(ii,jj) * origin(ii,jj); // right
//       }
//       if(flowdir(ii,jj) == (pi/4))
//       {
//         destination(ii-1,jj+1) += filler(ii,jj) * origin(ii,jj); // top right
//       }
//       if(flowdir(ii,jj) == (pi/2))
//       {
//         destination(ii-1,jj) += filler(ii,jj) * origin(ii,jj); // top 
//       }
//       if(flowdir(ii,jj) == (3*pi/4))
//       {
//         destination(ii-1,jj-1) += filler(ii,jj) * origin(ii,jj); // top left
//       }
//       if(flowdir(ii,jj) == pi)
//       {
//         destination(ii,jj-1) += filler(ii,jj) * origin(ii,jj); //  left
//       }
//       if(flowdir(ii,jj) == (5*pi/4))
//       {
//         destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj); //  bottom left
//       }
//       if(flowdir(ii,jj) == (3*pi/2))
//       {
//         destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj); //  bottom 
//       }
//       if(flowdir(ii,jj) == (7*pi/4))
//       {
//         destination(ii+1,jj+1) += filler(ii,jj) * origin(ii,jj); //  bottom right
//       }
//       
//       
//       if((flowdir(ii,jj)>0) & (flowdir(ii,jj)<(pi/4)))
//       {
//         if(flowdir(ii,jj)/((pi/2)-(pi/4))<0.5) {
//           destination(ii,jj+1) += filler(ii,jj) * origin(ii,jj)*(1-flowdir(ii,jj)/((pi/4)-0)); //// more to the right
//           destination(ii-1,jj+1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)/((pi/4)-0)));
//         }
//         else{
//           destination(ii-1,jj+1) += filler(ii,jj) * origin(ii,jj)*(1-(flowdir(ii,jj)/((pi/4)-0))); ////more to the top right
//           destination(ii,jj+1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)/((pi/4)-0)));
//         }}
//       
//       
//       if((flowdir(ii,jj)>(pi/4)) & (flowdir(ii,jj)<(pi/2)))
//       {
//         if((flowdir(ii,jj)-(pi/4))/((pi/4)-0)<0.5) {
//           destination(ii-1,jj+1) += filler(ii,jj) * origin(ii,jj)*((1-flowdir(ii,jj)-(pi/4))/((pi/4)-0)); //// more to the top right
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi/4))/((pi/4)-0));
//         }
//         else{
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(pi/4))/((pi/4)-0))); ////more to the top 
//           destination(ii-1,jj+1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi/4))/((pi/4)-0));
//         }}
//       
//       
//       if((flowdir(ii,jj)>(pi/2)) & (flowdir(ii,jj)<(3*pi/4))){
//         
//         if((flowdir(ii,jj)-(pi/2))/((pi/4)-0)<0.5) {
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(pi/2)))/((pi/4)-0)); //// more to the top 
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi/2))/((pi/4)-0));
//         }
//         else{
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(pi/2))/((pi/4)-0))); ////more to the top left
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi/2))/((pi/4)-0));
//         }}
//       
//       if((flowdir(ii,jj)>(3*pi/4)) & (flowdir(ii,jj)<pi))
//         
//       {
//         if((flowdir(ii,jj)-(3*pi/4))/((pi/4)-0)<0.5) {
//           destination(ii+1,jj-1) += ((1-(flowdir(ii,jj)-(3*pi/4)))/((pi/4)-0)); //// more to the top left
//           destination(ii,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(3*pi/4))/((pi/4)-0));
//         }
//         else{
//           destination(ii,jj-1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(3*pi/4))/((pi/4)-0))); ////more to the left
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(3*pi/4))/((pi/4)-0));
//         }}
//       
//       if((flowdir(ii,jj)>pi) & (flowdir(ii,jj)<(5*pi/4)))
//       {
//         if((flowdir(ii,jj)-(pi))/((pi/4)-0)<0.5) {
//           destination(ii,jj-1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(pi)))/((pi/4)-0)); //// more to the left
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi))/((pi/4)-0));
//         }
//         else{
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(pi))/((pi/4)-0))); ////more to the bottom left
//           destination(ii,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(pi))/((pi/4)-0));
//         }}
//       
//       if((flowdir(ii,jj)>(5*pi/4)) & (flowdir(ii,jj)<(3*pi/2)))
//       {
//         if((flowdir(ii,jj)-(5*pi/4))/((pi/4)-0)<0.5) {
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(5*pi/4)))/((pi/4)-0)); //// more to the bottom left
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(5*pi/4))/((pi/4)-0));
//         }
//         else{
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(5*pi/4))/((pi/4)-0))); ////more to the bottom
//           destination(ii+1,jj-1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(5*pi/4))/((pi/4)-0));
//         }
//       }
//       
//       
//       if((flowdir(ii,jj)>(3*pi/2)) & (flowdir(ii,jj)<(7*pi/4)))
//       {
//         if((flowdir(ii,jj)-(3*pi/2))/((pi/4)-0)<0.5) {
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(3*pi/2)))/((pi/4)-0)); //// more to the bottom 
//           destination(ii+1,jj+1) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(3*pi/2))/((pi/4)-0));
//         }
//         else{
//           destination(ii+1,jj+1) += filler(ii,jj) * origin(ii,jj)*((1-(flowdir(ii,jj)-(3*pi/2))/((pi/4)-0))); ////more to the bottom right
//           destination(ii+1,jj) += filler(ii,jj) * origin(ii,jj)*((flowdir(ii,jj)-(3*pi/2))/((pi/4)-0));
//         }
//       }
//       
//       
//       if((flowdir(ii,jj)>(7*pi/4)) & (flowdir(ii,jj)<(2*pi))){
//         
//         if((flowdir(ii,jj)-(7*pi/4))/((pi/4)-0)<0.5) {
//           destination(ii+1,jj+1) += origin(ii,jj)*((1-(flowdir(ii,jj)-(7*pi/4)))/((pi/4)-0)); //// more to the bottom right
//           destination(ii,jj+1) += origin(ii,jj)*((flowdir(ii,jj)-(7*pi/4))/((pi/4)-0));
//         }
//         else{
//           destination(ii,jj+1) += origin(ii,jj)*((1-(flowdir(ii,jj)-(7*pi/4))/((pi/4)-0))); ////more to the right
//           destination(ii+1,jj+1) += origin(ii,jj)*((flowdir(ii,jj)-(7*pi/4))/((pi/4)-0));
//         }
//         
//       }
//     }
//   }
//   
//   // }
//   
//   return  destination;
// }


/*** R

# flowdir <- read.table("flowdir.txt")
# flowdir <- as.matrix(flowdir)
# 
# flowdir <- apply(flowdir, 1, as.numeric)
# flowdir
# flowdirTable <-read.table("flowdirTable.txt")
# flowdirTable <- apply(flowdirTable, 1, as.numeric)
# flowdirTable
# qq
# whatever <- Surface(2,3,10,10,flowdir,flowdirTable,qq)
# whatever
# whatever <- matrix(0,nrow= 10, ncol =10)
# for (i in 2:8) { 
#   
#   for (j in 2:8){
#   whatever[i,j] <- Surface(i,j,10,10,flowdir,flowdirTable,qq)
#   }
# }
whatever <- Surface(10,10,flowdir,flowdirTable,qq)

 
*/
