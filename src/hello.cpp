#include <RcppArmadillo.h>   

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void hello_world() {
  Rcpp::Rcout << "Hello World!" << std::endl;  
}



// [[Rcpp::export]]
double fx(arma::colvec x, arma::mat Y,
          arma::colvec z) {
    // calculate the result
  double result = arma::as_scalar(arma::trans(x) * arma::inv(Y) * z);
  return result;
}


// After compile, this function will be immediately called using
// the below snippet and results will be sent to the R console.

/*** R
hello_world() 
*/
