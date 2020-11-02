//#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat
rpanet_cpp(int         nsteps,
	   const arma::vec   control,
	   const bool        directed = true)

{
  double alpha = control[0], beta   = control[1], gamma = control[2];
  double delin = control[3], delout = control[4];
  int batsize = 1000;
  
  int m = 1, nnode = 1, v1, v2; // change later
  double u;
  arma::mat outmat(nsteps, 4, arma::fill::zeros);
  arma::vec ins(nsteps, delin);
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < m; j++) {
      u = R::runif(0, 1);
      if (u < alpha) {
	v1 = nnode + 1;
	v2 = Rcpp::RcppArmadillo::sample(nnode, 1, true, ins.subvec(0, nnode - 1));
	nnode++;
      } else if (u < alpha + beta) {
	v1 = Rcpp::RcppArmadillo::sample(nnode, 1, true, ins.subvec(0, nnode - 1));
	v2 = Rcpp::RcppArmadillo::sample(nnode, 1, true, ins.subvec(0, nnode - 1));
      } else if (u < alpha + beta + gamma) {
	v1 = Rcpp::RcppArmadillo::sample(nnode, 1, true, ins.subvec(0, nnode - 1));
	v2 = nnode + 1;
	nnode++;
      } else {
	v1 = nnode + 1;
	v2 = nnode + 2;
	nnode = nnode + 2;
      }
      // continue with the rest
      
    }
    
  }
  return outmat;
}
  
  
			       
