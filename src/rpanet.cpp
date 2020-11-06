#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List rpanet_cpp(int       nsteps,
	                   arma::vec  control,
	                   const bool directed, 
	                   arma::vec  m,
	                   arma::vec  w,
	                   arma::vec  outstrength,
	                   arma::vec  instrength, 
	                   int        nnode,
	                   int        tnode) {
  double alpha = control[0], beta   = control[1], gamma = control[2];
  double u;
  int v1, v2, count = 0;
  arma::vec nodes = arma::linspace<arma::vec>(1, outstrength.size(), outstrength.size());
  arma::vec startnode(outstrength.size(), arma::fill::zeros); 
  arma::vec endnode(outstrength.size(), arma::fill::zeros);
  arma::vec subnode, subinstrength, suboutstrength;
  

  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < m[i]; j++) {
      tnode--;
      subnode = nodes.subvec(0, tnode);
      subinstrength = instrength.subvec(0, tnode);
      suboutstrength = outstrength.subvec(0, tnode);
      u = R::runif(0, 1);
      if (u < alpha) {
      	v1 = nnode + 1;
      	v2 = Rcpp::RcppArmadillo::sample(subnode, 1L, true, subinstrength)[0];
      	nnode++;
      } else if (u < alpha + beta) {
      	v1 = Rcpp::RcppArmadillo::sample(subnode, 1L, true, suboutstrength)[0];
      	v2 = Rcpp::RcppArmadillo::sample(subnode, 1L, true, subinstrength)[0];
      } else if (u < alpha + beta + gamma) {
      	v1 = Rcpp::RcppArmadillo::sample(subnode, 1L, true, suboutstrength)[0];
      	v2 = nnode + 1;
      	nnode++;
      } else {
      	v1 = nnode + 1;
      	v2 = nnode + 2;
      	nnode = nnode + 2;
      }
      startnode[count] = v1;
      endnode[count] = v2;
      outstrength[v1 - 1] += w[count];
      instrength[v2 - 1] += w[count];
      if (! directed) {
        outstrength[v2 - 1] += w[count];
        instrength[v1 - 1] += w[count];
      }
      count++;
      tnode = nnode;
    }
  }
  
  Rcpp::List ret;
  ret["startnode"] = startnode.subvec(0, count - 1);
  ret["endnode"] = endnode.subvec(0, count - 1);
  ret["instrength"] = instrength.subvec(0, tnode - 1);
  ret["outstrength"] = outstrength.subvec(0, tnode - 1);
  return ret;
}
  
  
			       
