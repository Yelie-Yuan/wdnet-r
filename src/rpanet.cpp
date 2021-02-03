#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
int my_sample(int tnode, double s_strength, arma::vec strength, double delta) {
  int i;
  double j = 0, rT;
  rT = unif_rand() * (s_strength + tnode * delta);
  for (i = 0; i < tnode; i++) {
    j = j + strength[i] + delta;
    if (j >= rT)
      break;
  }
  return i + 1;
}

// [[Rcpp::export]]
Rcpp::List rpanet_cpp(int       nsteps,
                      arma::vec  control,
                      const bool directed, 
                      arma::vec  m,
                      arma::vec  w,
                      arma::vec  u,
                      arma::vec  outstrength,
                      arma::vec  instrength, 
                      double     s_outstrength,
                      double     s_instrength,
                      int        nnode,
                      int        tnode) {
  double alpha = control[0], beta = control[1], gamma = control[2], xi = control[3];
  double delta_out = control[4], delta_in = control[5];
  double u2;
  int v1, v2, count = 0;
  arma::vec startnode(outstrength.size(), arma::fill::zeros); 
  arma::vec endnode(outstrength.size(), arma::fill::zeros);
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < m[i]; j++) {
      u2 = u[count];
      if (u2 < alpha) {
        v1 = nnode + 1;
        v2 = my_sample(tnode, s_instrength, instrength, delta_in);
        nnode++;
      } else if (u2 < alpha + beta) {
        v1 = my_sample(tnode, s_outstrength, outstrength, delta_out);
        v2 = my_sample(tnode, s_instrength, instrength, delta_in);
      } else if (u2 < alpha + beta + gamma) {
        v1 = my_sample(tnode, s_outstrength, outstrength, delta_out);
        v2 = nnode + 1;
        nnode++;
      } else if (u2 < alpha + beta + gamma + xi) {
        v1 = nnode + 1;
        v2 = nnode + 2;
        nnode += 2;
      } else {
        v1 = nnode + 1;
        v2 = nnode + 1;
        nnode++;
      }
      startnode[count] = v1;
      endnode[count] = v2;
      outstrength[v1 - 1] += w[count];
      instrength[v2 - 1] += w[count];
      s_outstrength += w[count];
      s_instrength += w[count];
      if (! directed) {
        outstrength[v2 - 1] += w[count];
        instrength[v1 - 1] += w[count];
        s_outstrength += w[count];
        s_instrength += w[count];
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



