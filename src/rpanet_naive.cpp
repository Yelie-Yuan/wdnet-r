#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Sample a node according to node strength.
//'
//' @param tnode Number of nodes at current step.
//' @param sumstrength Total strength at current step.
//' @param strength Vector of node strength.
//' @param delta The tuning parameter, delta_in or delta_out.
//' @return Sampled node.
// [[Rcpp::export]]
int sampleNode_naive_cpp(int tnode, double sumstrength, 
               arma::vec strength, double delta) {
  int i;
  double j = 0, v;
  v = unif_rand() * (sumstrength + tnode * delta);
  for (i = 0; i < tnode; i++) {
    j = j + strength[i] + delta;
    if (j >= v)
      break;
  }
  return i;
}


//' Sample a node according to node strength.
//'
//' @param nstep Number of steps.
//' @param control Vector of control parameters, i.e., alpha, beta,
//'   gamma, xi, delta_in and delta_out.
//' @param m Vector, number of edges at each step.
//' @param w Weight of new edges.
//' @param outstrength Vector of node out-strength.
//' @param instrength Vector of node in-strength.
//' @param sumstrength Total strength of inital network.
//' @param nnode Number of nodes of inital network.
//' @return A list of source nodes, target nodes, node out- and in-strength.
// [[Rcpp::export]]
Rcpp::List rpanet_naive_cpp(int        nstep,
                             arma::vec  control,
                             arma::vec  m,
                             arma::vec  w,
                             arma::vec  outstrength,
                             arma::vec  instrength, 
                             double     sumstrength,
                             int        nnode) {
  GetRNGstate();
  double alpha = control[0], beta = control[1], gamma = control[2], xi = control[3];
  double delta_out = control[4], delta_in = control[5];
  double u;
  int count = 0, tnode = nnode, max_m = arma::max(m);
  int i, j;
  arma::vec v1(max_m, arma::fill::zeros);
  arma::vec v2(max_m, arma::fill::zeros);
  arma::vec scenario(max_m, arma::fill::zeros);
  arma::vec start_node(outstrength.size(), arma::fill::zeros); 
  arma::vec end_node(outstrength.size(), arma::fill::zeros);
  arma::vec edgescenario(arma::accu(m), arma::fill::zeros); 
  for (i = 0; i < nstep; i++) {
    for (j = 0; j < m[i]; j++) {
      u = unif_rand();
      if (u <= alpha) {
        scenario[j] = 1;
        v1[j] = nnode;
        nnode++;
        v2[j] = sampleNode_naive_cpp(tnode, sumstrength, instrength, delta_in);
      } else if (u <= alpha + beta) {
        scenario[j] = 2;
        v1[j] = sampleNode_naive_cpp(tnode, sumstrength, outstrength, delta_out);
        v2[j] = sampleNode_naive_cpp(tnode, sumstrength, instrength, delta_in);
      } else if (u <= alpha + beta + gamma) {
        scenario[j] = 3;
        v1[j] = sampleNode_naive_cpp(tnode, sumstrength, outstrength, delta_out);
        v2[j] = nnode;
        nnode++;
      } else if (u <= alpha + beta + gamma + xi) {
        scenario[j] = 4;
        v1[j] = nnode;
        v2[j] = nnode + 1;
        nnode += 2;
      } else {
        scenario[j] = 5;
        v1[j] = nnode;
        v2[j] = nnode;
        nnode++;
      }
    }
    for (j = 0; j < m[i]; j++) {
      edgescenario[count] = scenario[j];
      start_node[count] = v1[j];
      end_node[count] = v2[j];
      outstrength[v1[j]] += w[count];
      instrength[v2[j]] += w[count];
      sumstrength += w[count];
      count++;
    }
    tnode = nnode;
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["start_node"] = start_node.subvec(0, count - 1);
  ret["end_node"] = end_node.subvec(0, count - 1);
  ret["instrength"] = instrength.subvec(0, nnode - 1);
  ret["outstrength"] = outstrength.subvec(0, nnode - 1);
  ret["scenario"] = edgescenario;
  return ret;
}