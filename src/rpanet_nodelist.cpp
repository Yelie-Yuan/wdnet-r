#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Preferential attachment algorithm for simple situations, 
//' i.e., edge weight equals to 1, number of new edges per step is 1.
//'
//' @param snode Source nodes.
//' @param tnode Target nodes.
//' @param scenario Sequence of alpha, beta, gamma, xi, rho scenarios.
//' @param nnode Number of nodes in seed network.
//' @param nedge Number of edges in seed network.
//' @param delta_out Tuning parameter.
//' @param delta_in Tuning parameter.
//' @param directed Whether the network is directed.
//' @return Number of nodes, sequences of source and target nodes.
//' 
// [[Rcpp::export]]
Rcpp::List rpanet_nodelist_cpp(arma::vec snode,
                               arma::vec tnode,
                               arma::vec scenario,
                               int nnode,
                               int nedge,
                               double delta_out,
                               double delta_in, 
                               bool directed) {
  GetRNGstate();
  int n = scenario.size();
  double u, v;
  int j;
  for (int i = 0; i < n; i++) {
    j = scenario[i];
    switch(j) {
      case 1: {
        u = unif_rand() * (nedge + nnode * delta_in);
        if (u < nedge) {
          if (directed) {
            tnode[nedge] = tnode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              tnode[nedge] = snode[floor(u)];
            } 
            else {
              tnode[nedge] = tnode[floor(u)];
            }
          }
        }
        else {
          tnode[nedge] = ceil((u - nedge) / delta_in);
        }
        nnode++;
        snode[nedge] = nnode;
        break;
      }
      case 2: {
        u = unif_rand() * (nedge + nnode * delta_out);
        if (u < nedge) {
          if (directed) {
            snode[nedge] = snode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              snode[nedge] = snode[floor(u)];
            } 
            else {
              snode[nedge] = tnode[floor(u)];
            }
          }
        } 
        else {
          snode[nedge] = ceil((u - nedge) / delta_out);
        }
        
        u = unif_rand() * (nedge + nnode * delta_in);
        if (u < nedge) {
          if (directed) {
            tnode[nedge] = tnode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              tnode[nedge] = snode[floor(u)];
            } 
            else {
              tnode[nedge] = tnode[floor(u)];
            }
          }
        } 
        else {
          tnode[nedge] = ceil((u - nedge) / delta_in);
        }
        break;
      }
      case 3: {
        u = unif_rand() * (nedge + nnode * delta_out);
        if (u < nedge) {
          if (directed) {
            snode[nedge] = snode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              snode[nedge] = snode[floor(u)];
            } 
            else {
              snode[nedge] = tnode[floor(u)];
            }
          }
        } 
        else {
          snode[nedge] = ceil((u - nedge) / delta_out);
        }
        nnode++;
        tnode[nedge] = nnode;
        break;
      }
      case 4: {
        nnode += 2;
        snode[nedge] = nnode - 1;
        tnode[nedge] = nnode;
        break;
      }
      case 5: {
        nnode++;
        snode[nedge] = nnode;
        tnode[nedge] = nnode;
        break;
      }
    }
    nedge++;
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["snode"] = snode;
  ret["tnode"] = tnode;
  ret["nnode"] = nnode;
  return ret;
}
