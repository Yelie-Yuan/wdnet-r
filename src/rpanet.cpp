#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Fill missing values in node sequence.
//'
//' @param nodes Sequence of source/target node of edges, missing values are denoted as 0.
//' @param edges Sampled edges according to preferential attachment.
//' @return Sequence of source/target node of edges.
// [[Rcpp::export]]
arma::vec findNode_cpp(arma::vec nodes, 
                       arma::vec edges) {
  int n = nodes.size(), n1 = 0;
  for (int j = 0; j < n; j++) {
    if (nodes[j] == 0) {
      nodes[j] = nodes[edges[n1] - 1];
      n1++;
    }
  }
  return nodes;
}

//' Fill missing values in node sequence.
//'
//' @param start_node Sequence of nodes from the first column of edgelist, i.e., \code{edgelist[, 1]}.
//' @param end_node Sequence of nodes from the first column of edgelist, i.e., \code{edgelist[, 2]}.
//' @param start_edge Index of sampled edges, corresponds to the missing nodes in start_node.
//' @param end_edge Index of sampled edges, corresponds to the missing nodes in end_node.
//' @return Sequence of source/target node of edges.
// [[Rcpp::export]]
Rcpp::List findNode_undirected_cpp(arma::vec start_node, 
                       arma::vec end_node, 
                       arma::vec start_edge, 
                       arma::vec end_edge) {
  GetRNGstate();
  int n = start_node.size(), n1 = 0, n2 = 0;
  double u;
  for (int j = 0; j < n; j++) {
    if (start_node[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        start_node[j] = start_node[start_edge[n1]  - 1];
      } else {
        start_node[j] = end_node[start_edge[n1] - 1];
      }
      n1++;
    }
    if (end_node[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        end_node[j] = start_node[end_edge[n2] - 1];
      } else {
        end_node[j] = end_node[end_edge[n2] - 1];
      }
      n2++;
    }
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["start_node"] = start_node;
  ret["end_node"] = end_node;
  return ret;
}

//' Aggregate edgeweight into nodes' strength.
//'
//' @param start_node Sequence of source nodes.
//' @param end_node Sequence of target nodes.
//' @param weight Sequence of edgeweight.
//' @param nnode Number of nodes of sampled network.
//' @param weighted Logical, ture if the edges are weighted, 
//'   false if not.
//' @return Sequence of outstrength/instrength.
// [[Rcpp::export]]
Rcpp::List nodeStrength_cpp(arma::vec start_node, 
                            arma::vec end_node,
                            arma::vec weight, 
                            int nnode, 
                            bool weighted = true) {
  int n = start_node.size();
  arma::vec outstrength(nnode, arma::fill::zeros);
  arma::vec instrength(nnode, arma::fill::zeros);
  if (weighted) {
    for (int i = 0; i < n; i++) {
      outstrength[start_node[i] - 1] += weight[i];
      instrength[end_node[i] - 1] += weight[i];
    }
  } else {
    for (int i = 0; i < n; i++) {
      outstrength[start_node[i] - 1] += 1;
      instrength[end_node[i] - 1] += 1;
    }
  }
  
  Rcpp::List ret;
  ret["outstrength"] = outstrength;
  ret["instrength"] = instrength;
  return ret;
}

//' Sample nodes with respect to the number of nodes at each step.
//'
//' @param total_node Number of nodes at each step.
//' @return Sequence of sampled nodes.
// [[Rcpp::export]]
arma::vec sampleNode_cpp(arma::vec total_node) {
  GetRNGstate();
  int n = total_node.size();
  arma::vec nodes(n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    nodes[i] = Rcpp::sample(total_node[i], 1)[0];
  }
  PutRNGstate();
  return nodes;
}


//' Preferential attachment algorithm for simple situations, 
//' e.g., edge weights are constant, 
//' number of new edges per step is 1.
//'
//' @param start_node Sequence of source nodes.
//' @param end_node Sequence of target nodes.
//' @param scenario Sequence of alpha, beta, gamma, xi, rho scenarios.
//' @param nnode Number of nodes at current step.
//' @param nedge Number of edges at current step.
//' @param delta_out Tuning parameter.
//' @param delta_in Tuning parameter.
//' @param directed Whether the network is directed.
//' @return Number of nodes, sequences of source and target nodes.
//' 
// [[Rcpp::export]]
Rcpp::List rpanet_cpp(arma::vec start_node,
                      arma::vec end_node,
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
            end_node[nedge] = end_node[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              end_node[nedge] = start_node[floor(u)];
            } 
            else {
              end_node[nedge] = end_node[floor(u)];
            }
          }
        }
        else {
          end_node[nedge] = ceil((u - nedge) / delta_in);
        }
        nnode++;
        start_node[nedge] = nnode;
        break;
      }
      case 2: {
        u = unif_rand() * (nedge + nnode * delta_out);
        if (u < nedge) {
          if (directed) {
            start_node[nedge] = start_node[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              start_node[nedge] = start_node[floor(u)];
            } 
            else {
              start_node[nedge] = end_node[floor(u)];
            }
          }
        } 
        else {
          start_node[nedge] = ceil((u - nedge) / delta_out);
        }
        
        u = unif_rand() * (nedge + nnode * delta_in);
        if (u < nedge) {
          if (directed) {
            end_node[nedge] = end_node[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              end_node[nedge] = start_node[floor(u)];
            } 
            else {
              end_node[nedge] = end_node[floor(u)];
            }
          }
        } 
        else {
          end_node[nedge] = ceil((u - nedge) / delta_in);
        }
        break;
      }
      case 3: {
        u = unif_rand() * (nedge + nnode * delta_out);
        if (u < nedge) {
          if (directed) {
            start_node[nedge] = start_node[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              start_node[nedge] = start_node[floor(u)];
            } 
            else {
              start_node[nedge] = end_node[floor(u)];
            }
          }
        } 
        else {
          start_node[nedge] = ceil((u - nedge) / delta_out);
        }
        nnode++;
        end_node[nedge] = nnode;
        break;
      }
      case 4: {
        nnode += 2;
        start_node[nedge] = nnode - 1;
        end_node[nedge] = nnode;
        break;
      }
      case 5: {
        nnode++;
        start_node[nedge] = nnode;
        end_node[nedge] = nnode;
        break;
      }
    }
    nedge++;
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["start_node"] = start_node;
  ret["end_node"] = end_node;
  ret["nnode"] = nnode;
  return ret;
}
