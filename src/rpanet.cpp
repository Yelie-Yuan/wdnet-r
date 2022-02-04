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
//' @param startNode Sequence of nodes from the first column of edgelist, i.e., \code{edgelist[, 1]}.
//' @param endNode Sequence of nodes from the first column of edgelist, i.e., \code{edgelist[, 2]}.
//' @param startEdge Index of sampled edges, corresponds to the missing nodes in startNode.
//' @param endEdge Index of sampled edges, corresponds to the missing nodes in endNode.
//' @return Sequence of source/target node of edges.
// [[Rcpp::export]]
Rcpp::List findNode_undirected_cpp(arma::vec startNode, 
                       arma::vec endNode, 
                       arma::vec startEdge, 
                       arma::vec endEdge) {
  GetRNGstate();
  int n = startNode.size(), n1 = 0, n2 = 0;
  double u;
  for (int j = 0; j < n; j++) {
    if (startNode[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        startNode[j] = startNode[startEdge[n1]  - 1];
      } else {
        startNode[j] = endNode[startEdge[n1] - 1];
      }
      n1++;
    }
    if (endNode[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        endNode[j] = startNode[endEdge[n2] - 1];
      } else {
        endNode[j] = endNode[endEdge[n2] - 1];
      }
      n2++;
    }
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["startNode"] = startNode;
  ret["endNode"] = endNode;
  return ret;
}

//' Aggregate edgeweight into nodes' strength.
//'
//' @param startNode Sequence of source nodes.
//' @param endNode Sequence of target nodes.
//' @param weight Sequence of edgeweight.
//' @param nNodes Number of nodes of sampled network.
//' @param weighted Logical, ture if the edges are weighted, 
//'   false if not.
//' @return Sequence of outstrength/instrength.
// [[Rcpp::export]]
Rcpp::List nodeStrength_cpp(arma::vec startNode, 
                            arma::vec endNode,
                            arma::vec weight, 
                            int nNodes, 
                            bool weighted = true) {
  int n = startNode.size();
  arma::vec outstrength(nNodes, arma::fill::zeros);
  arma::vec instrength(nNodes, arma::fill::zeros);
  if (weighted) {
    for (int i = 0; i < n; i++) {
      outstrength[startNode[i] - 1] += weight[i];
      instrength[endNode[i] - 1] += weight[i];
    }
  } else {
    for (int i = 0; i < n; i++) {
      outstrength[startNode[i] - 1] += 1;
      instrength[endNode[i] - 1] += 1;
    }
  }
  
  Rcpp::List ret;
  ret["outstrength"] = outstrength;
  ret["instrength"] = instrength;
  return ret;
}

//' Sample nodes with respect to the number of nodes at each step.
//'
//' @param totalNode Number of nodes at each step.
//' @return Sequence of sampled nodes.
// [[Rcpp::export]]
arma::vec sampleNode_cpp(arma::vec totalNode) {
  GetRNGstate();
  int n = totalNode.size();
  arma::vec nodes(n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    nodes[i] = Rcpp::sample(totalNode[i], 1)[0];
  }
  PutRNGstate();
  return nodes;
}


//' Preferential attachment algorithm for simple situations, 
//' e.g., edge weights are constant, 
//' number of new edges per step is 1.
//'
//' @param startNode Sequence of source nodes.
//' @param endNode Sequence of target nodes.
//' @param scenario Sequence of alpha, beta, gamma, xi, rho scenarios.
//' @param nNodes Number of nodes at current step.
//' @param nEdges Number of edges at current step.
//' @param delta_out Tuning parameter.
//' @param delta_in Tuning parameter.
//' @param directed Whether the network is directed.
//' @return Number of nodes, sequences of source and target nodes.
//' 
// [[Rcpp::export]]
Rcpp::List rpanet_cpp(arma::vec startNode,
                      arma::vec endNode,
                      arma::vec scenario,
                      int nNodes,
                      int nEdges,
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
        u = unif_rand() * (nEdges + nNodes * delta_in);
        if (u < nEdges) {
          if (directed) {
            endNode[nEdges] = endNode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              endNode[nEdges] = startNode[floor(u)];
            } 
            else {
              endNode[nEdges] = endNode[floor(u)];
            }
          }
        }
        else {
          endNode[nEdges] = ceil((u - nEdges) / delta_in);
        }
        nNodes++;
        startNode[nEdges] = nNodes;
        break;
      }
      case 2: {
        u = unif_rand() * (nEdges + nNodes * delta_out);
        if (u < nEdges) {
          if (directed) {
            startNode[nEdges] = startNode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              startNode[nEdges] = startNode[floor(u)];
            } 
            else {
              startNode[nEdges] = endNode[floor(u)];
            }
          }
        } 
        else {
          startNode[nEdges] = ceil((u - nEdges) / delta_out);
        }
        
        u = unif_rand() * (nEdges + nNodes * delta_in);
        if (u < nEdges) {
          if (directed) {
            endNode[nEdges] = endNode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              endNode[nEdges] = startNode[floor(u)];
            } 
            else {
              endNode[nEdges] = endNode[floor(u)];
            }
          }
        } 
        else {
          endNode[nEdges] = ceil((u - nEdges) / delta_in);
        }
        break;
      }
      case 3: {
        u = unif_rand() * (nEdges + nNodes * delta_out);
        if (u < nEdges) {
          if (directed) {
            startNode[nEdges] = startNode[floor(u)] ;
          }
          else {
            v = unif_rand();
            if (v <= 0.5) {
              startNode[nEdges] = startNode[floor(u)];
            } 
            else {
              startNode[nEdges] = endNode[floor(u)];
            }
          }
        } 
        else {
          startNode[nEdges] = ceil((u - nEdges) / delta_out);
        }
        nNodes++;
        endNode[nEdges] = nNodes;
        break;
      }
      case 4: {
        nNodes += 2;
        startNode[nEdges] = nNodes - 1;
        endNode[nEdges] = nNodes;
        break;
      }
      case 5: {
        nNodes++;
        startNode[nEdges] = nNodes;
        endNode[nEdges] = nNodes;
        break;
      }
    }
    nEdges++;
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["startNode"] = startNode;
  ret["endNode"] = endNode;
  ret["nNodes"] = nNodes;
  return ret;
}
