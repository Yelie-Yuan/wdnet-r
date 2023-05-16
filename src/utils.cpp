#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Fill missing nodes in the node sequence. Defined for \code{wdnet::rpanet}.
//'
//' @param nodes Source/target nodes, missing nodes are denoted as 0.
//' @param edges Sampled edges according to preferential attachment.
//' @return Source/target nodes.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec find_node_cpp(arma::vec nodes, 
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

//' Fill missing values in node sequence. Defined for \code{wdnet::rpanet}.
//'
//' @param node1 Nodes in the first column of edgelist, i.e., \code{edgelist[, 1]}.
//' @param node2 Nodes in the second column of edgelist, i.e., \code{edgelist[, 2]}.
//' @param start_edge Index of sampled edges, corresponds to the missing nodes in node1 and node2.
//' @param end_edge Index of sampled edges, corresponds to the missing nodes in node1 and node2.
//' @return Node sequence.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List find_node_undirected_cpp(arma::vec node1, 
                       arma::vec node2, 
                       arma::vec start_edge, 
                       arma::vec end_edge) {
  // GetRNGstate();
  int n = node1.size(), n1 = 0, n2 = 0;
  double u;
  for (int j = 0; j < n; j++) {
    if (node1[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        node1[j] = node1[start_edge[n1]  - 1];
      } else {
        node1[j] = node2[start_edge[n1] - 1];
      }
      n1++;
    }
    if (node2[j] == 0) {
      u = unif_rand();
      if (u <= 0.5) {
        node2[j] = node1[end_edge[n2] - 1];
      } else {
        node2[j] = node2[end_edge[n2] - 1];
      }
      n2++;
    }
  }
  // PutRNGstate();
  
  Rcpp::List ret;
  ret["node1"] = node1;
  ret["node2"] = node2;
  return ret;
}

//' Aggregate edgeweight into nodes' strength.
//'
//' @param snode Source nodes.
//' @param tnode Target nodes.
//' @param weight Edgeweight.
//' @param nnode Number of nodes.
//' @param weighted Logical, true if the edges are weighted, 
//'   false if not.
//' @return Out-strength and in-strength.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List node_strength_cpp(arma::vec snode, 
                            arma::vec tnode,
                            arma::vec weight, 
                            int nnode, 
                            bool weighted = true) {
  int n = snode.size();
  arma::vec outs(nnode, arma::fill::zeros);
  arma::vec ins(nnode, arma::fill::zeros);
  if (weighted) {
    for (int i = 0; i < n; i++) {
      outs[snode[i] - 1] += weight[i];
      ins[tnode[i] - 1] += weight[i];
    }
  } else {
    for (int i = 0; i < n; i++) {
      outs[snode[i] - 1] += 1;
      ins[tnode[i] - 1] += 1;
    }
  }
  
  Rcpp::List ret;
  ret["outs"] = outs;
  ret["ins"] = ins;
  return ret;
}

//' Uniformly draw a node from existing nodes for each time step.
//' Defined for \code{wdnet::rpanet()}.
//'
//' @param total_node Number of existing nodes at each time step.
//' @return Sampled nodes.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::vec sample_node_cpp(arma::vec total_node) {
  // GetRNGstate();
  int n = total_node.size();
  arma::vec nodes(n, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    nodes[i] = Rcpp::sample(total_node[i], 1)[0];
  }
  // PutRNGstate();
  return nodes;
}

//' Fill edgeweight into the adjacency matrix.
//' Defined for function \code{edgelist_to_adj}.
//'
//' @param adj An adjacency matrix.
//' @param edgelist A two column matrix represents the edgelist.
//' @param edgeweight A vector represents the weight of edges.
//' @return Adjacency matrix with edge weight.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat fill_weight_cpp(arma::mat adj, arma::mat edgelist, arma::vec edgeweight) {
  int n = edgelist.n_rows;
  for (int i = 0; i < n; i++) {
    adj(edgelist(i, 0), edgelist(i, 1)) += edgeweight(i);
  }
  return adj;
}
