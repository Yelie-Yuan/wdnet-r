#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Fill missing values in node sequence.
//'
//' @param nodes Sequence of source/target node of edges, missing values are denoted as 0.
//' @param edges Edges sampled according to preferential attachment.
//' @param index Index of missing values in nodes.
//' @return Sequence of source/target node of edges.
// [[Rcpp::export]]
arma::vec nodes_cpp(arma::vec nodes, 
                    arma::vec edges, 
                    arma::vec index) {
  int n = index.size();
  for (int j = 0; j < n; j++) {
    nodes[index[j] - 1] = nodes[edges[j] - 1];
  }
  return nodes;
}

//' Aggregate edgeweight into nodes' strength.
//'
//' @param node Sequence of source/target node of edges.
//' @param weight Sequence of edgeweight.
//' @param nNodes Number of nodes of sampled network.
//' @return Sequence of outstrength/instrength.
// [[Rcpp::export]]
arma::vec strength_cpp(arma::vec node, arma::vec weight, int nNodes) {
  int n = node.size();
  arma::vec strength(nNodes, arma::fill::zeros);
  for (int i = 0; i < n; i++) {
    strength[node[i] - 1] += weight[i];
  }
  return strength;
}