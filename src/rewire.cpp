#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Degree preserving rewiring process for directed networks.
//'
//' @param targetNode Target node sequence.
//' @param index_s Index sequence of source nodes' out- and in-degree. 
//'   index_s and index_t bridge the source/target node's and the 
//'   joint degree distribution matrix joint_e.
//' @param index_t Index sequence of target nodes' out- and in-degree. 
//' @param nattempts Integer, number of rewiring attempts.
//' @param joint_e Matrix, target joint distribution generated by
//'   wdnet::directed_edge_level_dist().
//' @return Target node sequence and number of 
//'   accepted rewiring attempts.
//'
// [[Rcpp::export]]
Rcpp::List directed_rewire_cpp(arma::vec targetNode,
                               arma::vec index_s,
                               arma::vec index_t,
                               int nattempts, 
                               arma::mat joint_e) {
  GetRNGstate();
  int accepted = 0, p1 = 0, nedge = targetNode.size();
  int e1, e2, temp;
  int s1, s2, t1, t2;
  double u, ratio;
  for (int i = 0; i < nattempts; i ++) {
    e1 = floor(unif_rand() * nedge);
    e2 = floor(unif_rand() * nedge);
    while (e1 == e2) {
      e2 = floor(unif_rand() * nedge);
    }
    s1 = index_s[e1];
    s2 = index_s[e2];
    t1 = index_t[e1];
    t2 = index_t[e2];
    ratio = joint_e(s1, t2) * joint_e(s2, t1) / 
      (joint_e(s1, t1) * joint_e(s2, t2));
    u = unif_rand();
    if (ratio == 1) {
      p1++;
    }
    if (u < ratio) {
      accepted++;
      temp = index_t[e1];
      index_t[e1] = index_t[e2];
      index_t[e2] = temp;
      temp = targetNode[e1];
      targetNode[e1] = targetNode[e2];
      targetNode[e2] = temp;
    }
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["targetNode"] = targetNode;
  ret["accepted"] = accepted;
  ret["index_t"] = index_t;
  ret["p1"] = p1;
  return ret;
}

//' Degree preserving rewiring process for undirected networks.
//'
//' @param node_s Node sequence of the first column of edgelist.
//' @param node_t Node sequence of the second column of edgelist.
//' @param index_s Index sequence of the first column of edgelist. 
//'   index_s and index_t bridge the nodes' degree and the 
//'   edge-level distribution e.
//' @param index_t Index sequence of the second column of edgelist.
//' @param nattempts Integer, number of rewiring attempts.
//' @param e Matrix, target edge-level distribution generated by
//'   wdnet::undirected_edge_level_dist().
//' @return Node sequences, node indexes and number of 
//'   accepted rewiring attempts.
//'
// [[Rcpp::export]]
Rcpp::List undirected_rewire_cpp(arma::vec node_s,
                                 arma::vec node_t,
                                 arma::vec index_s,
                                 arma::vec index_t,
                                 int nattempts, 
                                 arma::mat e) {
  GetRNGstate();
  int accepted = 0, p1 = 0, nedge = node_s.size();
  int e1, e2, temp;
  int s1, s2, t1, t2;
  double u, v, ratio;
  for (int i = 0; i < nattempts; i ++) {
    e1 = floor(unif_rand() * nedge);
    e2 = floor(unif_rand() * nedge);
    while (e1 == e2) {
      e2 = floor(unif_rand() * nedge);
    }
    s1 = index_s[e1];
    s2 = index_s[e2];
    t1 = index_t[e1];
    t2 = index_t[e2];
    v = unif_rand();
    u = unif_rand();
    if (v < 0.5) {
      ratio = e(s1, t2) * e(s2, t1) / 
        (e(s1, t1) * e(s2, t2));
      if (u < ratio) {
        accepted++;
        temp = index_t[e1];
        index_t[e1] = index_t[e2];
        index_t[e2] = temp;
        temp = node_t[e1];
        node_t[e1] = node_t[e2];
        node_t[e2] = temp;
      }
    } else {
      ratio = e(s1, s2) * e(t1, t2) / 
        (e(s1, t1) * e(s2, t2));
      if (u < ratio) {
        accepted++;
        temp = index_t[e1];
        index_t[e1] = index_s[e2];
        index_s[e2] = temp;
        temp = node_t[e1];
        node_t[e1] = node_s[e2];
        node_s[e2] = temp;
      }
    }
    if (ratio == 1) {
      p1++;
    }
  }
  PutRNGstate();
  
  Rcpp::List ret;
  ret["node_s"] = node_s;
  ret["node_t"] = node_t;
  ret["accepted"] = accepted;
  ret["index_s"] = index_s;
  ret["index_t"] = index_t;
  ret["p1"] = p1;
  return ret;
}