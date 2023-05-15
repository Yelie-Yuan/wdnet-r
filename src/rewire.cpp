#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//' Degree preserving rewiring process for directed networks.
//'
//' @param iteration Integer, number of iterations for rewiring attempts.
//' @param nattempts Integer, number of rewiring attempts per iteration.
//' @param tnode Vector, target node sequence - 1.
//' @param sout Vector, source nodes' out-degree.
//' @param sin Vector, source nodes' in-degree.
//' @param tout Vector, target nodes' out-degree.
//' @param tin Vector, target nodes' in-degree.
//' @param index_s Index of source nodes' out- and in-degree. 
//'   \code{index_s}/\code{index_t} bridges the indices of source/target nodes and the 
//'   target structure eta.
//' @param index_t Index of target nodes' out- and in-degree. 
//' @param eta Matrix, target structure eta generated by
//'   \code{wdnet::get_eta_directed()}.
//' @param rewire_history Logical, whether the rewiring history should be returned.
//' @return Returns target node sequence, four directed assortativity coefficients after each iteration, and rewire history.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List dprewire_directed_cpp(
    int iteration, 
    int nattempts,
    arma::uvec tnode,
    arma::vec sout,
    arma::vec sin, 
    arma::vec tout, 
    arma::vec tin,
    arma::uvec index_s,
    arma::uvec index_t,
    arma::mat eta, 
    bool rewire_history) {
  GetRNGstate();
  arma::vec outout(iteration, arma::fill::zeros);
  arma::vec outin(iteration, arma::fill::zeros);
  arma::vec inout(iteration, arma::fill::zeros);
  arma::vec inin(iteration, arma::fill::zeros);
  // arma::vec r_out_out(iteration, arma::fill::zeros);
  // arma::vec r_out_in(iteration, arma::fill::zeros);
  // arma::vec r_in_out(iteration, arma::fill::zeros);
  // arma::vec r_in_in(iteration, arma::fill::zeros);
  int nedge = tnode.size();
  int e1, e2, count = 0;
  int s1, s2, t1, t2, hist_row;
  double u, ratio, temp;
  if (rewire_history) {
    hist_row = iteration * nattempts;
  } else {
    hist_row = 1;
  }
  arma::mat history(hist_row, 4, arma::fill::zeros);
  for (int n = 0; n < iteration; n++) {
    for (int i = 0; i < nattempts; i++) {
      e1 = floor(unif_rand() * nedge);
      e2 = floor(unif_rand() * nedge);
      while (e1 == e2) {
        e2 = floor(unif_rand() * nedge);
      }
      if (rewire_history) {
        history(count, 0) = count;
        history(count, 1) = e1;
        history(count, 2) = e2;
      }
      s1 = index_s[e1];
      s2 = index_s[e2];
      t1 = index_t[e1];
      t2 = index_t[e2];
      if ((eta(s1, t2) * eta(s2, t1)) < (eta(s1, t1) * eta(s2, t2))) {
        ratio = eta(s1, t2) * eta(s2, t1) / 
          (eta(s1, t1) * eta(s2, t2));
      }
      else {
        ratio = 1;
      }
      u = unif_rand();
      if (u <= ratio) {
        temp = index_t[e1];
        index_t[e1] = index_t[e2];
        index_t[e2] = temp;
        temp = tnode[e1];
        tnode[e1] = tnode[e2];
        tnode[e2] = temp;
        temp = tout[e1];
        tout[e1] = tout[e2];
        tout[e2] = temp;
        temp = tin[e1];
        tin[e1] = tin[e2];
        tin[e2] = temp;
        // temp = r_targetOut[e1];
        // r_targetOut[e1] = r_targetOut[e2];
        // r_targetOut[e2] = temp;
        // temp = r_targetIn[e1];
        // r_targetIn[e1] = r_targetIn[e2];
        // r_targetIn[e2] = temp;
        if (rewire_history) {
          history(count, 3) = 1;
        }
      }
      count++;
    }
    outout[n] = (arma::cor(sout, tout)).eval()(0, 0);
    outin[n] = (arma::cor(sout, tin)).eval()(0, 0);
    inout[n] = (arma::cor(sin, tout)).eval()(0, 0);
    inin[n] = (arma::cor(sin, tin)).eval()(0, 0);
    // r_out_out[n] = (arma::cor(r_sourceOut, r_targetOut)).eval()(0, 0);
    // r_out_in[n] = (arma::cor(r_sourceOut, r_targetIn)).eval()(0, 0);
    // r_in_out[n] = (arma::cor(r_sourceIn, r_targetOut)).eval()(0, 0);
    // r_in_in[n] = (arma::cor(r_sourceIn, r_targetIn)).eval()(0, 0);
  }
  
  PutRNGstate();
  Rcpp::List ret;
  ret["tnode"] = tnode;
  if (rewire_history) {
    ret["history"] = history;
  }
  ret["outout"] = outout;
  ret["outin"] = outin;
  ret["inout"] = inout;
  ret["inin"] = inin;
  // ret["r_out_out"] = r_out_out;
  // ret["r_out_in"] = r_out_in;
  // ret["r_in_out"] = r_in_out;
  // ret["r_in_in"] = r_in_in;
  return ret;
}

//' Degree preserving rewiring process for undirected networks.
//' 
//' @param iteration Integer, number of iterations for rewiring attempts.
//' @param nattempts Integer, number of rewiring attempts per iteration.
//' @param node1 Vector, first column of edgelist.
//' @param node2 Vector, second column of edgelist.
//' @param degree1 Vector, degree of node1 and node2.
//' @param degree2 Vector, degree of node2 and node1. degree1 
//'   and degree2 are used to calculate assortativity coefficient,
//'   i.e., degree correlation.
//' @param index1 Index of the first column of edgelist. 
//'   \code{index1} and \code{index2} bridge the nodes' degree and the 
//'   structure \code{e}.
//' @param index2 Index of the second column of edgelist.
//' @param e Matrix, target structure (eta) generated by
//'   \code{wdnet::get_eta_undirected()}.
//' @param rewire_history Logical, whether the rewiring history should be returned.
//' @return Returns node sequences, assortativity coefficient after each iteration, and rewiring history.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List dprewire_undirected_cpp(
    int iteration,
    int nattempts,
    Rcpp::IntegerVector node1, 
    Rcpp::IntegerVector node2, 
    arma::vec degree1,
    arma::vec degree2,
    arma::vec index1,
    arma::vec index2,
    arma::mat e, 
    bool rewire_history) {
  GetRNGstate();
  arma::vec rho(iteration, arma::fill::zeros);
  int nedge = index1.size();
  int e1, e2, temp, count = 0;
  int s1, s2, t1, t2, hist_row;
  double u, v, ratio;
  if (rewire_history) {
    hist_row = iteration * nattempts;
  } else {
    hist_row = 1;
  }
  arma::mat history(hist_row, 5, arma::fill::zeros);
  
  for (int n = 0; n < iteration; n++) {
    for (int i = 0; i < nattempts; i++) {
      e1 = floor(unif_rand() * nedge);
      e2 = floor(unif_rand() * nedge);
      while (e1 == e2) {
        e2 = floor(unif_rand() * nedge);
      }
      if (rewire_history) {
        history(count, 0) = count;
        history(count, 1) = e1;
        history(count, 2) = e2;
      }
      s1 = index1[e1];
      s2 = index1[e2];
      t1 = index2[e1];
      t2 = index2[e2];
      v = unif_rand();
      u = unif_rand();
      if (v < 0.5) {
        // if (rewire_history) {
        //   history(count, 3) = 0;
        // }
        if ((e(s1, t2) * e(s2, t1)) < (e(s1, t1) * e(s2, t2))) {
          ratio = e(s1, t2) * e(s2, t1) / 
            (e(s1, t1) * e(s2, t2));
        }
        else {
          ratio = 1;
        }
        if (u <= ratio) {
          if (rewire_history) {
            history(count, 4) = 1;
          }
          temp = index2[e1];
          index2[e1] = index2[e2];
          index2[e2] = temp;
          temp = node2[e1];
          node2[e1] = node2[e2];
          node2[e2] = temp;
          temp = degree2[e1];
          degree2[e1] = degree2[e2];
          degree2[e2] = temp;
          temp = degree1[e1 + nedge];
          degree1[e1 + nedge] = degree1[e2 + nedge];
          degree1[e2 + nedge] = temp;
        }
      } else {
        if (rewire_history) {
          history(count, 3) = 1;
        }
        if ((e(s1, s2) * e(t1, t2)) < (e(s1, t1) * e(s2, t2))) {
          ratio = e(s1, s2) * e(t1, t2) / 
            (e(s1, t1) * e(s2, t2));
        }
        else {
          ratio = 1;
        }
        if (u <= ratio) {
          if (rewire_history) {
            history(count, 4) = 1;
          }
          temp = index2[e1];
          index2[e1] = index1[e2];
          index1[e2] = temp;
          temp = node2[e1];
          node2[e1] = node1[e2];
          node1[e2] = temp;
          temp = degree2[e1];
          degree2[e1] = degree1[e2];
          degree1[e2] = temp;
          temp = degree1[e1 + nedge];
          degree1[e1 + nedge] = degree2[e2 + nedge];
          degree2[e2 + nedge] = temp;
        }
      }
      count++;
    }
    rho[n] = (arma::cor(degree1, degree2)).eval()(0, 0);
  }
  PutRNGstate();
  Rcpp::List ret;
  if (rewire_history) {
    ret["history"] = history;
  }
  ret["node1"] = node1;
  ret["node2"] = node2;
  ret["rho"] = rho;
  // ret["degree1"] = degree1;
  // ret["degree2"] = degree2;
  return ret;
}
