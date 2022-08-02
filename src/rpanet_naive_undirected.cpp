#include<iostream>
#include<queue>
#include<R.h>
#include<deque>
#include<algorithm>
#include "funcPtrUnd.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;

/**
 * Defult preference function.
 *
 * @param strength Node strength.
 * @param params Parameters passed to the preference function.
 *
 * @return Preference of a node.
 */
double prefFuncDefaultNaive(double strength, Rcpp::NumericVector params) {
  return pow(strength, params[0]) + params[1];
}

/**
 * Calculate node preference.
 *
 * @param func_type Default or customized preference function.
 * @param strength Node strength.
 * @param params Parameters passed to the preference function.
 * @param prefFuncCppNaive Source preference function.
 *
 * @return Node preference.
 */
double calPrefNaive(int func_type, 
                    double strength,
                    Rcpp::NumericVector params, 
                    funcPtrUnd prefFuncCppNaive) {
  if (func_type == 1) {
    return prefFuncDefaultNaive(strength, params);
  }
  else {
    return prefFuncCppNaive(strength);
  }
}

/**
 * Sample an existing node.
 *
 * @param pref Sequence of node preference.
 * @param total_pref Total preference of existing nodes.
 * @param qm Nodes to be excluded from the sampling process.
 *
 * @return Sampled node.
 */
int sampleNodeUndNaive(Rcpp::NumericVector pref, double total_pref, deque<int> &qm) {
  double w;
  int i;
  while (true) {
    i = 0;
    w = 1;
    while (w == 1) {
      w = unif_rand();
    }
    w *= total_pref;
    while (w > 0) {
      w -= pref[i];
      i += 1;
    }
    i -= 1;
    if (find(qm.begin(), qm.end(), i) == qm.end()) {
      return i;
    }
  }
}

//' Preferential attachment algorithm.
//'
//' @param nstep Number of steps.
//' @param m Number of new edges in each step.
//' @param new_node_id New node ID.
//' @param new_edge_id New edge ID.
//' @param node_vec1 Sequence of nodes in the first column of edgelist.
//' @param node_vec2 Sequence of nodes in the second column of edgelist.
//' @param strength Sequence of node strength.
//' @param edgeweight Weight of existing and new edges.
//' @param scenario Scenario of existing and new edges.
//' @param pref Sequence of node preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
// [[Rcpp::export]]
Rcpp::List rpanet_naive_undirected_cpp(
    int nstep, 
    Rcpp::IntegerVector m,
    int new_node_id, 
    int new_edge_id, 
    Rcpp::IntegerVector node_vec1, 
    Rcpp::IntegerVector node_vec2, 
    Rcpp::NumericVector strength, 
    Rcpp::NumericVector edgeweight, 
    Rcpp::IntegerVector scenario,
    Rcpp::NumericVector pref, 
    Rcpp::List control) {
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  Rcpp::List newedge_ctl = control["newedge"];
  bool node_unique = ! newedge_ctl["node.replace"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector params(2);
  funcPtrUnd prefFuncCppNaive;
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type) {
  case 1: 
    params = preference_ctl["params"];
    break;
  case 2: {
      SEXP pref_func_ptr = preference_ctl["pref.pointer"];
      Rcpp::XPtr<funcPtrUnd> xpfun(pref_func_ptr);
      prefFuncCppNaive = *xpfun;
      break;
    }
  }

  double u, total_pref = 0;
  bool m_error;
  int i, j, k, n_existing, current_scenario;
  int node1, node2, temp_node;
  queue<int> q1;
  deque<int> qm;
  for (i = 0; i < new_node_id; i++) {
    pref[i] = calPrefNaive(func_type, strength[i], params, prefFuncCppNaive);
    total_pref += pref[i];
  }
  // sample edges
  GetRNGstate();
  for (i = 0; i < nstep; i++) {
    m_error = false;
    n_existing = new_node_id;
    for (j = 0; j < m[i]; j++) {
      u = unif_rand();
      if (u <= alpha) {
        current_scenario = 1;
      }
      else if (u <= alpha + beta) {
        current_scenario = 2;
      }
      else if (u <= alpha + beta + gamma) {
        current_scenario = 3;
      }
      else if (u <= alpha + beta + gamma + xi) {
        current_scenario = 4;
      }
      else {
        current_scenario = 5;
      }
      if (node_unique) {
        k = qm.size();
        switch (current_scenario) {
          case 1:
            if (k + 1 > n_existing) {
              m_error = true;
            }
            break;
          case 2:
            if (k + 2 - int(beta_loop) > n_existing) {
              m_error = true;
            }
            break;
          case 3:
            if (k + 1 > n_existing) {
              m_error = true;
            }
            break;
        }
      }
      if (m_error) {
        break;
      }
      switch (current_scenario) {
        case 1:
          node1 = new_node_id;
          new_node_id++;
          node2 = sampleNodeUndNaive(pref, total_pref, qm);
          break;
        case 2:
          node1 = sampleNodeUndNaive(pref, total_pref, qm);
          if (! beta_loop) {
            qm.push_back(node1);
            node2 = sampleNodeUndNaive(pref, total_pref, qm);
            qm.pop_back();
          }
          else {
            node2 = sampleNodeUndNaive(pref, total_pref, qm);
          }
          break;
        case 3:
          node1 = sampleNodeUndNaive(pref, total_pref, qm);
          node2 = new_node_id;
          new_node_id++;
          break;
        case 4:
          node1 = new_node_id;
          new_node_id++;
          node2 = new_node_id;
          new_node_id++;
          break;
        case 5:
          node1 = node2 = new_node_id;
          new_node_id++;
          break;
      }
      // handle duplicate nodes
      if (node_unique) {
        if (node1 < n_existing) {
          qm.push_back(node1);
        }
        if ((node2 < n_existing) && (node1 != node2)) {
          qm.push_back(node2);
        }
      }
      strength[node1] += edgeweight[new_edge_id];
      strength[node2] += edgeweight[new_edge_id];
      node_vec1[new_edge_id] = node1;
      node_vec2[new_edge_id] = node2;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      new_edge_id++;
    }
    if (m_error) {
      m[i] = j;
      // need to print this info
      Rprintf("Unique nodes exhausted at step %u. Set the value of m at current step to %u.\n", i + 1, j);
    }
    while (! q1.empty()) {
      temp_node = q1.front();
      total_pref -= pref[temp_node];
      pref[temp_node] = calPrefNaive(func_type, strength[temp_node], params, prefFuncCppNaive);
      total_pref += pref[temp_node];
      q1.pop();
    }
    qm.clear();
  }
  PutRNGstate();
  // check total preference = sum of node preference
  // Rprintf("Total pref %f.\n", total_pref);
  // for (i = 0; i < new_node_id; i++) {
  //   total_pref -= pref[i];
  // }
  // Rprintf("Total pref %f.\n", total_pref * pow(10, 10));

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = node_vec1;
  ret["node_vec2"] = node_vec2;
  ret["pref"] = pref;
  ret["strength"] = strength;
  ret["scenario"] = scenario;
  return ret;
}