#include<iostream>
#include<queue>
#include<R.h>
#include "funcPtrUnd.h"
#include<Rcpp.h>
using namespace std;
funcPtrUnd prefFuncCppNaive;

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
 * @param params Parameters passed to the default preference function.
 * @param prefFuncCppNaive Pointer of the customized source preference function.
 *
 * @return Node preference.
 */
double calcPrefNaive(int func_type, 
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
 * @param n_existing Number of existing nodes.
 * @param pref Sequence of node preference.
 * @param total_pref Total preference of existing nodes.
 *
 * @return Sampled node.
 */
int sampleNodeUndNaive(int n_existing, Rcpp::NumericVector pref, double total_pref) {
  double w = 1;
  int i = 0;
  while (w == 1) {
    w = unif_rand();
  }
  w *= total_pref;
  while ((w > 0) && (i < n_existing)) {
    w -= pref[i];
    i += 1;
  }
  if (w > 0) {
    Rprintf("Numerical error! Returning the last node (node %d) as the sampled node. \n", n_existing);
    // i = n_existing;
  }
  return i - 1;
}

// /**
//  * Calculate total preference.
//  * 
//  * @param pref Preference vector.
//  * @param n_exising Number of existing nodes.
//  * 
//  * @return Total preference.
//  * 
//  */
// double calcTotalprefUnd(Rcpp::NumericVector pref, int n_existing) {
//   int k;
//   double temp = 0;
//   for (k = 0; k < n_existing; k++) {
//     temp += pref[k];
//   }
//   return temp;
// }

// /**
//  * Check difference.
//  * 
//  * @param total_pref Total preference.
//  * @param pref Preference vector.
//  * 
//  */
// void checkDiffUnd(Rcpp::NumericVector pref, double total_pref) {
//   int k;
//   double temp = 0, tol = 0.00000001;
//   for (k = 0; k < pref.size(); k++) {
//     temp += pref[k];
//   }
//   if ((total_pref - temp > tol) || (temp - total_pref) > tol) {
//     Rprintf("Total pref warning, diff = %f. \n", total_pref - temp);
//   }
// }

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
//' @keywords internal
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
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type) {
  case 1: 
    params = preference_ctl["params"];
    break;
  case 2: {
      SEXP pref_func_ptr = preference_ctl["pref.pointer"];
      prefFuncCppNaive = *Rcpp::XPtr<funcPtrUnd>(pref_func_ptr);
      break;
    }
  }

  double u, total_pref = 0, temp_p;
  bool m_error;
  int i, j, k, n_existing, current_scenario;
  int node1, node2, temp_node;
  queue<int> q1;
  for (i = 0; i < new_node_id; i++) {
    pref[i] = calcPrefNaive(func_type, strength[i], params, prefFuncCppNaive);
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
        if (current_scenario <= 3) {
          // check whether sum(pref) == 0
          for (k = 0; k < n_existing; k++) {
            if (pref[k] > 0) {
              break;
            }
          }
          if (k == n_existing) {
            total_pref = 0;
            m_error = true;
            break;
          }
        }
      }
      switch (current_scenario) {
        case 1:
          node1 = new_node_id;
          new_node_id++;
          node2 = sampleNodeUndNaive(n_existing, pref, total_pref);
          break;
        case 2:
          node1 = sampleNodeUndNaive(n_existing, pref, total_pref);
          if (! beta_loop) {
            if (pref[node1] == total_pref) {
              m_error = true;
              break;
            }
            temp_p = pref[node1];
            pref[node1] = 0;
            total_pref -= temp_p;
            // check whether sum(pref) == 0
            for (k = 0; k < n_existing; k++) {
              if (pref[k] > 0) {
                break;
              }
            }
            if (k == n_existing) {
              total_pref = 0;
              m_error = true;
              break;
            }

            node2 = sampleNodeUndNaive(n_existing, pref, total_pref);
            pref[node1] = temp_p;
            total_pref += temp_p;
          }
          else {
            node2 = sampleNodeUndNaive(n_existing, pref, total_pref);
          }
          break;
        case 3:
          node1 = sampleNodeUndNaive(n_existing, pref, total_pref);
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
      if (m_error) {
        break;
      }
      // sample without replacement
      if (node_unique) {
        if (node1 < n_existing) {
          total_pref -= pref[node1];
          pref[node1] = 0;
        }
        if ((node2 < n_existing) && (node1 != node2)) {
          total_pref -= pref[node2];
          pref[node2] = 0;
        }
      }
      // checkDiffUnd(pref, total_pref);
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
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n", current_scenario, i + 1, j);
    }
    while (! q1.empty()) {
      temp_node = q1.front();
      total_pref -= pref[temp_node];
      pref[temp_node] = calcPrefNaive(func_type, strength[temp_node], params, prefFuncCppNaive);
      total_pref += pref[temp_node];
      q1.pop();
    }
    // checkDiffUnd(pref, total_pref);
  }
  PutRNGstate();

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