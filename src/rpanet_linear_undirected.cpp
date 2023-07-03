#include <iostream>
#include <queue>
#include <R.h>
#include <Rcpp.h>
#include "rpanet_binary_linear.h"

using namespace std;
funcPtrUnd custmPrefLinear;

/**
 * Calculate node preference.
 *
 * @param func_type Default or customized preference function.
 * @param s Node strength.
 * @param params Parameters passed to the default preference function.
 * @param custmPrefLinear Pointer of the customized source preference function.
 *
 * @return Node preference.
 */
double calcPrefLinearUnd(int func_type,
                      double s,
                      double *params,
                      funcPtrUnd custmPrefLinear)
{
  double ret;
  if (func_type == 1)
  {
    ret = prefFuncUnd(s, params);
  }
  else
  {
    ret = custmPrefLinear(s);
  }
  if (ret < 0)
  {
    Rcpp::stop("Negative preference score returned, please check your preference function(s).");
  }
  return ret;
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
//' @param s Sequence of node strength.
//' @param edgeweight Weight of existing and new edges.
//' @param scenario Scenario of existing and new edges.
//' @param pref Sequence of node preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List rpanet_linear_undirected_cpp(
    int nstep,
    Rcpp::IntegerVector m,
    int new_node_id,
    int new_edge_id,
    Rcpp::IntegerVector node_vec1,
    Rcpp::IntegerVector node_vec2,
    Rcpp::NumericVector s,
    Rcpp::NumericVector edgeweight,
    Rcpp::IntegerVector scenario,
    Rcpp::NumericVector pref_vec,
    Rcpp::List control)
{
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  Rcpp::List newedge_ctl = control["newedge"];
  bool node_unique = !newedge_ctl["node.replace"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector params_vec(2);
  double *params = nullptr;
  double *pref = &(pref_vec[0]);
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type)
  {
  case 1:
    params_vec = preference_ctl["params"];
    params = &(params_vec[0]);
    break;
  case 2:
  {
    SEXP pref_func_ptr = preference_ctl["pref.pointer"];
    custmPrefLinear = *Rcpp::XPtr<funcPtrUnd>(pref_func_ptr);
    break;
  }
  }

  double u, total_pref = 0, temp_p;
  bool m_error;
  int i, j, k, n_existing, current_scenario;
  int node1, node2, temp_node, n_seednode = new_node_id;

  // update node id; from R to c++
  for (i = 0; i < new_edge_id; i++)
  {
    node_vec1[i] = node_vec1[i] - 1;
    node_vec2[i] = node_vec2[i] - 1;
  }

  // sort nodes according to node preference
  Rcpp::IntegerVector sorted_node_vec = Rcpp::seq(0, n_seednode - 1);
  for (i = 0; i < new_node_id; i++)
  {
    pref[i] = calcPrefLinearUnd(func_type, s[i], params, custmPrefLinear);
    total_pref += pref[i];
  }
  sort(sorted_node_vec.begin(), sorted_node_vec.end(),
       [&](int k, int l){ return pref[k] > pref[l]; });
  int *sorted_node = &(sorted_node_vec[0]);

  // sample edges
  queue<int> q1;
  // GetRNGstate();
  for (i = 0; i < nstep; i++)
  {
    m_error = false;
    n_existing = new_node_id;
    for (j = 0; j < m[i]; j++)
    {
      u = unif_rand();
      if (u <= alpha)
      {
        current_scenario = 1;
      }
      else if (u <= alpha + beta)
      {
        current_scenario = 2;
      }
      else if (u <= alpha + beta + gamma)
      {
        current_scenario = 3;
      }
      else if (u <= alpha + beta + gamma + xi)
      {
        current_scenario = 4;
      }
      else
      {
        current_scenario = 5;
      }
      if (node_unique)
      {
        if (current_scenario <= 3)
        {
          // check whether sum(pref) == 0
          for (k = 0; k < n_existing; k++)
          {
            if (pref[k] > 0)
            {
              break;
            }
          }
          if (k == n_existing)
          {
            total_pref = 0;
            m_error = true;
            break;
          }
        }
      }
      switch (current_scenario)
      {
      case 1:
        node1 = new_node_id;
        new_node_id++;
        node2 = sampleNodeLinear(n_existing, n_seednode, pref, total_pref, sorted_node);
        break;
      case 2:
        node1 = sampleNodeLinear(n_existing, n_seednode, pref, total_pref, sorted_node);
        if (!beta_loop)
        {
          if (pref[node1] == total_pref)
          {
            m_error = true;
            break;
          }
          temp_p = pref[node1];
          pref[node1] = 0;
          total_pref -= temp_p;
          // check whether sum(pref) == 0
          for (k = 0; k < n_existing; k++)
          {
            if (pref[k] > 0)
            {
              break;
            }
          }
          if (k == n_existing)
          {
            total_pref = 0;
            m_error = true;
            break;
          }

          node2 = sampleNodeLinear(n_existing, n_seednode, pref, total_pref, sorted_node);
          pref[node1] = temp_p;
          total_pref += temp_p;
        }
        else
        {
          node2 = sampleNodeLinear(n_existing, n_seednode, pref, total_pref, sorted_node);
        }
        break;
      case 3:
        node1 = sampleNodeLinear(n_existing, n_seednode, pref, total_pref, sorted_node);
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
      if (m_error)
      {
        break;
      }
      // sample without replacement
      if (node_unique)
      {
        if (node1 < n_existing)
        {
          total_pref -= pref[node1];
          pref[node1] = 0;
        }
        if ((node2 < n_existing) && (node1 != node2))
        {
          total_pref -= pref[node2];
          pref[node2] = 0;
        }
      }
      // checkDiffUnd(pref, total_pref);
      s[node1] += edgeweight[new_edge_id];
      s[node2] += edgeweight[new_edge_id];
      node_vec1[new_edge_id] = node1;
      node_vec2[new_edge_id] = node2;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      new_edge_id++;
    }
    if (m_error)
    {
      m[i] = j;
      // need to print this info
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n", current_scenario, i + 1, j);
    }
    while (!q1.empty())
    {
      temp_node = q1.front();
      total_pref -= pref[temp_node];
      pref[temp_node] = calcPrefLinearUnd(func_type, s[temp_node], params, custmPrefLinear);
      total_pref += pref[temp_node];
      q1.pop();
    }
    // checkDiffUnd(pref, total_pref);
  }
  // PutRNGstate();

  // update node id; from c++ to R
  for (i = 0; i < new_edge_id; i++)
  {
    node_vec1[i] = node_vec1[i] + 1;
    node_vec2[i] = node_vec2[i] + 1;
  }

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = node_vec1;
  ret["node_vec2"] = node_vec2;
  ret["pref"] = pref_vec;
  ret["s"] = s;
  ret["scenario"] = scenario;
  return ret;
}
