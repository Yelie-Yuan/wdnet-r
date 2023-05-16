#include <iostream>
#include <queue>
#include <R.h>
#include <Rcpp.h>
#include "rpanet_binary_linear.h"

using namespace std;
funcPtrD custmSourcePrefLinear;
funcPtrD custmTargetPrefLinear;

/**
 *  Calculate node source preference.
 *
 * @param func_type Default or customized preference function.
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param params Parameters passed to the default source/target preference function.
 * @param custmPrefLinear Pointer of the customized source/target preference function.
 *
 * @return Node source preference.
 */
double calcPrefLinearD(int func_type,
                       double outs,
                       double ins,
                       double *params,
                       funcPtrD custmPrefLinear)
{
  double ret;
  if (func_type == 1)
  {
    ret = prefFuncD(outs, ins, params);
  }
  else
  {
    ret = custmPrefLinear(outs, ins);
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
// double calcTotalprefD(Rcpp::NumericVector pref, int n_existing) {
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
// void checkDiffD(Rcpp::NumericVector pref, double total_pref) {
//   int k;
//   double temp = 0, tol = 0.00000001;
//   for (k = 0; k < pref.size(); k++) {
//     temp += pref[k];
//   }
//   if ((total_pref - temp > tol) || (temp - total_pref) > tol) {
//     Rprintf("Total pref warning, diff = %f. \n", total_pref - temp);
//   }
// }

//'  Preferential attachment algorithm.
//'
//' @param nstep Number of steps.
//' @param m Number of new edges in each step.
//' @param new_node_id New node ID.
//' @param new_edge_id New edge ID.
//' @param source_node Sequence of source nodes.
//' @param target_node Sequence of target nodes.
//' @param outs Sequence of out-strength.
//' @param ins Sequence of in-strength.
//' @param edgeweight Weight of existing and new edges.
//' @param scenario Scenario of existing and new edges.
//' @param sample_recip Logical, whether reciprocal edges will be added.
//' @param node_group Sequence of node group.
//' @param spref Sequence of node source preference.
//' @param tpref Sequence of node target preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List rpanet_linear_directed_cpp(
    int nstep,
    Rcpp::IntegerVector m,
    int new_node_id,
    int new_edge_id,
    Rcpp::IntegerVector source_node,
    Rcpp::IntegerVector target_node,
    Rcpp::NumericVector outs,
    Rcpp::NumericVector ins,
    Rcpp::NumericVector edgeweight,
    Rcpp::IntegerVector scenario,
    bool sample_recip,
    Rcpp::IntegerVector node_group,
    Rcpp::NumericVector spref_vec,
    Rcpp::NumericVector tpref_vec,
    Rcpp::List control)
{
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  bool source_first = scenario_ctl["source.first"];
  Rcpp::List newedge_ctl = control["newedge"];
  // bool node_unique = ! newedge_ctl["node.replace"];
  bool snode_unique = !newedge_ctl["snode.replace"];
  bool tnode_unique = !newedge_ctl["tnode.replace"];
  Rcpp::List reciprocal_ctl = control["reciprocal"];
  bool selfloop_recip = reciprocal_ctl["selfloop.recip"];
  Rcpp::NumericVector group_prob_vec = reciprocal_ctl["group.prob"];
  double *group_prob = &(group_prob_vec[0]);
  Rcpp::NumericMatrix recip_prob = reciprocal_ctl["recip.prob"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector sparams_vec(5);
  Rcpp::NumericVector tparams_vec(5);
  double *sparams, *tparams;
  double *spref = &(spref_vec[0]);
  double *tpref = &(tpref_vec[0]);
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type)
  {
  case 1:
    sparams_vec = preference_ctl["sparams"];
    tparams_vec = preference_ctl["tparams"];
    sparams = &(sparams_vec[0]);
    tparams = &(tparams_vec[0]);
    break;
  case 2:
  {
    SEXP source_pref_func_ptr = preference_ctl["spref.pointer"];
    custmSourcePrefLinear = *Rcpp::XPtr<funcPtrD>(source_pref_func_ptr);
    SEXP target_pref_func_ptr = preference_ctl["tpref.pointer"];
    custmTargetPrefLinear = *Rcpp::XPtr<funcPtrD>(target_pref_func_ptr);
    break;
  }
  }

  double u, p, temp_p, total_spref = 0, total_tpref = 0;
  bool m_error;
  int i, j, k, n_existing, current_scenario, n_reciprocal;
  int node1, node2, temp_node, n_seednode = new_node_id;

  // sort nodes according to node preference
  Rcpp::IntegerVector sorted_snode_vec = Rcpp::seq(0, n_seednode - 1);
  Rcpp::IntegerVector sorted_tnode_vec = Rcpp::seq(0, n_seednode - 1);
  for (int i = 0; i < new_node_id; i++)
  {
    spref[i] = calcPrefLinearD(func_type, outs[i], ins[i], sparams, custmSourcePrefLinear);
    tpref[i] = calcPrefLinearD(func_type, outs[i], ins[i], tparams, custmTargetPrefLinear);
    total_spref += spref[i];
    total_tpref += tpref[i];
  }
  sort(sorted_snode_vec.begin(), sorted_snode_vec.end(),
       [&](int k, int l){ return spref[k] > spref[l]; });
  sort(sorted_tnode_vec.begin(), sorted_tnode_vec.end(),
       [&](int k, int l){ return tpref[k] > tpref[l]; });
  int *sorted_snode = &(sorted_snode_vec[0]);
  int *sorted_tnode = &(sorted_tnode_vec[0]);

  // sample edges
  queue<int> q1;
  // GetRNGstate();
  for (i = 0; i < nstep; i++)
  {
    n_reciprocal = 0;
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
      if (snode_unique)
      {
        if ((current_scenario == 2) || (current_scenario == 3))
        {
          for (k = 0; k < n_existing; k++)
          {
            if (spref[k] > 0)
            {
              break;
            }
          }
          if (k == n_existing)
          {
            total_spref = 0;
            m_error = true;
            break;
          }
        }
      }
      if (tnode_unique)
      {
        if ((current_scenario == 1) || (current_scenario == 2))
        {
          for (k = 0; k < n_existing; k++)
          {
            if (tpref[k] > 0)
            {
              break;
            }
          }
          if (k == n_existing)
          {
            total_tpref = 0;
            m_error = true;
            break;
          }
        }
      }
      switch (current_scenario)
      {
      case 1:
        node1 = new_node_id;
        if (sample_recip)
        {
          node_group[node1] = sampleGroup(group_prob);
        }
        new_node_id++;
        node2 = sampleNodeLinear(n_existing, n_seednode, tpref, total_tpref, sorted_tnode);
        break;
      case 2:
        if (source_first)
        {
          node1 = sampleNodeLinear(n_existing, n_seednode, spref, total_spref, sorted_snode);
          if (beta_loop)
          {
            node2 = sampleNodeLinear(n_existing, n_seednode, tpref, total_tpref, sorted_tnode);
          }
          else
          {
            if (tpref[node1] == total_tpref)
            {
              m_error = true;
              break;
            }
            if (tpref[node1] == 0)
            {
              node2 = sampleNodeLinear(n_existing, n_seednode, tpref, total_tpref, sorted_tnode);
            }
            else
            {
              temp_p = tpref[node1];
              tpref[node1] = 0;
              total_tpref -= temp_p;
              // check whether sum(tpref) == 0
              for (k = 0; k < n_existing; k++)
              {
                if (tpref[k] > 0)
                {
                  break;
                }
              }
              if (k == n_existing)
              {
                total_tpref = 0;
                m_error = true;
                break;
              }

              node2 = sampleNodeLinear(n_existing, n_seednode, tpref, total_tpref, sorted_tnode);
              tpref[node1] = temp_p;
              total_tpref += temp_p;
            }
          }
        }
        else
        {
          node2 = sampleNodeLinear(n_existing, n_seednode, tpref, total_tpref, sorted_tnode);
          if (beta_loop)
          {
            node1 = sampleNodeLinear(n_existing, n_seednode, spref, total_spref, sorted_snode);
          }
          else
          {
            if (spref[node2] == total_spref)
            {
              m_error = true;
              break;
            }
            if (spref[node2] == 0)
            {
              node1 = sampleNodeLinear(n_existing, n_seednode, spref, total_spref, sorted_snode);
            }
            else
            {
              temp_p = spref[node2];
              spref[node2] = 0;
              total_spref -= temp_p;
              // check whether sum(spref) == 0
              for (k = 0; k < n_existing; k++)
              {
                if (spref[k] > 0)
                {
                  break;
                }
              }
              if (k == n_existing)
              {
                total_spref = 0;
                m_error = true;
                break;
              }

              node1 = sampleNodeLinear(n_existing, n_seednode, spref, total_spref, sorted_snode);
              spref[node2] = temp_p;
              total_spref += temp_p;
            }
          }
        }
        break;
      case 3:
        node1 = sampleNodeLinear(n_existing, n_seednode, spref, total_spref, sorted_snode);
        node2 = new_node_id;
        if (sample_recip)
        {
          node_group[node2] = sampleGroup(group_prob);
        }
        new_node_id++;
        break;
      case 4:
        node1 = new_node_id;
        new_node_id++;
        node2 = new_node_id;
        new_node_id++;
        if (sample_recip)
        {
          node_group[node1] = sampleGroup(group_prob);
          node_group[node2] = sampleGroup(group_prob);
        }
        break;
      case 5:
        node1 = node2 = new_node_id;
        if (sample_recip)
        {
          node_group[node1] = sampleGroup(group_prob);
        }
        new_node_id++;
        break;
      }
      if (m_error)
      {
        break;
      }
      // sample without replacement
      if (snode_unique && (node1 < n_existing))
      {
        total_spref -= spref[node1];
        spref[node1] = 0;
      }
      if (tnode_unique && (node2 < n_existing))
      {
        total_tpref -= tpref[node2];
        tpref[node2] = 0;
      }
      // checkDiffD(spref, total_spref);
      // checkDiffD(tpref, total_tpref);
      outs[node1] += edgeweight[new_edge_id];
      ins[node2] += edgeweight[new_edge_id];
      source_node[new_edge_id] = node1;
      target_node[new_edge_id] = node2;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      // handel reciprocal
      if (sample_recip)
      {
        if ((node1 != node2) || selfloop_recip)
        {
          p = unif_rand();
          if (p <= recip_prob(node_group[node2], node_group[node1]))
          {
            new_edge_id++;
            n_reciprocal++;
            outs[node2] += edgeweight[new_edge_id];
            ins[node1] += edgeweight[new_edge_id];
            source_node[new_edge_id] = node2;
            target_node[new_edge_id] = node1;
            scenario[new_edge_id] = 6;
          }
        }
      }
      new_edge_id++;
    }
    m[i] += n_reciprocal;
    if (m_error)
    {
      m[i] = j + n_reciprocal;
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n",
              current_scenario, i + 1, m[i]);
    }
    while (!q1.empty())
    {
      temp_node = q1.front();
      total_spref -= spref[temp_node];
      total_tpref -= tpref[temp_node];
      spref[temp_node] = calcPrefLinearD(func_type, outs[temp_node], ins[temp_node], sparams, custmSourcePrefLinear);
      tpref[temp_node] = calcPrefLinearD(func_type, outs[temp_node], ins[temp_node], tparams, custmTargetPrefLinear);
      total_spref += spref[temp_node];
      total_tpref += tpref[temp_node];
      q1.pop();
    }
    // checkDiffD(spref, total_spref);
    // checkDiffD(tpref, total_tpref);
  }
  // PutRNGstate();

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = source_node;
  ret["node_vec2"] = target_node;
  ret["outs"] = outs;
  ret["ins"] = ins;
  ret["scenario"] = scenario;
  ret["nodegroup"] = node_group;
  ret["spref"] = spref_vec;
  ret["tpref"] = tpref_vec;
  return ret;
}
