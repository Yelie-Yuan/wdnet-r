#include<iostream>
#include<queue>
#include<R.h>
#include "funcPtrD.h"
#include<Rcpp.h>

using namespace std;
funcPtrD sourcePrefFuncCppNaive;
funcPtrD targetPrefFuncCppNaive;

/**
 * Default source preference function.
 *
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param sparams Parameters passed to the source preference function.
 *  
 * @return Source preference of a node.
 */
double sourcePrefFuncDefaultNaive(double outs, double ins, Rcpp::NumericVector sparams) {
  return sparams[0] * pow(outs, sparams[1]) + 
    sparams[2] * pow(ins, sparams[3]) + sparams[4];
}

/**
 *  Default target preference function.
 * 
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param tparams Parameters passed to the target preference function.
 *  
 * @return Target preference of a node.
 */
double targetPrefFuncDefaultNaive(double outs, double ins, Rcpp::NumericVector tparams) {
  return tparams[0] * pow(outs, tparams[1]) + 
    tparams[2] * pow(ins, tparams[3]) + tparams[4];
}

/**
 *  Calculate node source preference.
 * 
 * @param func_type Default or customized preference function.
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param sparams Parameters passed to the default source preference function.
 * @param sourcePrefFuncCppNaive Pointer of the customized source preference function.
 *  
 * @return Node source preference.
 */
double calcSourcePrefNaive(int func_type, 
                          double outs,
                          double ins,
                          Rcpp::NumericVector sparams, 
                          funcPtrD sourcePrefFuncCppNaive) {
  if (func_type == 1) {
    return sourcePrefFuncDefaultNaive(outs, ins, sparams);
  }
  else {
    return sourcePrefFuncCppNaive(outs, ins);
  }
}

/**
 *  Calculate node target preference.
 * 
 * @param func_type Default or customized preference function.
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param tparams Parameters passed to the default target preference function.
 * @param targetPrefFuncCppNaive Pointer of the customized target preference function.
 *  
 * @return Node target preference.
 */
double calcTargetPrefNaive(int func_type, 
                          double outs,
                          double ins,
                          Rcpp::NumericVector tparams, 
                          funcPtrD targetPrefFuncCppNaive) {
  if (func_type == 1) {
    return targetPrefFuncDefaultNaive(outs, ins, tparams);
  }
  else {
    return targetPrefFuncCppNaive(outs, ins);
  }
}

/**
 *  Sample a source/target node.
 * 
 * @param n_existing Number of existing nodes.
 * @param pref Sequence of node source/target preference.
 * @param total_pref Total source/target preference of existing nodes.
 *  
 * @return Sampled source/target node.
 */
int sampleNodeDNaive(int n_existing, Rcpp::NumericVector pref, double total_pref) {
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
    Rprintf("Numerical error! Returning the last node (node %d) as the sampled node.\n", n_existing);
    // i = n_existing;
  }
  return i - 1;
}

/**
 *  Sample a node group.
 * 
 * @param group_prob Probability weights for sampling the group of new nodes.
 *  
 * @return Sampled group for the new node.
 */
int sampleGroupNaive(Rcpp::NumericVector group_prob) {
  double g = 0;
  int i = 0;
  while ((g == 0) || (g == 1)) {
    g = unif_rand();
  }
  while (g > 0) {
    g -= group_prob[i];
    i++;
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
//' @param source_pref Sequence of node source preference.
//' @param target_pref Sequence of node target preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List rpanet_naive_directed_cpp(
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
    Rcpp::NumericVector source_pref, 
    Rcpp::NumericVector target_pref, 
    Rcpp::List control) {
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  bool source_first = scenario_ctl["source.first"];
  Rcpp::List newedge_ctl = control["newedge"];
  // bool node_unique = ! newedge_ctl["node.replace"];
  bool snode_unique = ! newedge_ctl["snode.replace"];
  bool tnode_unique = ! newedge_ctl["tnode.replace"];
  Rcpp::List reciprocal_ctl = control["reciprocal"];
  bool selfloop_recip = reciprocal_ctl["selfloop.recip"];
  Rcpp::NumericVector group_prob = reciprocal_ctl["group.prob"];
  Rcpp::NumericMatrix recip_prob = reciprocal_ctl["recip.prob"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector sparams(5);
  Rcpp::NumericVector tparams(5);
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type) {
  case 1: 
    sparams = preference_ctl["sparams"];
    tparams = preference_ctl["tparams"];
    break;
  case 2: {
      SEXP source_pref_func_ptr = preference_ctl["spref.pointer"];
      sourcePrefFuncCppNaive = *Rcpp::XPtr<funcPtrD>(source_pref_func_ptr);
      SEXP target_pref_func_ptr = preference_ctl["tpref.pointer"];
      targetPrefFuncCppNaive = *Rcpp::XPtr<funcPtrD>(target_pref_func_ptr);
      break;
    }
  }

  double u, p, temp_p;
  bool m_error;
  int i, j, k, n_existing, current_scenario;
  int node1, node2, temp_node;
  double total_source_pref = 0, total_target_pref = 0;
  queue<int> q1;
  for (int i = 0; i < new_node_id; i++) {
    source_pref[i] = calcSourcePrefNaive(func_type, outs[i], ins[i], sparams, sourcePrefFuncCppNaive);
    target_pref[i] = calcTargetPrefNaive(func_type, outs[i], ins[i], tparams, targetPrefFuncCppNaive);
    total_source_pref += source_pref[i];
    total_target_pref += target_pref[i];
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
      if (snode_unique) {
        if ((current_scenario == 2) || (current_scenario == 3)) {
          for (k = 0; k < n_existing; k++) {
            if (source_pref[k] > 0) {
              break;
            }
          }
          if (k == n_existing) {
            total_source_pref = 0;
            m_error = true;
            break;
          }
        }
      }
      if (tnode_unique) {
        if ((current_scenario == 1) || (current_scenario == 2)) {
          for (k = 0; k < n_existing; k++) {
            if (target_pref[k] > 0) {
              break;
            }
          }
          if (k == n_existing) {
            total_target_pref = 0;
            m_error = true;
            break;
          }
        }
      }
      switch (current_scenario) {
        case 1:
          node1 = new_node_id;
          if (sample_recip) {
            node_group[node1] = sampleGroupNaive(group_prob);
          }
          new_node_id++;
          node2 = sampleNodeDNaive(n_existing, target_pref, total_target_pref);
          break;
        case 2:
          if (source_first) {
            node1 = sampleNodeDNaive(n_existing, source_pref, total_source_pref);
            if (beta_loop) {
              node2 = sampleNodeDNaive(n_existing, target_pref, total_target_pref);
            }
            else {
              if (target_pref[node1] == total_target_pref) {
                m_error = true;
                break;
              }
              if (target_pref[node1] == 0) {
                node2 = sampleNodeDNaive(n_existing, target_pref, total_target_pref);
              }
              else {
                temp_p = target_pref[node1];
                target_pref[node1] = 0;
                total_target_pref -= temp_p;
                // check whether sum(target_pref) == 0
                for (k = 0; k < n_existing; k++) {
                  if (target_pref[k] > 0) {
                    break;
                  }
                }
                if (k == n_existing) {
                  total_target_pref = 0;
                  m_error = true;
                  break;
                }

                node2 = sampleNodeDNaive(n_existing, target_pref, total_target_pref);
                target_pref[node1] = temp_p;
                total_target_pref += temp_p;
              }
            }
          }
          else {
            node2 = sampleNodeDNaive(n_existing, target_pref, total_target_pref);
            if (beta_loop) {
              node1 = sampleNodeDNaive(n_existing, source_pref, total_source_pref);
            }
            else {
              if (source_pref[node2] == total_source_pref) {
                m_error = true;
                break;
              }
              if (source_pref[node2] == 0) {
                node1 = sampleNodeDNaive(n_existing, source_pref, total_source_pref);
              }
              else {
                temp_p = source_pref[node2];
                source_pref[node2] = 0;
                total_source_pref -= temp_p;
                // check whether sum(source_pref) == 0
                for (k = 0; k < n_existing; k++) {
                  if (source_pref[k] > 0) {
                    break;
                  }
                }
                if (k == n_existing) {
                  total_source_pref = 0;
                  m_error = true;
                  break;
                }

                node1 = sampleNodeDNaive(n_existing, source_pref, total_source_pref);
                source_pref[node2] = temp_p;
                total_source_pref += temp_p;
              }
            }
          }
          break;
        case 3:
          node1 = sampleNodeDNaive(n_existing, source_pref, total_source_pref);
          node2 = new_node_id;
          if (sample_recip) {
            node_group[node2] = sampleGroupNaive(group_prob);
          }
          new_node_id++;
          break;
        case 4:
          node1 = new_node_id;
          new_node_id++;
          node2 = new_node_id;
          new_node_id++;
          if (sample_recip) {
            node_group[node1] = sampleGroupNaive(group_prob);
            node_group[node2] = sampleGroupNaive(group_prob);
          }
          break;
        case 5:
          node1 = node2 = new_node_id;
          if (sample_recip) {
            node_group[node1] = sampleGroupNaive(group_prob);
          }
          new_node_id++;
          break;
      }
      if (m_error) {
        break;
      }
      // sample without replacement
      if (snode_unique && (node1 < n_existing)) {
        total_source_pref -= source_pref[node1];
        source_pref[node1] = 0;
      }
      if (tnode_unique && (node2 < n_existing)) {
        total_target_pref -= target_pref[node2];
        target_pref[node2] = 0;
      }
      // checkDiffD(source_pref, total_source_pref);
      // checkDiffD(target_pref, total_target_pref);
      outs[node1] += edgeweight[new_edge_id];
      ins[node2] += edgeweight[new_edge_id];
      source_node[new_edge_id] = node1;
      target_node[new_edge_id] = node2;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      // handel reciprocal
      if (sample_recip) {
        if ((node1 != node2) || selfloop_recip) {
          p = unif_rand();
          if (p <= recip_prob(node_group[node2], node_group[node1])) {
            new_edge_id++;
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
    if (m_error) {
      m[i] = j;
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n", current_scenario, i + 1, j);
    }
    while(! q1.empty()) {
      temp_node = q1.front();
      total_source_pref -= source_pref[temp_node];
      total_target_pref -= target_pref[temp_node];
      source_pref[temp_node] = calcSourcePrefNaive(func_type, outs[temp_node], ins[temp_node], sparams, sourcePrefFuncCppNaive);
      target_pref[temp_node] = calcTargetPrefNaive(func_type, outs[temp_node], ins[temp_node], tparams, targetPrefFuncCppNaive);
      total_source_pref += source_pref[temp_node];
      total_target_pref += target_pref[temp_node];
      q1.pop();
    }
    // checkDiffD(source_pref, total_source_pref);
    // checkDiffD(target_pref, total_target_pref);
  }
  PutRNGstate();

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = source_node;
  ret["node_vec2"] = target_node;
  ret["outstrength"] = outs;
  ret["instrength"] = ins;
  ret["scenario"] = scenario;
  ret["nodegroup"] = node_group;
  ret["source_pref"] = source_pref;
  ret["target_pref"] = target_pref;
  return ret;
}