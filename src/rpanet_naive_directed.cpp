#include<iostream>
#include<queue>
#include<deque>
#include<algorithm>
#include<R.h>
#include "funcPtrD.h"
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

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
 * @param sparams Parameters passed to the source preference function.
 * @param sourcePrefFuncCppNaive Source preference function.
 *  
 * @return Node source preference.
 */
double calSourcePrefNaive(int func_type, 
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
 * @param tparams Parameters passed to the source preference function.
 * @param targetPrefFuncCppNaive Source preference function.
 *  
 * @return Node target preference.
 */
double calTargetPrefNaive(int func_type, 
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
 * @param pref Sequence of node source/target preference.
 * @param total_pref Total source/target preference of existing nodes.
 * @param qm Nodes to be excluded from the sampling process.
 *  
 * @return Sampled source/target node.
 */
int sampleNodeDNaive(Rcpp::NumericVector pref, double total_pref, deque<int> &qm) {
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
  bool node_unique = ! newedge_ctl["node.replace"];
  bool snode_unique = ! newedge_ctl["snode.replace"];
  bool tnode_unique = ! newedge_ctl["tnode.replace"];
  Rcpp::List reciprocal_ctl = control["reciprocal"];
  bool selfloop_recip = reciprocal_ctl["selfloop.recip"];
  Rcpp::NumericVector group_prob = reciprocal_ctl["group.prob"];
  Rcpp::NumericMatrix recip_prob = reciprocal_ctl["recip.prob"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector sparams(5);
  Rcpp::NumericVector tparams(5);
  funcPtrD sourcePrefFuncCppNaive;
  funcPtrD targetPrefFuncCppNaive;
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type) {
  case 1: 
    sparams = preference_ctl["sparams"];
    tparams = preference_ctl["tparams"];
    break;
  case 2: {
      SEXP source_pref_func_ptr = preference_ctl["spref.pointer"];
      Rcpp::XPtr<funcPtrD> xpfunSource(source_pref_func_ptr);
      sourcePrefFuncCppNaive = *xpfunSource;
      SEXP target_pref_func_ptr = preference_ctl["tpref.pointer"];
      Rcpp::XPtr<funcPtrD> xpfunTarget(target_pref_func_ptr);
      targetPrefFuncCppNaive = *xpfunTarget;
      break;
    }
  }

  double u, p;
  bool check_unique = node_unique || snode_unique || tnode_unique;
  bool m_error;
  int i, j, ks, kt, n_existing, current_scenario;
  int node1, node2, temp_node;
  double total_source_pref = 0, total_target_pref = 0;
  queue<int> q1;
  deque<int> qm_source, qm_target;
  for (int i = 0; i < new_node_id; i++) {
    source_pref[i] = calSourcePrefNaive(func_type, outs[i], ins[i], sparams, sourcePrefFuncCppNaive);
    target_pref[i] = calTargetPrefNaive(func_type, outs[i], ins[i], tparams, targetPrefFuncCppNaive);
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
      ks = qm_source.size();
      kt = qm_target.size();
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
      if (check_unique) {
        switch (current_scenario) {
          case 1:
            if (kt + 1 > n_existing) {
              m_error = true;
            }
            break;
          case 2:
            if (node_unique) {
              if (ks + 2 - int(beta_loop) > n_existing) {
                m_error = true;
              }
            }
            else {
              if (snode_unique) {
                if (ks + 1 > n_existing) {
                  m_error = true;
                }
              }
              if (tnode_unique) {
                if (kt + 1 > n_existing) {
                  m_error = true;
                }
              }
            }
            break;
          case 3:
            if (ks + 1 > n_existing) {
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
          if (sample_recip) {
            node_group[node1] = sampleGroupNaive(group_prob);
          }
          new_node_id++;
          node2 = sampleNodeDNaive(target_pref, total_target_pref, qm_target);
          break;
        case 2:
          if (source_first) {
            node1 = sampleNodeDNaive(source_pref, total_source_pref, qm_source);
            if (beta_loop) {
              node2 = sampleNodeDNaive(target_pref, total_target_pref, qm_target);
            }
            else {
              if (find(qm_target.begin(), qm_target.end(), node1) != qm_target.end()) {
                node2 = sampleNodeDNaive(target_pref, total_target_pref, qm_target);
              }
              else {
                if (kt + 2 > n_existing) {
                  m_error = true;
                  break;
                }
                qm_target.push_back(node1);
                node2 = sampleNodeDNaive(target_pref, total_target_pref, qm_target);
                qm_target.pop_back();
              }
            }
          }
          else {
            node2 = sampleNodeDNaive(target_pref, total_target_pref, qm_target);
            if (beta_loop) {
              node1 = sampleNodeDNaive(source_pref, total_source_pref, qm_source);
            }
            else {
              if (find(qm_source.begin(), qm_source.end(), node2) != qm_source.end()) {
                node1 = sampleNodeDNaive(source_pref, total_source_pref, qm_source);
              }
              else {
                if (ks + 2 > n_existing) {
                  m_error = true;
                  break;
                }
                qm_source.push_back(node2);
                node1 = sampleNodeDNaive(source_pref, total_source_pref, qm_source);
                qm_source.pop_back();
              }
            }
          }
          break;
        case 3:
          node1 = sampleNodeDNaive(source_pref, total_source_pref, qm_source);
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
      // handle duplicate nodes
      if (node_unique) {
        if (node1 < n_existing) {
          qm_source.push_back(node1);
          qm_target.push_back(node1);
        }
        if ((node2 < n_existing) && (node1 != node2)) {
          qm_source.push_back(node2);
          qm_target.push_back(node2);
        }
      }
      else {
        if (snode_unique && (node1 < n_existing)) {
          qm_source.push_back(node1);
        }
        if (tnode_unique && (node2 < n_existing)) {
          qm_target.push_back(node2);
        }
      }
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
      Rprintf("Unique nodes exhausted at step %u. Set the value of m at current step to %u.\n", i + 1, j);
    }
    while(! q1.empty()) {
      temp_node = q1.front();
      total_source_pref -= source_pref[temp_node];
      total_target_pref -= target_pref[temp_node];
      source_pref[temp_node] = calSourcePrefNaive(func_type, outs[temp_node], ins[temp_node], sparams, sourcePrefFuncCppNaive);
      target_pref[temp_node] = calTargetPrefNaive(func_type, outs[temp_node], ins[temp_node], tparams, targetPrefFuncCppNaive);
      total_source_pref += source_pref[temp_node];
      total_target_pref += target_pref[temp_node];
      q1.pop();
    }
    qm_source.clear();
    qm_target.clear();
  }
  PutRNGstate();
  // check total preference = sum of node preference
  // Rprintf("Total source pref %f.\n", total_source_pref);
  // Rprintf("Total target pref %f.\n", total_target_pref);
  // for (i = 0; i < new_node_id; i++) {
  //   total_source_pref -= source_pref[i];
  //   total_target_pref -= target_pref[i];
  // }
  // Rprintf("Total source pref %f.\n", total_source_pref * pow(10, 10));
  // Rprintf("Total target pref %f.\n", total_target_pref * pow(10, 10));

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