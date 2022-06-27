#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
#include<deque>
#include<algorithm>
using namespace std;

/**
 * Preference function.
 *
 * @param strength Node strength.
 * @param params Parameters passed to the preference function.
 * 
 * @return Preference of a node.
 * 
 */
double preferenceFuncNaive(double strength, double *params) {
  return pow(strength, params[0]) + params[1];
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
int sampleNodeNaiveUnd(double *pref, double total_pref, deque<int> &qm) {
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

extern "C" {
  /**
   * Preferential attachment algorithm.
   *
   * @param nstep_ptr Number of steps.
   * @param m Number of new edges in each step.
   * @param new_node_id_ptr New node ID.
   * @param new_edge_id_ptr New edge ID.
   * @param node_vec1 Sequence of nodes in the first column of edgelist.
   * @param node_vec2 Sequence of nodes in the second column of edgelist.
   * @param strength Sequence of node strength.
   * @param edgeweight Weight of existing and new edges.
   * @param scenario Scenario of existing and new edges.
   * @param alpha_ptr Probability of alpha acenario.
   * @param beta_ptr Probability of beta acenario.
   * @param gamma_ptr Probability of gamma acenario.
   * @param xi_ptr Probability of xi acenario.
   * @param beta_loop_ptr Whether self loops are allowed under beta scenario.
   * @param node_unique_ptr Logical, whether the nodes in the same step should bedifferent from
   *   each other.
   * @param params Parameters of the preference function for undirected networks. 
   *   Probability of choosing an existing node is proportional to strength^param[1] + param[2].
   * @param pref Sequence of node preference.
   * 
   */
  void rpanet_naive_undirected_cpp(
      int *nstep_ptr, int *m,
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *node_vec1, int *node_vec2, 
      double *strength, double *edgeweight, 
      int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr, int *node_unique_ptr,
      double *params, double *pref) {
    double u, total_pref = 0;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, 
      new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, 
      gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr, node_unique = *node_unique_ptr,
      m_error;
    int i, j, k, n_existing, current_scenario;
    int node1, node2, temp_node;
    queue<int> q1;
    deque<int> qm;
    for (i = 0; i < new_node_id; i++) {
      pref[i] = preferenceFuncNaive(strength[i], params);
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
            node2 = sampleNodeNaiveUnd(pref, total_pref, qm);
            break;
          case 2:
            node1 = sampleNodeNaiveUnd(pref, total_pref, qm);
            if (! beta_loop) {
              qm.push_back(node1);
              node2 = sampleNodeNaiveUnd(pref, total_pref, qm);
              qm.pop_back();
            }
            else {
              node2 = sampleNodeNaiveUnd(pref, total_pref, qm);
            }
            break;
          case 3:
            node1 = sampleNodeNaiveUnd(pref, total_pref, qm);
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
        pref[temp_node] = preferenceFuncNaive(strength[temp_node], params);
        total_pref += pref[temp_node];
        q1.pop();
      }
      qm.clear();
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
    // check total preference = sum of node preference
    // Rprintf("Total pref %f.\n", total_pref);
    // for (i = 0; i < new_node_id; i++) {
    //   total_pref -= pref[i];
    // }
    // Rprintf("Total pref %f.\n", total_pref * pow(10, 10));
  }
}