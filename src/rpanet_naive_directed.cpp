#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
// #include<bits/stdc++.h>
#include<deque>
#include<algorithm>
using namespace std;

// preference functions
double sourcePreferenceFuncNaive(double outs, double ins, double *source_params) {
  return source_params[0] * pow(outs, source_params[1]) + 
    source_params[2] * pow(ins, source_params[3]) + source_params[4];
}
double targetPreferenceFuncNaive(double outs, double ins, double *target_params) {
  return target_params[0] * pow(outs, target_params[1]) + 
    target_params[2] * pow(ins, target_params[3]) + target_params[4];
}

// sample a node from the tree
int sampleNodeNaive2(double *pref, double total_pref, deque<int> &qm) {
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

// sample a node group
int sampleGroupNaive(double *group_dist) {
  double g = 0;
  int i = 0;
  while ((g == 0) | (g == 1)) {
    g = unif_rand();
  }
  while (g > 0) {
    g -= group_dist[i];
    i++;
  }
  return i - 1;
}

extern "C" {
  void rpanet_naive_directed_cpp(
      int *nstep_ptr, int *m, 
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *source_node, int *target_node, 
      double *outs, double *ins, 
      double *edgeweight, int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr, int *source_first_ptr,
      int *node_unique_ptr,
      int *snode_unique_ptr, int *tnode_unique_ptr,
      double *source_params, double *target_params, 
      int *sample_recip_ptr,
      double *group_dist, double *recip, 
      int *node_group, int *ngroup_ptr, 
      double *source_pref, double *target_pref) {
    double u, p;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr,
      new_edge_id = *new_edge_id_ptr, ngroup = *ngroup_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr, 
      source_first = *source_first_ptr,
      node_unique = *node_unique_ptr, 
      snode_unique = *snode_unique_ptr,
      tnode_unique = *tnode_unique_ptr, 
      m_error, sample_recip = *sample_recip_ptr, 
      check_unique = node_unique | snode_unique | tnode_unique;
    int i, j, ks, kt, n_existing, current_scenario;
    int node1, node2, temp_node;
    double total_source_pref = 0, total_target_pref = 0;
    queue<int> q1;
    deque<int> qm_source, qm_target;
    for (int i = 0; i < new_node_id; i++) {
      source_pref[i] = sourcePreferenceFuncNaive(outs[i], ins[i], source_params);
      target_pref[i] = targetPreferenceFuncNaive(outs[i], ins[i], target_params);
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
              node_group[node1] = sampleGroupNaive(group_dist);
            }
            new_node_id++;
            node2 = sampleNodeNaive2(target_pref, total_target_pref, qm_target);
            break;
          case 2:
            if (source_first) {
              node1 = sampleNodeNaive2(source_pref, total_source_pref, qm_source);
              if (beta_loop) {
                node2 = sampleNodeNaive2(target_pref, total_target_pref, qm_target);
              }
              else {
                if (find(qm_target.begin(), qm_target.end(), node1) != qm_target.end()) {
                  node2 = sampleNodeNaive2(target_pref, total_target_pref, qm_target);
                }
                else {
                  if (kt + 2 > n_existing) {
                    m_error = true;
                    break;
                  }
                  qm_target.push_back(node1);
                  node2 = sampleNodeNaive2(target_pref, total_target_pref, qm_target);
                  qm_target.pop_back();
                }
              }
            }
            else {
              node2 = sampleNodeNaive2(target_pref, total_target_pref, qm_target);
              if (beta_loop) {
                node1 = sampleNodeNaive2(source_pref, total_source_pref, qm_source);
              }
              else {
                if (find(qm_source.begin(), qm_source.end(), node2) != qm_source.end()) {
                  node1 = sampleNodeNaive2(source_pref, total_source_pref, qm_source);
                }
                else {
                  if (ks + 2 > n_existing) {
                    m_error = true;
                    break;
                  }
                  qm_source.push_back(node2);
                  node1 = sampleNodeNaive2(source_pref, total_source_pref, qm_source);
                  qm_source.pop_back();
                }
              }
            }
            break;
          case 3:
            node1 = sampleNodeNaive2(source_pref, total_source_pref, qm_source);
            node2 = new_node_id;
            if (sample_recip) {
              node_group[node2] = sampleGroupNaive(group_dist);
            }
            new_node_id++;
            break;
          case 4:
            node1 = new_node_id;
            new_node_id++;
            node2 = new_node_id;
            new_node_id++;
            if (sample_recip) {
              node_group[node1] = sampleGroupNaive(group_dist);
              node_group[node2] = sampleGroupNaive(group_dist);
            }
            break;
          case 5:
            node1 = node2 = new_node_id;
            if (sample_recip) {
              node_group[node1] = sampleGroupNaive(group_dist);
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
          if ((node2 < n_existing) & (node1 != node2)) {
            qm_source.push_back(node2);
            qm_target.push_back(node2);
          }
        }
        else {
          if (snode_unique & (node1 < n_existing)) {
            qm_source.push_back(node1);
          }
          if (tnode_unique & (node2 < n_existing)) {
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
          p = unif_rand();
          if (scenario[new_edge_id] != 5) {
            if (p <= recip[node_group[node2] * ngroup + node_group[node1]]) {
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
        source_pref[temp_node] = sourcePreferenceFuncNaive(outs[temp_node], 
          ins[temp_node], source_params);
        target_pref[temp_node] = targetPreferenceFuncNaive(outs[temp_node],  
          ins[temp_node], target_params);
        total_source_pref += source_pref[temp_node];
        total_target_pref += target_pref[temp_node];
        q1.pop();
      }
      qm_source.clear();
      qm_target.clear();
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
  }
}