#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
#include<bits/stdc++.h>
using namespace std;

// 1. user defined preference functions (sourcePreferenceFunc and 
//    targetPreferenceFunc)
// 2. w can not equal to 1 in function SampleNode, otherwise findNode returns error 
//    because of numeric precision
// 3. add a parameter m to control number of new edges per step

// node structure
// id: node id
// outs, ins: out- and in-strength
// sourcep: preference of being chosen as a source node
// targetp: preference of being chosed as a target node
// total_sourcep: sum of sourcep of current node and its children
// total_targetp: sum of targetp of current node and its children
// *left, *right, *parent: pointers to its left, right and parent
struct node {
  int id, group;
  double outs, ins;
  double sourcep, targetp, total_sourcep, total_targetp;
  node *left, *right, *parent;
};

// preference functions
double sourcePreferenceFunc(double outs, double ins, double *source_params) {
  return source_params[0] * pow(outs, source_params[1]) + 
    source_params[2] * pow(ins, source_params[3]) + source_params[4];
}
double targetPreferenceFunc(double outs, double ins, double *target_params) {
  return target_params[0] * pow(outs, target_params[1]) + 
    target_params[2] * pow(ins, target_params[3]) + target_params[4];
}

// update total source preference from current node to root
void addSourceIncrement(node *current_node, double increment) {
  current_node->total_sourcep += increment;
  while(current_node->id > 0) {
    return addSourceIncrement(current_node->parent, increment);
  }
}
// update total source preference from current node to root
void addTargetIncrement(node *current_node, double increment) {
  current_node->total_targetp += increment;
  while(current_node->id > 0) {
    return addTargetIncrement(current_node->parent, increment);
  }
}

// update strength, preference and total preference from the sampled node to root
void updatePreference2(node *temp_node, 
    double *source_params, double *target_params) {
  double tp = temp_node->sourcep;
  temp_node->sourcep = sourcePreferenceFunc(temp_node->outs, temp_node->ins, 
    source_params);
  if (temp_node->sourcep != tp) {
    addSourceIncrement(temp_node, temp_node->sourcep - tp);
  }
  tp = temp_node->targetp;
  temp_node->targetp = targetPreferenceFunc(temp_node->outs, temp_node->ins, 
    target_params);
  if (temp_node->targetp != tp) {
    addTargetIncrement(temp_node, temp_node->targetp - tp);
  }
}

// create a new node
node *createNode2(int id) {
  node *new_node = new node();
  new_node->id = id;
  new_node->group = -1;
  new_node->outs = new_node->ins = 0;
  new_node->sourcep = new_node->total_sourcep = 0;
  new_node->targetp = new_node->total_targetp = 0;
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
}

// insert a new node to the tree
node *insertNode2(queue<node*> &q, int new_node_id) {
  node *new_node = createNode2(new_node_id);
  node *temp_node = q.front();
  if(temp_node->left == NULL) {
    temp_node->left = new_node;
  }
  else if (temp_node->right == NULL) {
    temp_node->right = new_node;
    q.pop();
  }
  new_node->parent = temp_node;
  q.push(new_node);
  return new_node;
}

// find source node with a given critical point w
node *findSourceNode(node *root, double w) {
  w -= root->sourcep;
  if (w <= 0) {
    return root;
  } 
  else {
    if (w > root->left->total_sourcep) {
      return findSourceNode(root->right, w - root->left->total_sourcep);
    }
    else {
      return findSourceNode(root->left, w);
    }
  }
}
// find target node with a given critical point w
node *findTargetNode(node *root, double w) {
  w -= root->targetp;
  if (w <= 0) {
    return root;
  }
  else {
    if (w > root->left->total_targetp) {
      return findTargetNode(root->right, w - root->left->total_targetp);
    }
    else {
      return findTargetNode(root->left, w);
    }
  }
}

// sample a node from the tree
node* sampleNode2(node *root, char type, deque<node*> &qm) {
  double w;
  node *temp_node;
  while (true) {
    w = 1;
    while (w == 1) {
      w = unif_rand();
    }
    if (type == 's') {
      w *= root->total_sourcep;
      temp_node = findSourceNode(root, w);
    }
    else {
      w *= root->total_targetp;
      temp_node = findTargetNode(root, w);
    }
    if (find(qm.begin(), qm.end(), temp_node) == qm.end()) {
      // if temp_node not in qm
      return temp_node;
    }
  }
}

// sample a node group
int sampleGroup(double *group_dist) {
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
  void rpanet_binary_directed_cpp(
      int *nstep_ptr, int *m, 
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *source_node, int *target_node, 
      double *outs, double *ins, 
      double *edgeweight, int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr, int *m_unique_ptr,
      int *m_source_unique_ptr, int *m_target_unique_ptr,
      double *source_params, double *target_params, 
      int *sample_recip_ptr,
      double *group_dist, double *recip, 
      int *node_group, int *ngroup_ptr, 
      double *source_pref, double *target_pref) {
    double u, p;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr,
      new_edge_id = *new_edge_id_ptr, ngroup = *ngroup_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr, m_unique = *m_unique_ptr, 
      m_source_unique = *m_source_unique_ptr,
      m_target_unique = *m_target_unique_ptr, 
      m_error, sample_recip = *sample_recip_ptr, 
      check_unique = m_unique | m_source_unique | m_target_unique;
    int i, j, ks, kt, n_existing, current_scenario;
    node *node1, *node2;
    // initialize a tree from the seed graph
    node *root = createNode2(0);
    root->outs = outs[0];
    root->ins = ins[0];
    root->group = node_group[0];
    updatePreference2(root, source_params, target_params);
    queue<node*> q, q1;
    deque<node*> qm_source, qm_target;
    q.push(root);
    for (int i = 1; i < new_node_id; i++) {
      node1 = insertNode2(q, i);
      node1->outs = outs[i];
      node1->ins = ins[i];
      node1->group = node_group[i];
      updatePreference2(node1, source_params, target_params);
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
              if (m_unique) {
                if (ks + 2 - int(beta_loop) > n_existing) {
                  m_error = true;
                }
              }
              else {
                if (m_source_unique) {
                  if (ks + 1 > n_existing) {
                    m_error = true;
                  }
                }
                if (m_target_unique) {
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
            node1 = insertNode2(q, new_node_id);
            if (sample_recip) {
              node1->group = sampleGroup(group_dist);
            }
            new_node_id++;
            node2 = sampleNode2(root, 't', qm_target);
            break;
          case 2:
            node1 = sampleNode2(root, 's', qm_source);
            if (beta_loop) {
              node2 = sampleNode2(root, 't', qm_target);
            }
            else {
              if (find(qm_target.begin(), qm_target.end(), node1) != qm_target.end()) {
                node2 = sampleNode2(root, 't', qm_target);
              }
              else {
                if (kt + 2 > n_existing) {
                  m_error = true;
                  break;
                }
                qm_target.push_back(node1);
                node2 = sampleNode2(root, 't', qm_target);
                qm_target.pop_back();
              }
            }
            break;
          case 3:
            node1 = sampleNode2(root, 's', qm_source);
            node2 = insertNode2(q, new_node_id);
            if (sample_recip) {
              node2->group = sampleGroup(group_dist);
            }
            new_node_id++;
            break;
          case 4:
            node1 = insertNode2(q, new_node_id);
            new_node_id++;
            node2 = insertNode2(q, new_node_id);
            new_node_id++;
            if (sample_recip) {
              node1->group = sampleGroup(group_dist);
              node2->group = sampleGroup(group_dist);
            }
            break;
          case 5:
            node1 = node2 = insertNode2(q, new_node_id);
            if (sample_recip) {
              node1->group = sampleGroup(group_dist);
            }
            new_node_id++;
            break;
        }
        if (m_error) {
          break;
        }
        // handle duplicate nodes
        if (m_unique) {
          if (node1->id < n_existing) {
            qm_source.push_back(node1);
            qm_target.push_back(node1);
          }
          if ((node2->id < n_existing) & (node1 != node2)) {
            qm_source.push_back(node2);
            qm_target.push_back(node2);
          }
        }
        else {
          if (m_source_unique & (node1->id < n_existing)) {
            qm_source.push_back(node1);
          }
          if (m_target_unique & (node2->id < n_existing)) {
            qm_target.push_back(node2);
          }
        }
        node1->outs += edgeweight[new_edge_id];
        node2->ins += edgeweight[new_edge_id];
        source_node[new_edge_id] = node1->id;
        target_node[new_edge_id] = node2->id;
        scenario[new_edge_id] = current_scenario;
        q1.push(node1);
        q1.push(node2);
        // handel reciprocal
        if (sample_recip) {
          p = unif_rand();
          if (scenario[new_edge_id] != 5) {
            if (p <= recip[node2->group * ngroup + node1->group]) {
              new_edge_id++;
              node2->outs += edgeweight[new_edge_id];
              node1->ins += edgeweight[new_edge_id];
              source_node[new_edge_id] = node2->id;
              target_node[new_edge_id] = node1->id;
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
        updatePreference2(q1.front(), source_params, target_params);
        q1.pop();
      }
      qm_source.clear();
      qm_target.clear();
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
    // free memory (queue)
    queue<node*>().swap(q);
    queue<node*>().swap(q1);
    // save strength and preference
    q.push(root);
    node *temp_node;
    j = 0;
    while (! q.empty())
    {
      temp_node = q.front();
      q.pop();
      if (temp_node->left != NULL) {
        q.push(temp_node->left);
      }
      if (temp_node->right != NULL) {
        q.push(temp_node->right);
      }
      outs[j] = temp_node->outs;
      ins[j] = temp_node->ins;
      node_group[j] = temp_node->group;
      source_pref[j] = temp_node->sourcep;
      target_pref[j] = temp_node->targetp;
      // free memory (node and tree)
      delete(temp_node);
      j++;
    }
    // free memory (queue)
    queue<node*>().swap(q);
  }
}