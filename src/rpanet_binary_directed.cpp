#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
#include<deque>
#include<algorithm>
using namespace std;

/**
 * Node structure.
 * id: node id
 * outs, ins: out- and in-strength
 * sourcep: preference of being chosen as a source node
 * targetp: preference of being chosed as a target node
 * total_sourcep: sum of sourcep of current node and its children
 * total_targetp: sum of targetp of current node and its children
 * *left, *right, *parent: pointers to its left, right and parent
 */
struct node {
  int id, group;
  double outs, ins;
  double sourcep, targetp, total_sourcep, total_targetp;
  node *left, *right, *parent;
};


/**
 * Source preference function.
 *
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param source_params Parameters passed to the source preference function.
 * 
 * @return Source preference of a node.
 * 
 */
double sourcePreferenceFunc(double outs, double ins, double *source_params) {
  return source_params[0] * pow(outs, source_params[1]) + 
    source_params[2] * pow(ins, source_params[3]) + source_params[4];
}

/**
 * Target preference function.
 *
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param target_params Parameters passed to the target preference function.
 * 
 * @return Target preference of a node.
 * 
 */
double targetPreferenceFunc(double outs, double ins, double *target_params) {
  return target_params[0] * pow(outs, target_params[1]) + 
    target_params[2] * pow(ins, target_params[3]) + target_params[4];
}

/**
 * Update total source preference from current node to root.
 *
 * @param current_node The current node.
 * @param increment Value to be added to the total source preference.
 * 
 */
void addSourceIncrement(node *current_node, double increment) {
  current_node->total_sourcep += increment;
  while(current_node->id > 0) {
    return addSourceIncrement(current_node->parent, increment);
  }
}

/**
 * Update total target preference from current node to root.
 *
 * @param current_node The current node.
 * @param increment Value to be added to the total target preference.
 * 
 */
void addTargetIncrement(node *current_node, double increment) {
  current_node->total_targetp += increment;
  while(current_node->id > 0) {
    return addTargetIncrement(current_node->parent, increment);
  }
}

/**
 * Update node strength, preference and total preference from the sampled node to root.
 *
 * @param temp_node The sampled node.
 * @param source_params Parameters passed to the source preference function.
 * @param target_params Parameters passed to the target preference function.
 * 
 */
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

/**
 * Create a new node.
 *
 * @param id Node ID.
 * 
 * @return The new node.
 */
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

/**
 * Insert a new node to the tree.
 *
 * @param q Sequence of nodes that have less than 2 children.
 * @param new_node_id New node ID.
 * 
 * @return The new node.
 */
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

/**
 * Find a source node with a given cutoff point w.
 *
 * @param root Root node of the tree.
 * @param w A cutoff point.
 * 
 * @return Sampled source/target node.
 */
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

/**
 * Find a target node with a given cutoff point w.
 *
 * @param root Root node of the tree.
 * @param w A cutoff point.
 * 
 * @return Sampled source/target node.
 */
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

/**
 * Sample a source/target node from the tree.
 *
 * @param root Root node of the tree.
 * @param type Represent source node or target node.
 * @param qm Nodes to be excluded from the sampling process.
 * 
 * @return Sampled source/target node.
 */
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

/**
 * Sample a node group.
 *
 * @param group_dist Probability weights for sampling the group of new nodes.
 * 
 * @return Sampled group for the new node.
 */
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
  /**
   * Preferential attachment algorithm.
   *
   * @param nstep_ptr Number of steps.
   * @param m Number of new edges in each step.
   * @param new_node_id_ptr New node ID.
   * @param new_edge_id_ptr New edge ID.
   * @param source_node Sequence of source nodes.
   * @param target_node Sequence of target nodes.
   * @param outs Sequence of out-strength.
   * @param ins Sequence of in-strength.
   * @param edgeweight Weight of existing and new edges.
   * @param scenario Scenario of existing and new edges.
   * @param alpha_ptr Probability of alpha acenario.
   * @param beta_ptr Probability of beta acenario.
   * @param gamma_ptr Probability of gamma acenario.
   * @param xi_ptr Probability of xi acenario.
   * @param beta_loop_ptr Whether self loops are allowed under beta scenario.
   * @param source_first_ptr Logical, wheter the source node is sampled prior to the target
   *   node when adding beta scenario edges.
   * @param node_unique_ptr Logical, whether the nodes in the same step should bedifferent from
   *   each other. Defined for undirected and directed networks. For directed networks, when 
   *   node.unique is TRUE, sampled source and target nodes in the same step are all different 
   *   from each other, beta.loop will be FALSE, snode.unique and tnode.unique will 
   *   be TRUE.
   * @param snode_unique_ptr Logical, whether the source nodes in the same step should 
   *   be sampled different from each other. Defined for directed networks.
   * @param tnode_unique_ptr Logical, whether the target nodes in the same step should 
   *   be sampled different from each other. Defined for directed networks.
   * @param source_params Parameters of the source preference function for directed networks. 
   *   Probability of choosing an existing node as the source node is proportional 
   *   to sparams[1] * out-strength^sparams[2] + sparams[3] * in-strength^sparams[4] + sparams[5].
   * @param target_params  Parameters of the target preference function for directed networks. 
   *   Probability of choosing an existing node as the source node is proportional to 
   *   tparams[1] * out-strength^tparams[2] + tparams[3] * in-strength^tparams[4] + tparams[5].
   * @param sample_recip_ptr Logical, whether reciprocal edges will be added.
   * @param selfloop_recip_ptr Logical, whether reciprocal of self loops are allowed.
   * @param group_dist Probability weights for sampling the group of new nodes. Defined for 
   *   directed networks. Groups are 1:length(group_dist) in R, and 0:(length(group_dist) - 1)
   *   in c. length(group_dist) must equal to the square root of length(recip).
   * @param recip The probability of adding a reciprocal edge after a new edge is introduced.
   *   Vectorized from the matrix recip.prob.
   * @param node_group Sequence of node group.
   * @param ngroup_ptr Number of groups.
   * @param source_pref Sequence of node source preference.
   * @param target_pref Sequence of node target preference.
   * 
   */
  void rpanet_binary_directed_cpp(
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
      int *sample_recip_ptr, int *selfloop_recip_ptr,
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
      selfloop_recip = *selfloop_recip_ptr,
      check_unique = node_unique | snode_unique | tnode_unique;
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
            node1 = insertNode2(q, new_node_id);
            if (sample_recip) {
              node1->group = sampleGroup(group_dist);
            }
            new_node_id++;
            node2 = sampleNode2(root, 't', qm_target);
            break;
          case 2:
            if (source_first) {
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
            }
            else {
              node2 = sampleNode2(root, 't', qm_target);
              if (beta_loop) {
                node1 = sampleNode2(root, 's', qm_source);
              }
              else {
                if (find(qm_source.begin(), qm_source.end(), node2) != qm_source.end()) {
                  node1 = sampleNode2(root, 's', qm_source);
                }
                else {
                  if (ks + 2 > n_existing) {
                    m_error = true;
                    break;
                  }
                  qm_source.push_back(node2);
                  node1 = sampleNode2(root, 's', qm_source);
                  qm_source.pop_back();
                }
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
        if (node_unique) {
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
          if (snode_unique & (node1->id < n_existing)) {
            qm_source.push_back(node1);
          }
          if (tnode_unique & (node2->id < n_existing)) {
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
        // handle reciprocal
        if (sample_recip) {
          if ((node1->id != node2->id) | selfloop_recip) {
            p = unif_rand();
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