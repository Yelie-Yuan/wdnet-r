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
 * strength: node strength
 * p: preference of being chosen from the existing nodes
 * totalp: sum of preference of current node and its children
 * *left, *right, *parent: pointers to its left, right and parent
 */
struct node {
  int id;
  double strength;
  double p, totalp;
  node *left, *right, *parent;
};

/**
 * Preference function.
 *
 * @param strength Node strength.
 * @param params Parameters passed to the preference function.
 * 
 * @return Preference of a node.
 * 
 */
double preferenceFunc(double strength, double *params) {
  return pow(strength, params[0]) + params[1];
}

/**
 * Update total preference from current node to root.
 *
 * @param current_node The current node.
 * @param increment Value to be added to the total preference.
 * 
 */
void addIncrement(node *current_node, double increment) {
  current_node->totalp += increment;
  while(current_node->id > 0) {
    return addIncrement(current_node->parent, increment);
  }
}

/**
 * Update node preference and total preference from the sampled node to root.
 *
 * @param temp_node The sampled node.
 * @param params Parameters passed to the preference function.
 * 
 */
void updatePreference(node *temp_node, double *params) {
  double temp_p = temp_node->p;
  temp_node->p = preferenceFunc(temp_node->strength, params);
  addIncrement(temp_node, temp_node->p - temp_p);
}

/**
 * Create a new node.
 *
 * @param id Node ID.
 * 
 * @return The new node.
 */
node *createNode(int id) {
  node *new_node = new node();
  new_node->id = id;
  new_node->strength = 0;
  new_node->p = new_node->totalp = 0;
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
}

/**
 * Insert a new node to the tree.
 *
 * @param q Sequence of nodes that have less than 2 children.
 * @param id New node ID.
 * 
 * @return The new node.
 */
node *insertNode(queue<node*> &q, int id) {
  node *new_node = createNode(id);
  node *temp_node = q.front();
  // check left
  if(temp_node->left == NULL) {
    temp_node->left = new_node;
  }
  // check right
  else if (temp_node->right == NULL) {
    temp_node->right = new_node;
    q.pop();
  }
  new_node->parent = temp_node;
  q.push(new_node);
  return new_node;
}

/**
 * Find a node with a given cutoff point w.
 *
 * @param root Root node of the tree.
 * @param w A cutoff point.
 * 
 * @return Sampled node.
 */
node *findNode(node *root, double w) {
  w -= root->p;
  if (w <= 0) {
    return root;
  }
  else {
    if (w > root->left->totalp) {
      return findNode(root->right, w - root->left->totalp);
    }
    else {
      return findNode(root->left, w);
    }
  }
}

/**
 * Sample a node from the tree.
 *
 * @param root Root node of the tree.
 * @param qm Nodes to be excluded from the sampling process.
 * 
 * @return Sampled node.
 */
node *sampleNode(node *root, deque<node*> &qm) {
  double w;
  node *temp_node;
  while (true) {
    w = 1;
    while (w == 1) {
      w = unif_rand();
    }
    w *= root->totalp;
    temp_node = findNode(root, w);
    if (find(qm.begin(), qm.end(), temp_node) == qm.end()) {
      // if temp_node not in qm
      return temp_node;
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
  void rpanet_binary_undirected_cpp(
      int *nstep_ptr, int *m,
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *node_vec1, int *node_vec2, 
      double *strength, double *edgeweight, 
      int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr, int *node_unique_ptr,
      double *params, double *pref) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, 
      new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, 
      gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr, node_unique = *node_unique_ptr,
      m_error;
    int i, j, k, n_existing, current_scenario;
    node *node1, *node2;
    // initialize a tree from seed graph
    node *root = createNode(0);
    root->strength = strength[0];
    updatePreference(root, params);
    queue<node*> q, q1;
    deque<node*> qm;
    q.push(root);
    for (i = 1; i < new_node_id; i++) {
      node1 = insertNode(q, i);
      node1->strength = strength[i];
      updatePreference(node1, params);
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
            node1 = insertNode(q, new_node_id);
            new_node_id++;
            node2 = sampleNode(root, qm);
            break;
          case 2:
            node1 = sampleNode(root, qm);
            if (! beta_loop) {
              qm.push_back(node1);
              node2 = sampleNode(root, qm);
              qm.pop_back();
            }
            else {
              node2 = sampleNode(root, qm);
            }
            break;
          case 3:
            node1 = sampleNode(root, qm);
            node2 = insertNode(q, new_node_id);
            new_node_id++;
            break;
          case 4:
            node1 = insertNode(q, new_node_id);
            new_node_id++;
            node2 = insertNode(q, new_node_id);
            new_node_id++;
            break;
          case 5:
            node1 = node2 = insertNode(q, new_node_id);
            new_node_id++;
            break;
        }
        // handle duplicate nodes
        if (node_unique) {
          if (node1->id < n_existing) {
            qm.push_back(node1);
          }
          if ((node2->id < n_existing) & (node1 != node2)) {
            qm.push_back(node2);
          }
        }
        node1->strength += edgeweight[new_edge_id];
        node2->strength += edgeweight[new_edge_id];
        node_vec1[new_edge_id] = node1->id;
        node_vec2[new_edge_id] = node2->id;
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
        updatePreference(q1.front(), params);
        q1.pop();
      }
      qm.clear();
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
    while (! q.empty()) {
      temp_node = q.front();
      q.pop();
      if (temp_node->left != NULL) {
        q.push(temp_node->left);
      }
      if (temp_node->right != NULL) {
        q.push(temp_node->right);
      }
      strength[j] = temp_node->strength;
      pref[j] = temp_node->p;
      // free memory (node and tree)
      delete(temp_node);
      j++;
    }
    // free memory (queue)
    queue<node*>().swap(q);
  }
}