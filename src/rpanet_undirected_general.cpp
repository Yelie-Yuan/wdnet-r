#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
#include<bits/stdc++.h>
using namespace std;

// 1. user defined preference functions (sourcePreferenceFunc and 
//    targetPreferenceFunc); how to pass R functions in c++?
// 2. w can not equal to 1 in function SampleNode, otherwise FindNode returns error 
//    because of numeric precision
// 3. add a parameter m to control number of new edges per step

// node structure
struct node {
  int id;
  double strength;
  double p, totalp;
  node *left, *right, *parent;
};

// preference function; import from R
double PreferenceFunc(double strength, double *params) {
  return pow(strength, params[0]) + params[1];
}

// update total preference from current node to root
void AddIncrement(node *current_node, double increment) {
  current_node->totalp += increment;
  while(current_node->id > 0) {
    return AddIncrement(current_node->parent, increment);
  }
}

// update strength and preference of from the sampled node to root
void UpdatePreference(node *temp_node, double *params) {
  double temp_p = temp_node->p;
  temp_node->p = PreferenceFunc(temp_node->strength, params);
  AddIncrement(temp_node, temp_node->p - temp_p);
}

// create a new node
node *CreateNode(int id) {
  node *new_node = new node();
  new_node->id = id;
  new_node->strength = 0;
  new_node->p = new_node->totalp = 0;
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
}

// may need a better way to construct the complete tree
node *InsertNode(queue<node*> &q, int id) {
  node *new_node = CreateNode(id);
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

// find node with a given critical point w
node *FindNode(node *root, double w) {
  w -= root->p;
  if (w <= 0) {
    return root;
  }
  else {
    if (w > root->left->totalp) {
      return FindNode(root->right, w - root->left->totalp);
    }
    else {
      return FindNode(root->left, w);
    }
  }
}

// sample a node from the tree
node *SampleNode(node *root, deque<node*> &qm) {
  double w;
  node *temp_node;
  while (true) {
    w = 1;
    while (w == 1) {
      w = unif_rand();
    }
    w *= root->totalp;
    temp_node = FindNode(root, w);
    if (find(qm.begin(), qm.end(), temp_node) == qm.end()) {
      // if temp_node not in qm
      return temp_node;
    }
  }
}

// sample edges with a given seed graph
// alpha and gamma scenarios are kept for interfacing with R
// alpha: (new, existing); gamma: (existing, new)
// if (directed) {
//   rpanet_directed_general_cpp(alpha, beta, gamma, ...)
// }
// else {
//   rpanet_undirected_general_cpp(alpha, beta, gamma, ...)
// }
extern "C" {
  void rpanet_undirected_general_cpp(
      int *nstep_ptr, int *m,
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *node_vec1, int *node_vec2, 
      double *strength, double *edgeweight, 
      int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr, int *m_unique_ptr,
      double *params, double *pref) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, 
      new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, 
      gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr, m_unique = *m_unique_ptr,
      m_error;
    int i, j, k, n_existing, current_scenario;
    node *node1, *node2;
    // initialize a tree from seed graph
    node *root = CreateNode(0);
    root->strength = strength[0];
    UpdatePreference(root, params);
    queue<node*> q, q1;
    deque<node*> qm;
    q.push(root);
    for (i = 1; i < new_node_id; i++) {
      node1 = InsertNode(q, i);
      node1->strength = strength[i];
      UpdatePreference(node1, params);
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
        if (m_unique) {
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
            node1 = InsertNode(q, new_node_id);
            new_node_id++;
            node2 = SampleNode(root, qm);
            break;
          case 2:
            node1 = SampleNode(root, qm);
            if (! beta_loop) {
              qm.push_back(node1);
              node2 = SampleNode(root, qm);
              qm.pop_back();
            }
            else {
              node2 = SampleNode(root, qm);
            }
            break;
          case 3:
            node1 = SampleNode(root, qm);
            node2 = InsertNode(q, new_node_id);
            new_node_id++;
            break;
          case 4:
            node1 = InsertNode(q, new_node_id);
            new_node_id++;
            node2 = InsertNode(q, new_node_id);
            new_node_id++;
            break;
          case 5:
            node1 = node2 = InsertNode(q, new_node_id);
            new_node_id++;
            break;
        }
        // handle duplicate nodes
        if (m_unique) {
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
        UpdatePreference(q1.front(), params);
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