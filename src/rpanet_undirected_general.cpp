#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
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
node* SampleNode(node *root){
  double w = 1;
  while (w == 1) {
    w = unif_rand();
  }
  w *= root->totalp;
  node *temp_node = FindNode(root, w);
  return temp_node;
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
  void rpanet_undirected_general_cpp(int *nstep_ptr, int *m,
      int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *node_vec1, int *node_vec2, 
      double *strength, double *edgeweight, int *scenario,
      double *alpha_ptr, double *beta_ptr, 
      double *gamma_ptr, double *xi_ptr, 
      int *beta_loop_ptr,
      double *params, double *pref) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    bool beta_loop = *beta_loop_ptr;
    int i, j;
    node *node1, *node2;
    // initialize a tree from seed graph
    node *root = CreateNode(0);
    root->strength = strength[0];
    UpdatePreference(root, params);
    queue<node*> q;
    queue<node*> qm;
    q.push(root);
    for (i = 1; i < new_node_id; i++) {
      node1 = InsertNode(q, i);
      node1->strength = strength[i];
      UpdatePreference(node1, params);
    }
    // sample edges
    GetRNGstate();
    for (i = 0; i < nstep; i++) {
      for (j = 0; j < m[i]; j++) {
        u = unif_rand();
        if (u <= alpha) {
          node1 = InsertNode(q, new_node_id);
          new_node_id++;
          node2 = SampleNode(root);
          scenario[new_edge_id] = 1;
        }
        else if (u <= alpha + beta) {
          node1 = SampleNode(root);
          node2 = SampleNode(root);
          if (! beta_loop) {
            while (node1 == node2) {
              node2 = SampleNode(root);
            }
          }
          scenario[new_edge_id] = 2;
        }
        else if (u <= alpha + beta + gamma) {
          node1 = SampleNode(root);
          node2 = InsertNode(q, new_node_id);
          new_node_id++;
          scenario[new_edge_id] = 3;
        }
        else if (u <= alpha + beta + gamma + xi) {
          node1 = InsertNode(q, new_node_id);
          new_node_id++;
          node2 = InsertNode(q, new_node_id);
          new_node_id++;
          scenario[new_edge_id] = 4;
        }
        else {
          node1 = node2 = InsertNode(q, new_node_id);
          new_node_id++;
          scenario[new_edge_id] = 5;
        }
        // handle duplicate nodes (sampled nodes from the same step should be unique)
        // to be added
        node1->strength += edgeweight[new_edge_id];
        node2->strength += edgeweight[new_edge_id];
        node_vec1[new_edge_id] = node1->id;
        node_vec2[new_edge_id] = node2->id;
        qm.push(node1);
        qm.push(node2);
        new_edge_id++;
      }
      while (! qm.empty()) {
        UpdatePreference(qm.front(), params);
        qm.pop();
      }
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
    // save strength and preference
    queue<node*> q2;
    q2.push(root);
    node *temp_node;
    j = 0;
    while (! q2.empty()) {
      temp_node = q2.front();
      q2.pop();
      if (temp_node->left != NULL) {
        q2.push(temp_node->left);
      }
      if (temp_node->right != NULL) {
        q2.push(temp_node->right);
      }
      strength[j] = temp_node->strength;
      pref[j] = temp_node->p;
      j++;
    }
  }
}