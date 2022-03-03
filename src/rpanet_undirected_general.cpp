#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
using namespace std;

// 1. user defined preference functions (sourcePreferenceFunc and 
//    targetPreferenceFunc); how to pass R functions in c++?
// 2. w can not equal to 1 in function sampleNode, otherwise findNode returns error 
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
double preferenceFunc(double strength, double *params) {
  return pow(strength, params[0]) + params[1];
}

// create a new node
node *createNode(int id, double strength, double *params) {
  node *new_node = new node();
  new_node->id = id;
  new_node->strength = strength;
  new_node->p = new_node->totalp = preferenceFunc(strength, params);
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
}

// update total preference from current node to root
void updateTotalp(node *currentNode, double increment) {
  currentNode->totalp += increment;
  while(currentNode->id > 0) {
    return updateTotalp(currentNode->parent, increment);
  }
}

// may need a better way to construct the complete tree
void insert(queue<node*> &q, int id, double strength, double *params) {
  node *new_node = createNode(id, strength, params);
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
  updateTotalp(new_node->parent, new_node->totalp);
}

// find node with a given critical point w
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

// sample a node from the tree
node* sampleNode(node *root){
  double w = 1;
  while (w == 1) {
    w = unif_rand();
  }
  w *= root->totalp;
  node *temp_node = findNode(root, w);
  return temp_node;
}

// update strength and preference of from the sampled node to root
void updateTree(node *temp_node, double edgeweight, double *params) {
  temp_node->strength += edgeweight;
  double tp = temp_node->p;
  temp_node->p = preferenceFunc(temp_node->strength, params);
  updateTotalp(temp_node, temp_node->p - tp);
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
  void rpanet_undirected_general_cpp(int *nstep_ptr, int *new_node_id_ptr, int *new_edge_id_ptr, 
      int *node_vec1, int *node_vec2, double *strength, double *edgeweight, int *scenario,
      double *alpha_ptr, double *beta_ptr, double *gamma_ptr, double *xi_ptr, double *params, 
      double *pref) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    node *temp_node1, *temp_node2;
    // initialize a tree from seed graph
    node *root = createNode(0, strength[0], params);
    queue<node*> q;
    queue<node*> q2;
    q.push(root);
    q2.push(root);
    for (int i = 1; i < new_node_id; i++) {
      insert(q, i, strength[i], params);
    }
    // sample edges
    GetRNGstate();
    for (int i = 0; i < nstep; i++) {
      u = unif_rand();
      if (u <= alpha) {
        // sample an existing node before inserting a new node
        temp_node2 = sampleNode(root);
        updateTree(temp_node2, edgeweight[new_edge_id], params);
        node_vec2[new_edge_id] = temp_node2->id;
        insert(q, new_node_id, edgeweight[new_edge_id], params);
        node_vec1[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 1;
      }
      else if (u <= alpha + beta) {
        temp_node1 = sampleNode(root);
        temp_node2 = sampleNode(root);
        updateTree(temp_node1, edgeweight[new_edge_id], params);
        updateTree(temp_node2, edgeweight[new_edge_id], params);
        node_vec1[new_edge_id] = temp_node1->id;
        node_vec2[new_edge_id] = temp_node2->id;
        scenario[new_edge_id] = 2;
      }
      else if (u <= alpha + beta + gamma) {
        temp_node1 = sampleNode(root);
        updateTree(temp_node1, edgeweight[new_edge_id], params);
        node_vec1[new_edge_id] = temp_node1->id;
        insert(q, new_node_id, edgeweight[new_edge_id], params);
        node_vec2[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 3;
      }
      else if (u <= alpha + beta + gamma + xi) {
        insert(q, new_node_id, edgeweight[new_edge_id], params);
        node_vec1[new_edge_id] = new_node_id;
        new_node_id++;
        insert(q, new_node_id, edgeweight[new_edge_id], params);
        node_vec2[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 4;
      }
      else {
        insert(q, new_node_id, edgeweight[new_edge_id] * 2, params);
        node_vec1[new_edge_id] = new_node_id;
        node_vec2[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 5;
      }
      new_edge_id++;
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
    // save strength and preference
    node *temp_node;
    int j = 0;
    while (! q2.empty())
    {
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