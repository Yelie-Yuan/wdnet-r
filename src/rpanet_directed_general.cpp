#include<iostream>
#include<queue>
#include<math.h>
#include<R.h>
using namespace std;

// 1. user defined preference functions (sourcePreferenceFunc and 
//    targetPreferenceFunc); how to pass R functions in c++?
// 2. w can not equal to 1 in function sampleNode, otherwise findNode returns error 
//    because of numeric precision
// 3. add a parameter p to control reciprocal edges
// 4. add a parameter m to control number of new edges per step

// node structure
// id: node id
// outs, ins: out- and in-strength
// sourcep: preference of being chosen as a source node
// targetp: preference of being chosed as a target node
// totalSourcep: sum of sourcep of current node and its children
// totalTargetp: sum of targetp of current node and its children
// *left, *right, *parent: pointers to its left, right and parent
struct node {
  int id;
  double outs, ins;
  double sourcep, targetp, totalSourcep, totalTargetp;
  node *left, *right, *parent;
};

// preference functions
double sourcePreferenceFunc(double outs, double ins, double *source_params) {
  return source_params[0] * pow(outs, source_params[1]) + 
    source_params[2] * pow(ins, source_params[3]) + source_params[4];
}
double targetPreferenceFunc(double outs, double ins, double *target_params) {
  return target_params[0] * pow(ins, target_params[1]) + 
    target_params[2] * pow(outs, target_params[3]) + target_params[4];
}

// create a new node
node *createNode(int id, double outs, double ins, 
      double *source_params, double *target_params) {
  node *newNode = new node();
  newNode->id = id;
  newNode->outs = outs;
  newNode->ins = ins;
  newNode->sourcep = newNode->totalSourcep = sourcePreferenceFunc(outs, ins, source_params);
  newNode->targetp = newNode->totalTargetp = targetPreferenceFunc(outs, ins, target_params);
  newNode->left = newNode->right = newNode->parent = NULL;
  return newNode;
}

// update total source preference from current node to root
void updateTotalSourcep(node *currentNode, double increment) {
  currentNode->totalSourcep += increment;
  while(currentNode->id > 0) {
    return updateTotalSourcep(currentNode->parent, increment);
  }
}
// update total source preference from current node to root
void updateTotalTargetp(node *currentNode, double increment) {
  currentNode->totalTargetp += increment;
  while(currentNode->id > 0) {
    return updateTotalTargetp(currentNode->parent, increment);
  }
}

// insert a new node to the tree
void insert(queue<node*> &q, int new_node_id, double outs, double ins, 
      double *source_params, double *target_params) {
  node *newNode = createNode(new_node_id, outs, ins, source_params, target_params);
  node *tempNode = q.front();
  if(tempNode->left == NULL) {
    tempNode->left = newNode;
  }
  else if (tempNode->right == NULL) {
    tempNode->right = newNode;
    q.pop();
  }
  newNode->parent = tempNode;
  q.push(newNode);
  if (newNode->totalSourcep != 0) {
    updateTotalSourcep(newNode->parent, newNode->totalSourcep);
  }
  if (newNode->totalTargetp != 0) {
    updateTotalTargetp(newNode->parent, newNode->totalTargetp);
  }
}

// find source node with a given critical point w
node *findSourceNode(node *root, double w) {
  w -= root->sourcep;
  if (w <= 0) {
    return root;
  } 
  else {
    if (w > root->left->totalSourcep) {
      return findSourceNode(root->right, w - root->left->totalSourcep);
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
    if (w > root->left->totalTargetp) {
      return findTargetNode(root->right, w - root->left->totalTargetp);
    }
    else {
      return findTargetNode(root->left, w);
    }
  }
}

// sample a node from the tree
node* sampleNode(node *root, char type) {
  GetRNGstate();
  double w = 1;
  node *tempNode;
  while (w == 1) {
    w = unif_rand();
  }
  if (type == 's') {
    w *= root->totalSourcep;
    tempNode = findSourceNode(root, w);
  }
  else {
    w *= root->totalTargetp;
    tempNode = findTargetNode(root, w);
  }
  PutRNGstate();
  return tempNode;
}

// update strength, preference and total preference from the sampled node to root
void updateTree(node *tempNode, double weight, char type, 
      double *source_params, double *target_params) {
  double increment;
  if (type == 's') {
    tempNode->outs += weight;
  }
  else {
    tempNode->ins += weight;
  }
  increment = sourcePreferenceFunc(tempNode->outs, tempNode->ins,
    source_params) - tempNode->sourcep;
  if (increment != 0) {
    tempNode->sourcep += increment;
    updateTotalSourcep(tempNode, increment);
  }
  increment = targetPreferenceFunc(tempNode->outs, tempNode->ins,
    target_params) - tempNode->targetp;
  if (increment != 0) {
    tempNode->targetp += increment;
    updateTotalTargetp(tempNode, increment);
  }
}

// sample edges with a given seed graph
extern "C" {
  void rpanet_directed_general_cpp(int *nstep_ptr, int *new_node_id_ptr, int *new_edge_id_ptr, 
        int *source_node, int *target_node, double *outs, double *ins, double *weight, int *scenario,
        double *alpha_ptr, double *beta_ptr, double *gamma_ptr, double *xi_ptr, 
        double *source_params, double *target_params) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    node *temp_node1, *temp_node2;
    // initialize a tree from the seed graph
    node *root = createNode(0, outs[0], ins[0], source_params, target_params);
    queue<node*> q;
    q.push(root);
    for (int i = 1; i < new_node_id; i++) {
      insert(q, i, outs[i], ins[i], source_params, target_params);
    }
    // sample edges
    GetRNGstate();
    for (int i = 0; i < nstep; i++) {
      u = unif_rand();
      if (u <= alpha) {
        // sample an existing node before inserting a new node
        temp_node2 = sampleNode(root, 't');
        updateTree(temp_node2, weight[new_edge_id], 't', source_params, target_params);
        target_node[new_edge_id] = temp_node2->id;
        insert(q, new_node_id, weight[new_edge_id], 0, source_params, target_params);
        source_node[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 1;
      }
      else if (u <= alpha + beta) {
        temp_node1 = sampleNode(root, 's');
        temp_node2 = sampleNode(root, 't');
        updateTree(temp_node1, weight[new_edge_id], 's', source_params, target_params);
        updateTree(temp_node2, weight[new_edge_id], 't', source_params, target_params);
        source_node[new_edge_id] = temp_node1->id;
        target_node[new_edge_id] = temp_node2->id;
        scenario[new_edge_id] = 2;
      }
      else if (u <= alpha + beta + gamma) {
        temp_node1 = sampleNode(root, 's');
        updateTree(temp_node1, weight[new_edge_id], 's', source_params, target_params);
        source_node[new_edge_id] = temp_node1->id;
        insert(q, new_node_id, 0, weight[new_edge_id], source_params, target_params);
        target_node[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 3;
      }
      else if (u <= alpha + beta + gamma + xi) {
        insert(q, new_node_id, weight[new_edge_id], 0, source_params, target_params);
        source_node[new_edge_id] = new_node_id;
        new_node_id++;
        insert(q, new_node_id, 0, weight[new_edge_id], source_params, target_params);
        target_node[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 4;
      }
      else {
        insert(q, new_node_id, weight[new_edge_id], weight[new_edge_id], source_params, target_params);
        source_node[new_edge_id] = target_node[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 5;
      }
      new_edge_id++;
    }
    PutRNGstate();
    *new_node_id_ptr = new_node_id;
    *new_edge_id_ptr = new_edge_id;
  }
}