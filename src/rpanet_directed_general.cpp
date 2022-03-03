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
  node *new_node = new node();
  new_node->id = id;
  new_node->outs = outs;
  new_node->ins = ins;
  new_node->sourcep = new_node->totalSourcep = sourcePreferenceFunc(outs, ins, source_params);
  new_node->targetp = new_node->totalTargetp = targetPreferenceFunc(outs, ins, target_params);
  new_node->left = new_node->right = new_node->parent = NULL;
  return new_node;
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
  node *new_node = createNode(new_node_id, outs, ins, source_params, target_params);
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
  if (new_node->totalSourcep != 0) {
    updateTotalSourcep(new_node->parent, new_node->totalSourcep);
  }
  if (new_node->totalTargetp != 0) {
    updateTotalTargetp(new_node->parent, new_node->totalTargetp);
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
  double w = 1;
  node *temp_node;
  while (w == 1) {
    w = unif_rand();
  }
  if (type == 's') {
    w *= root->totalSourcep;
    temp_node = findSourceNode(root, w);
  }
  else {
    w *= root->totalTargetp;
    temp_node = findTargetNode(root, w);
  }
  return temp_node;
}

// update strength, preference and total preference from the sampled node to root
void updateTree(node *temp_node, double out_increase, double in_increase, 
      double *source_params, double *target_params) {
  double increment;
  temp_node->outs += out_increase;
  temp_node->ins += in_increase;
  increment = sourcePreferenceFunc(temp_node->outs, temp_node->ins,
    source_params) - temp_node->sourcep;
  if (increment != 0) {
    temp_node->sourcep += increment;
    updateTotalSourcep(temp_node, increment);
  }
  increment = targetPreferenceFunc(temp_node->outs, temp_node->ins,
    target_params) - temp_node->targetp;
  if (increment != 0) {
    temp_node->targetp += increment;
    updateTotalTargetp(temp_node, increment);
  }
}

// sample edges with a given seed graph
extern "C" {
  void rpanet_directed_general_cpp(int *nstep_ptr, int *new_node_id_ptr, int *new_edge_id_ptr, 
        int *source_node, int *target_node, double *outs, double *ins, double *weight, int *scenario,
        double *alpha_ptr, double *beta_ptr, double *gamma_ptr, double *xi_ptr, 
        double *source_params, double *target_params, 
        double *source_pref, double *target_pref) {
    double u;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr, new_edge_id = *new_edge_id_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    node *temp_node1, *temp_node2;
    // initialize a tree from the seed graph
    node *root = createNode(0, outs[0], ins[0], source_params, target_params);
    queue<node*> q;
    queue<node*> q2;
    q.push(root);
    q2.push(root);
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
        updateTree(temp_node2, 0, weight[new_edge_id], source_params, target_params);
        target_node[new_edge_id] = temp_node2->id;
        insert(q, new_node_id, weight[new_edge_id], 0, source_params, target_params);
        source_node[new_edge_id] = new_node_id;
        new_node_id++;
        scenario[new_edge_id] = 1;
      }
      else if (u <= alpha + beta) {
        temp_node1 = sampleNode(root, 's');
        temp_node2 = sampleNode(root, 't');
        updateTree(temp_node1, weight[new_edge_id], 0, source_params, target_params);
        updateTree(temp_node2, 0, weight[new_edge_id], source_params, target_params);
        source_node[new_edge_id] = temp_node1->id;
        target_node[new_edge_id] = temp_node2->id;
        scenario[new_edge_id] = 2;
      }
      else if (u <= alpha + beta + gamma) {
        temp_node1 = sampleNode(root, 's');
        updateTree(temp_node1, weight[new_edge_id], 0, source_params, target_params);
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
      outs[j] = temp_node->outs;
      ins[j] = temp_node->ins;
      source_pref[j] = temp_node->sourcep;
      target_pref[j] = temp_node->targetp;
      j++;
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

// add reciprocal edges

extern "C" {
  void rpanet_directed_general_recip_cpp(int *nstep_ptr, int *new_node_id_ptr, int *new_edge_id_ptr, 
        int *source_node, int *target_node, double *outs, double *ins, double *weight, int *scenario,
        double *alpha_ptr, double *beta_ptr, double *gamma_ptr, double *xi_ptr, 
        double *source_params, double *target_params, double *group_dist, double *recip, int *group, 
        int *ngroup_ptr, double *source_pref, double *target_pref) {
    double u, p;
    int nstep = *nstep_ptr, new_node_id = *new_node_id_ptr,
        new_edge_id = *new_edge_id_ptr, ngroup = *ngroup_ptr;
    double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr, xi = *xi_ptr;
    int temp_group1, temp_group2;
    node *temp_node1, *temp_node2;
    // initialize a tree from the seed graph
    node *root = createNode(0, outs[0], ins[0], source_params, target_params);
    queue<node*> q;
    queue<node*> q2;
    q.push(root);
    q2.push(root);
    for (int i = 1; i < new_node_id; i++) {
      insert(q, i, outs[i], ins[i], source_params, target_params);
    }
    // sample edges
    GetRNGstate();
    for (int i = 0; i < nstep; i++) {
      u = unif_rand();
      p = unif_rand();
      if (u <= alpha) {
        // sample an existing node before inserting a new node
        temp_node2 = sampleNode(root, 't');
        temp_group1 = sampleGroup(group_dist);
        group[new_node_id] = temp_group1;
        target_node[new_edge_id] = temp_node2->id;
        source_node[new_edge_id] = new_node_id;
        scenario[new_edge_id] = 1;
        if (p <= recip[group[temp_node2->id] * ngroup + temp_group1]) {
          updateTree(temp_node2, weight[new_edge_id + 1], weight[new_edge_id], source_params, target_params);
          insert(q, new_node_id, weight[new_edge_id], weight[new_edge_id + 1], source_params, target_params);
          new_edge_id++;
          scenario[new_edge_id] = 6;
          target_node[new_edge_id] = new_node_id;
          source_node[new_edge_id] = temp_node2->id;
        }
        else {
          updateTree(temp_node2, 0, weight[new_edge_id], source_params, target_params);
          insert(q, new_node_id, weight[new_edge_id], 0, source_params, target_params);
        }
        new_node_id++;
      }
      else if (u <= alpha + beta) {
        temp_node1 = sampleNode(root, 's');
        temp_node2 = sampleNode(root, 't');
        source_node[new_edge_id] = temp_node1->id;
        target_node[new_edge_id] = temp_node2->id;
        scenario[new_edge_id] = 2;
        if (p <= recip[group[temp_node2->id] * ngroup + group[temp_node1->id]]) {
          updateTree(temp_node1, weight[new_edge_id], weight[new_edge_id + 1], source_params, target_params);
          updateTree(temp_node2, weight[new_edge_id + 1], weight[new_edge_id], source_params, target_params);
          new_edge_id++;
          scenario[new_edge_id] = 6;
          source_node[new_edge_id] = temp_node2->id;
          target_node[new_edge_id] = temp_node1->id;
        }
        else {
          updateTree(temp_node1, weight[new_edge_id], 0, source_params, target_params);
          updateTree(temp_node2, 0, weight[new_edge_id], source_params, target_params);
        }
      }
      else if (u <= alpha + beta + gamma) {
        temp_node1 = sampleNode(root, 's');
        temp_group2 = sampleGroup(group_dist);
        group[new_node_id] = temp_group2;
        source_node[new_edge_id] = temp_node1->id;
        target_node[new_edge_id] = new_node_id;
        scenario[new_edge_id] = 3;
        if (p <= recip[temp_group2 * ngroup + group[temp_node1->id]]) {
          updateTree(temp_node1, weight[new_edge_id], weight[new_edge_id + 1], source_params, target_params);
          insert(q, new_node_id, weight[new_edge_id + 1], weight[new_edge_id], source_params, target_params);
          new_edge_id++;
          source_node[new_edge_id] = new_node_id;
          target_node[new_edge_id] = temp_node1->id;
          scenario[new_edge_id] = 6;
        }
        else {
          updateTree(temp_node1, weight[new_edge_id], 0, source_params, target_params);
          insert(q, new_node_id, 0, weight[new_edge_id], source_params, target_params);
        }
        new_node_id++;
      }
      else if (u <= alpha + beta + gamma + xi) {
        temp_group1 = sampleGroup(group_dist);
        temp_group2 = sampleGroup(group_dist);
        group[new_node_id] = temp_group1;
        source_node[new_edge_id] = new_node_id;
        new_node_id++;
        group[new_node_id] = temp_group2;
        target_node[new_edge_id] = new_node_id;
        scenario[new_edge_id] = 4;
        if (p <= recip[temp_group2 * ngroup + temp_group1]) {
          insert(q, new_node_id - 1, weight[new_edge_id], weight[new_edge_id + 1], source_params, target_params);
          insert(q, new_node_id, weight[new_edge_id + 1], weight[new_edge_id], source_params, target_params);
          new_edge_id++;
          source_node[new_edge_id] = new_node_id;
          target_node[new_edge_id] = new_node_id - 1;
          scenario[new_edge_id] = 6;
        }
        else {
          insert(q, new_node_id - 1, weight[new_edge_id], 0, source_params, target_params);
          insert(q, new_node_id, 0, weight[new_edge_id], source_params, target_params);
        }
        new_node_id++;
      }
      else {
        temp_group1 = sampleGroup(group_dist);
        group[new_node_id] = temp_group1;
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
      outs[j] = temp_node->outs;
      ins[j] = temp_node->ins;
      source_pref[j] = temp_node->sourcep;
      target_pref[j] = temp_node->targetp;
      j++;
    }
  }
}