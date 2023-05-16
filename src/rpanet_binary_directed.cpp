#include <iostream>
#include <queue>
#include <R.h>
#include <Rcpp.h>
#include "rpanet_binary_linear.h"

using namespace std;
funcPtrD custmSourcePref;
funcPtrD custmTargetPref;

/**
 * Node structure in directed networks.
 * id: node id
 * outs, ins: out- and in-strength
 * sourcep: preference of being chosen as a source node
 * targetp: preference of being chosen as a target node
 * total_sourcep: sum of sourcep of current node and its children
 * total_targetp: sum of targetp of current node and its children
 * *left, *right, *parent: pointers to its left, right and parent
 */
struct node_d
{
  int id, group;
  double outs, ins;
  double sourcep, targetp, total_sourcep, total_targetp;
  node_d *left, *right, *parent;
};

/**
 * Update total source preference from current node to root.
 *
 * @param current_node The current node.
 */
void updateTotalSourcep(node_d *current_node)
{
  if (current_node->left == NULL)
  {
    current_node->total_sourcep = current_node->sourcep;
  }
  else if (current_node->right == NULL)
  {
    current_node->total_sourcep = current_node->sourcep + current_node->left->total_sourcep;
  }
  else
  {
    current_node->total_sourcep = current_node->sourcep + current_node->left->total_sourcep + current_node->right->total_sourcep;
  }
  while (current_node->parent != NULL)
  {
    return updateTotalSourcep(current_node->parent);
  }
}

/**
 * Update total target preference from current node to root.
 *
 * @param current_node The current node.
 */
void updateTotalTargetp(node_d *current_node)
{
  if (current_node->left == NULL)
  {
    current_node->total_targetp = current_node->targetp;
  }
  else if (current_node->right == NULL)
  {
    current_node->total_targetp = current_node->targetp + current_node->left->total_targetp;
  }
  else
  {
    current_node->total_targetp = current_node->targetp + current_node->left->total_targetp + current_node->right->total_targetp;
  }
  while (current_node->parent != NULL)
  {
    return updateTotalTargetp(current_node->parent);
  }
}

/**
 * Update node preference and total preference from the sampled node to root.
 *
 * @param temp_node The sampled node.
 * @param func_type Default or customized preference function.
 * @param sparams Parameters passed to the default source preference function.
 * @param tparams Parameters passed to the default target preference function.
 * @param custmSourcePref Pointer of customized source preference function.
 * @param custmTargetPref Pointer of customized target preference function.
 */
void updatePrefD(node_d *temp_node, int func_type,
                 double *sparams, double *tparams,
                 funcPtrD custmSourcePref,
                 funcPtrD custmTargetPref)
{
  double temp_sourcep = temp_node->sourcep, temp_targetp = temp_node->targetp;
  if (func_type == 1)
  {
    temp_node->sourcep = prefFuncD(temp_node->outs, temp_node->ins, sparams);
    temp_node->targetp = prefFuncD(temp_node->outs, temp_node->ins, tparams);
  }
  else
  {
    temp_node->sourcep = custmSourcePref(temp_node->outs, temp_node->ins);
    temp_node->targetp = custmTargetPref(temp_node->outs, temp_node->ins);
  }

  if ((temp_node->sourcep < 0) || (temp_node->targetp < 0))
  {
    Rcpp::stop("Negative preference score returned, please check your preference function(s).");
  }

  if (temp_node->sourcep != temp_sourcep)
  {
    updateTotalSourcep(temp_node);
  }
  if (temp_node->targetp != temp_targetp)
  {
    updateTotalTargetp(temp_node);
  }
}

/**
 * Create a new node.
 *
 * @param id Node ID.
 *
 * @return The new node.
 */
node_d *createNodeD(int id)
{
  node_d *new_node = new node_d();
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
node_d *insertNodeD(queue<node_d *> &q, int new_node_id)
{
  node_d *new_node = createNodeD(new_node_id);
  node_d *temp_node = q.front();
  if (temp_node->left == NULL)
  {
    temp_node->left = new_node;
  }
  else if (temp_node->right == NULL)
  {
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
node_d *findSourceNode(node_d *root, double w)
{
  if (w > root->total_sourcep)
  {
    // numerical error
    // Rprintf("Numerical error. Node %d. Diff %f.\n", root->id, (w - root->total_sourcep) * pow(10, 10));
    w = root->total_sourcep;
  }
  w -= root->sourcep;
  if (w <= 0)
  {
    return root;
  }
  else
  {
    if (w > root->left->total_sourcep)
    {
      return findSourceNode(root->right, w - root->left->total_sourcep);
    }
    else
    {
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
node_d *findTargetNode(node_d *root, double w)
{
  if (w > root->total_targetp)
  {
    // numerical error
    // Rprintf("Numerical error. Node %d. Diff %f.\n", root->id, (w - root->total_targetp) * pow(10, 10));
    w = root->total_targetp;
  }
  w -= root->targetp;
  if (w <= 0)
  {
    return root;
  }
  else
  {
    if (w > root->left->total_targetp)
    {
      return findTargetNode(root->right, w - root->left->total_targetp);
    }
    else
    {
      return findTargetNode(root->left, w);
    }
  }
}

/**
 * Sample a source/target node from the tree.
 *
 * @param root Root node of the tree.
 * @param type Source node or target node.
 *
 * @return Sampled source/target node.
 */
node_d *sampleNodeD(node_d *root, char type)
{
  double w;
  w = 1;
  while (w == 1)
  {
    w = unif_rand();
  }
  if (type == 's')
  {
    w *= root->total_sourcep;
    return findSourceNode(root, w);
  }
  else
  {
    w *= root->total_targetp;
    return findTargetNode(root, w);
  }
}

//' Preferential attachment algorithm.
//'
//' @param nstep Number of steps.
//' @param m Number of new edges in each step.
//' @param new_node_id New node ID.
//' @param new_edge_id New edge ID.
//' @param source_node Sequence of source nodes.
//' @param target_node Sequence of target nodes.
//' @param outs Sequence of out-strength.
//' @param ins Sequence of in-strength.
//' @param edgeweight Weight of existing and new edges.
//' @param scenario Scenario of existing and new edges.
//' @param sample_recip Logical, whether reciprocal edges will be added.
//' @param node_group Sequence of node group.
//' @param spref Sequence of node source preference.
//' @param tpref Sequence of node target preference.
//' @param control List of controlling arguments.
//' @return Sampled network.
//'
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List rpanet_binary_directed(
    int nstep,
    Rcpp::IntegerVector m,
    int new_node_id,
    int new_edge_id,
    Rcpp::IntegerVector source_node,
    Rcpp::IntegerVector target_node,
    Rcpp::NumericVector outs,
    Rcpp::NumericVector ins,
    Rcpp::NumericVector edgeweight,
    Rcpp::IntegerVector scenario,
    bool sample_recip,
    Rcpp::IntegerVector node_group,
    Rcpp::NumericVector spref_vec,
    Rcpp::NumericVector tpref_vec,
    Rcpp::List control)
{
  Rcpp::List scenario_ctl = control["scenario"];
  double alpha = scenario_ctl["alpha"];
  double beta = scenario_ctl["beta"];
  double gamma = scenario_ctl["gamma"];
  double xi = scenario_ctl["xi"];
  bool beta_loop = scenario_ctl["beta.loop"];
  bool source_first = scenario_ctl["source.first"];
  Rcpp::List newedge_ctl = control["newedge"];
  // bool node_unique = ! newedge_ctl["node.replace"];
  bool snode_unique = !newedge_ctl["snode.replace"];
  bool tnode_unique = !newedge_ctl["tnode.replace"];
  Rcpp::List reciprocal_ctl = control["reciprocal"];
  bool selfloop_recip = reciprocal_ctl["selfloop.recip"];
  Rcpp::NumericVector group_prob_vec = reciprocal_ctl["group.prob"];
  double *group_prob = &(group_prob_vec[0]);
  Rcpp::NumericMatrix recip_prob = reciprocal_ctl["recip.prob"];
  Rcpp::List preference_ctl = control["preference"];
  Rcpp::NumericVector sparams_vec(5);
  Rcpp::NumericVector tparams_vec(5);
  double *sparams, *tparams;
  double *spref = &(spref_vec[0]);
  double *tpref = &(tpref_vec[0]);
  // different types of preference functions
  int func_type = preference_ctl["ftype.temp"];
  switch (func_type)
  {
  case 1:
    sparams_vec = preference_ctl["sparams"];
    tparams_vec = preference_ctl["tparams"];
    sparams = &(sparams_vec[0]);
    tparams = &(tparams_vec[0]);
    break;
  case 2:
  {
    SEXP source_pref_func_ptr = preference_ctl["spref.pointer"];
    custmSourcePref = *Rcpp::XPtr<funcPtrD>(source_pref_func_ptr);
    SEXP target_pref_func_ptr = preference_ctl["tpref.pointer"];
    custmTargetPref = *Rcpp::XPtr<funcPtrD>(target_pref_func_ptr);
    break;
  }
  }

  double u, p, temp_p;
  bool m_error;
  int i, j, n_existing, current_scenario, n_reciprocal;
  node_d *node1, *node2;

  // re-order label nodes according to source preference and target preference
  Rcpp::NumericVector temp_source_pref(new_node_id);
  Rcpp::NumericVector temp_target_pref(new_node_id);
  Rcpp::IntegerVector sorted_node = Rcpp::seq(0, new_node_id - 1);
  if (func_type == 1)
  {
    for (i = 0; i < new_node_id; i++)
    {
      temp_source_pref[i] = prefFuncD(outs[i], ins[i], sparams);
      temp_target_pref[i] = prefFuncD(outs[i], ins[i], tparams);
    }
  }
  else
  {
    for (i = 0; i < new_node_id; i++)
    {
      temp_source_pref[i] = custmSourcePref(outs[i], ins[i]);
      temp_target_pref[i] = custmTargetPref(outs[i], ins[i]);
    }
  }
  if (alpha < gamma)
  {
    sort(sorted_node.begin(), sorted_node.end(),
         [&](int k, int l){ return temp_source_pref[k] > temp_source_pref[l]; });
  }
  else
  {
    sort(sorted_node.begin(), sorted_node.end(),
         [&](int k, int l){ return temp_target_pref[k] > temp_target_pref[l]; });
  }

  // initialize a tree from the seed graph
  j = sorted_node[0];
  node_d *root = createNodeD(j);
  root->outs = outs[j];
  root->ins = ins[j];
  root->group = node_group[j];
  updatePrefD(root, func_type, sparams, tparams, custmSourcePref, custmTargetPref);
  queue<node_d *> q, q1;
  q.push(root);
  for (int i = 1; i < new_node_id; i++)
  {
    j = sorted_node[i];
    node1 = insertNodeD(q, j);
    node1->outs = outs[j];
    node1->ins = ins[j];
    node1->group = node_group[j];
    updatePrefD(node1, func_type, sparams, tparams, custmSourcePref, custmTargetPref);
  }
  // sample edges
  // GetRNGstate();
  for (i = 0; i < nstep; i++)
  {
    n_reciprocal = 0;
    m_error = false;
    n_existing = new_node_id;
    for (j = 0; j < m[i]; j++)
    {
      u = unif_rand();
      if (u <= alpha)
      {
        current_scenario = 1;
      }
      else if (u <= alpha + beta)
      {
        current_scenario = 2;
      }
      else if (u <= alpha + beta + gamma)
      {
        current_scenario = 3;
      }
      else if (u <= alpha + beta + gamma + xi)
      {
        current_scenario = 4;
      }
      else
      {
        current_scenario = 5;
      }
      switch (current_scenario)
      {
      case 1:
        if (root->total_targetp == 0)
        {
          m_error = true;
          break;
        }
        node1 = insertNodeD(q, new_node_id);
        if (sample_recip)
        {
          node1->group = sampleGroup(group_prob);
        }
        new_node_id++;
        node2 = sampleNodeD(root, 't');
        break;
      case 2:
        if ((root->total_targetp == 0) || (root->total_sourcep == 0))
        {
          m_error = true;
          break;
        }
        if (source_first)
        {
          node1 = sampleNodeD(root, 's');
          if (beta_loop)
          {
            node2 = sampleNodeD(root, 't');
          }
          else
          {
            if (node1->targetp == root->total_targetp)
            {
              m_error = true;
              break;
            }
            if (node1->targetp == 0)
            {
              node2 = sampleNodeD(root, 't');
            }
            else
            {
              temp_p = node1->targetp;
              node1->targetp = 0;
              updateTotalTargetp(node1);
              node2 = sampleNodeD(root, 't');
              node1->targetp = temp_p;
              updateTotalTargetp(node1);
            }
          }
        }
        else
        {
          node2 = sampleNodeD(root, 't');
          if (beta_loop)
          {
            node1 = sampleNodeD(root, 's');
          }
          else
          {
            if (node2->sourcep == root->total_sourcep)
            {
              m_error = true;
              break;
            }
            if (node2->sourcep == 0)
            {
              node1 = sampleNodeD(root, 's');
            }
            else
            {
              temp_p = node2->sourcep;
              node2->sourcep = 0;
              updateTotalSourcep(node2);
              node1 = sampleNodeD(root, 's');
              node2->sourcep = temp_p;
              updateTotalSourcep(node2);
            }
          }
        }
        break;
      case 3:
        if (root->total_sourcep == 0)
        {
          m_error = true;
          break;
        }
        node1 = sampleNodeD(root, 's');
        node2 = insertNodeD(q, new_node_id);
        if (sample_recip)
        {
          node2->group = sampleGroup(group_prob);
        }
        new_node_id++;
        break;
      case 4:
        node1 = insertNodeD(q, new_node_id);
        new_node_id++;
        node2 = insertNodeD(q, new_node_id);
        new_node_id++;
        if (sample_recip)
        {
          node1->group = sampleGroup(group_prob);
          node2->group = sampleGroup(group_prob);
        }
        break;
      case 5:
        node1 = node2 = insertNodeD(q, new_node_id);
        if (sample_recip)
        {
          node1->group = sampleGroup(group_prob);
        }
        new_node_id++;
        break;
      }
      if (m_error)
      {
        break;
      }
      // if (node_unique) {
      //   if (node1->id < n_existing) {
      //     node1->sourcep = 0;
      //     node1->targetp = 0;
      //     updateTotalSourcep(node1);
      //     updateTotalTargetp(node1);
      //   }
      //   if ((node2->id < n_existing) && (node1 != node2)) {
      //     node2->sourcep = 0;
      //     node2->targetp = 0;
      //     updateTotalSourcep(node2);
      //     updateTotalTargetp(node2);
      //   }
      // }
      // else {
      //   if (snode_unique && (node1->id < n_existing)) {
      //     node1->sourcep = 0;
      //     updateTotalSourcep(node1);
      //   }
      //   if (tnode_unique && (node2->id < n_existing)) {
      //     node2->targetp = 0;
      //     updateTotalTargetp(node2);
      //   }
      // }
      // sample without replacement
      if (snode_unique && (node1->id < n_existing))
      {
        node1->sourcep = 0;
        updateTotalSourcep(node1);
      }
      if (tnode_unique && (node2->id < n_existing))
      {
        node2->targetp = 0;
        updateTotalTargetp(node2);
      }
      node1->outs += edgeweight[new_edge_id];
      node2->ins += edgeweight[new_edge_id];
      source_node[new_edge_id] = node1->id;
      target_node[new_edge_id] = node2->id;
      scenario[new_edge_id] = current_scenario;
      q1.push(node1);
      q1.push(node2);
      // handle reciprocal
      if (sample_recip)
      {
        if ((node1->id != node2->id) || selfloop_recip)
        {
          p = unif_rand();
          if (p <= recip_prob(node2->group, node1->group))
          {
            new_edge_id++;
            n_reciprocal++;
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
    m[i] += n_reciprocal;
    if (m_error)
    {
      m[i] = j + n_reciprocal;
      Rprintf("No enough unique nodes for a scenario %d edge at step %d. Added %d edge(s) at current step.\n",
              current_scenario, i + 1, m[i]);
    }
    while (!q1.empty())
    {
      updatePrefD(q1.front(), func_type, sparams, tparams, custmSourcePref, custmTargetPref);
      q1.pop();
    }
  }
  // PutRNGstate();
  // free memory (queue)
  queue<node_d *>().swap(q);
  queue<node_d *>().swap(q1);
  // save strength and preference
  q.push(root);
  while (!q.empty())
  {
    node1 = q.front();
    q.pop();
    if (node1->right != NULL)
    {
      q.push(node1->left);
      q.push(node1->right);
    }
    else if (node1->left != NULL)
    {
      q.push(node1->left);
    }
    j = node1->id;
    outs[j] = node1->outs;
    ins[j] = node1->ins;
    node_group[j] = node1->group;
    spref[j] = node1->sourcep;
    tpref[j] = node1->targetp;
    // free memory (node and tree)
    delete node1;
  }
  // free memory (queue)
  queue<node_d *>().swap(q);

  Rcpp::List ret;
  ret["m"] = m;
  ret["nnode"] = new_node_id;
  ret["nedge"] = new_edge_id;
  ret["node_vec1"] = source_node;
  ret["node_vec2"] = target_node;
  ret["outs"] = outs;
  ret["ins"] = ins;
  ret["scenario"] = scenario;
  ret["nodegroup"] = node_group;
  ret["spref"] = spref_vec;
  ret["tpref"] = tpref_vec;
  return ret;
}
