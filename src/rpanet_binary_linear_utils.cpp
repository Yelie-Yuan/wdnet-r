# include <math.h>
# include <R.h>


/**
 * Default source/target preference function.
 *
 * @param outs Node out-strength.
 * @param ins Node in-strength.
 * @param params Parameters passed to the source/target preference function.
 *
 * @return Source preference of a node.
 */
double prefFuncD(double outs, double ins, double *params)
{
  return params[0] * pow(outs, params[1]) +
         params[2] * pow(ins, params[3]) + params[4];
}

/**
 * Default preference function.
 *
 * @param s Node strength.
 * @param params Parameters passed to the preference function.
 *
 * @return Preference of a node.
 */
double prefFuncUnd(double s, double *params)
{
  return pow(s, params[0]) + params[1];
}

/**
 * Sample a node group. Directed networks only.
 *
 * @param group_prob Probability weights for sampling the group of new nodes.
 *
 * @return Sampled group for the new node.
 */
int sampleGroup(double *group_prob)
{
  double g = 0;
  int i = 0;
  while ((g == 0) || (g == 1))
  {
    g = unif_rand();
  }
  while (g > 0)
  {
    g -= group_prob[i];
    i++;
  }
  return i - 1;
}

/**
 * Sample a node. Linear method.
 *
 * @param n_existing Number of existing nodes.
 * @param n_seednode Number of nodes in the initial network.
 * @param pref Sequence of node source/target preference.
 * @param total_pref Total source/target preference of existing nodes.
 * @param sorted_node Node sequence.
 *
 * @return Sampled source/target node.
 */
int sampleNodeLinear(int n_existing, int n_seednode, double *pref,
                      double total_pref, int *sorted_node)
{
  double w = 1;
  int i = 0, j = 0;
  while (w == 1)
  {
    w = unif_rand();
  }
  w *= total_pref;
  while ((w > 0) && (i < n_existing))
  {
    if (i < n_seednode)
    {
      j = sorted_node[i];
    }
    else
    {
      j = i;
    }
    w -= pref[j];
    i++;
  }
  if (w > 0)
  {
    Rprintf("Numerical error! Returning the last node (%d) as the sampled node.\n", n_existing);
    // i = n_existing;
  }
  return j;
}