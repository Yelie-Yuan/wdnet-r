#pragma once


typedef double (*funcPtrUnd)(double x);

typedef double (*funcPtrD)(double x, double y);


double prefFuncD(double outs, double ins, double *params);

double prefFuncUnd(double strength, double *params);

int sampleGroup(double *group_prob);

int sampleNodeLinear(int n_existing, int n_seednode, double *pref,
                     double total_pref, int *sorted_node);
