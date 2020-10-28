// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
#include <R.h>

void netSim(double* alpha_ptr,
	    double* beta_ptr,
	    double* gamma_ptr,
	    double* xi_ptr,
	    double* delta_in_ptr,
	    double* delta_out_ptr,
	    int* target_n_edges_double_ptr,
	    //output
	    int *in_deg_counts,
	    int *out_deg_counts,
	    int *edge_start,
	    int *edge_end,
	    int *evolution
	    )
{
  double alpha = *alpha_ptr, beta = *beta_ptr, gamma = *gamma_ptr,
    xi = *xi_ptr, delta_in = *delta_in_ptr, delta_out = *delta_out_ptr;

    
    int target_n_edges = *target_n_edges_double_ptr;
    int current_n_edges;
    int current_n_nodes;
    int current_scenario;
    /* int* in_deg_counts = (int*) malloc(sizeof(int) *target_n_edges*2); */
    /* int* out_deg_counts = (int*) malloc(sizeof(int) *target_n_edges*2); */
    /* int* edge_start = (int*) malloc(sizeof(int) *target_n_edges); */
    /* int* edge_end = (int*) malloc(sizeof(int) *target_n_edges); */
    /* int* evolution = (int*) malloc(sizeof(int) *target_n_edges); */
    int this_edge_start;
    int this_edge_end;
    
    edge_start[0] = 0;
    edge_end[0] = 1;
    evolution[0] = 1;
    current_n_edges = 1;
    current_n_nodes = 2;
    current_scenario = 1;
    in_deg_counts[0] = 0;
    out_deg_counts[0] = 1;
    in_deg_counts[1] = 1;
    out_deg_counts[1] = 0;
    
    double u,v,w;
    
    // srand((time(NULL) & 0xFFFF) | (getpid() << 16));
    GetRNGstate();
    
    while (current_n_edges < target_n_edges)
    {
        // u = (double)rand() / (double)RAND_MAX;
        // u = (double) ( rand() % 100000 / (double) 100000 );
        u = unif_rand();
        
        if (u <= alpha)
        {
            // Scenario alpha
            // v = (double) ( rand() % 100000 / (double) 100000 ) * (current_n_edges + current_n_nodes*delta_in);
	    v = unif_rand() * (current_n_edges + current_n_nodes*delta_in);
            if (v < current_n_edges)
            {
                this_edge_end = edge_end[(int) ceil(v)-1];
            }
            else
            {
                this_edge_end = (int) ceil((v-current_n_edges)/delta_in) -1;
            }
            this_edge_start = current_n_nodes;
            in_deg_counts[this_edge_end] ++;
            out_deg_counts[this_edge_start] = 1;
            in_deg_counts[this_edge_start] = 0;
            current_n_nodes++;
            current_scenario = 1;
        }
        
        else if ((u > alpha) & (u <= alpha+beta))
        {
            // Scenario beta
            // v = (double) ( rand() % 100000 / (double) 100000 )* (current_n_edges + current_n_nodes*delta_in);
	    v = unif_rand() * (current_n_edges + current_n_nodes*delta_in);
            if (v < current_n_edges)
            {
                this_edge_end = edge_end[(int) ceil(v)-1];
            }
            else
            {
                this_edge_end = (int) ceil((v-current_n_edges)/delta_in) -1;
            }
            
            // w = (double) ( rand() % 100000 / (double) 100000 )* (current_n_edges + current_n_nodes*delta_out);
	    w = unif_rand() * (current_n_edges + current_n_nodes*delta_out);
            if (w < current_n_edges)
            {
                this_edge_start = edge_start[(int) ceil(w)-1];
            }
            else
            {
                this_edge_start = (int) ceil((w-current_n_edges)/delta_out) -1;
            }
            in_deg_counts[this_edge_end] ++;
            out_deg_counts[this_edge_start] ++;
            current_scenario = 2;
        }
        
        else if ((u > alpha+beta) & (u <= alpha+beta+gamma))
        {
            // Scenario gamma
            // w = (double) ( rand() % 100000 / (double) 100000 )* (current_n_edges + current_n_nodes*delta_out);
	    w = unif_rand() * (current_n_edges + current_n_nodes*delta_out);
            if (w < current_n_edges)
            {
                this_edge_start = edge_start[(int) ceil(w)-1];
            }
            else
            {
                this_edge_start = (int) ceil((w-current_n_edges)/delta_out) -1;
            }
            this_edge_end = current_n_nodes;
            in_deg_counts[this_edge_end] = 1;
            out_deg_counts[this_edge_end] = 0;
            out_deg_counts[this_edge_start] ++;
            current_n_nodes++;
            current_scenario = 3;
        }
        else if ((u > alpha+beta+gamma) & (u <= alpha+beta+gamma+xi))
        {
            // Scenario xi
            this_edge_start = (int) current_n_nodes;
            in_deg_counts[this_edge_start] = 0;
            out_deg_counts[this_edge_start] = 1;
            current_n_nodes++;
            this_edge_end = (int) current_n_nodes;
            in_deg_counts[this_edge_end] = 1;
            out_deg_counts[this_edge_end] = 0;
            current_n_nodes++;
            current_scenario = 4;
        }
        
        else // if (u > alpha+beta+gamma+xi & u <= 1)
        {
            // Scenario loop
            this_edge_start = (int) current_n_nodes;
            this_edge_end = (int) current_n_nodes;
            in_deg_counts[this_edge_end] = 1;
            out_deg_counts[this_edge_end] = 1;
            current_n_nodes++;
            current_scenario = 5;
        }
        
        evolution[current_n_edges] = current_scenario;
        edge_end[current_n_edges] = this_edge_end;
        edge_start[current_n_edges] = this_edge_start;
        current_n_edges ++;
    }
    PutRNGstate();
}

/* // Registering */
/* #include <stdlib.h> */
/* #include <R_ext/Rdynload.h> */
/* #include <R_ext/Visibility.h>  // optional */

/* #define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n} */

/* static const R_CallMethodDef R_CallDef[] = { */
/*     CALLDEF(netSim, 12), */
/*     {NULL, NULL, 0} */
/* }; */

/* void */
/* attribute_visible  // optional */
/* R_init_wdnet(DllInfo *dll) */
/* { */
/*     R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL); */
/*     R_useDynamicSymbols(dll, FALSE); */
/*     R_forceSymbols(dll, TRUE); */
/* } */
