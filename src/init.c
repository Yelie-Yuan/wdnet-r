#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void netSim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _wdnet_dprewire_directed_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_dprewire_undirected_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_fill_weight_cpp(SEXP, SEXP, SEXP);
extern SEXP _wdnet_find_node_cpp(SEXP, SEXP);
extern SEXP _wdnet_find_node_undirected_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_fx(SEXP, SEXP, SEXP);
extern SEXP _wdnet_hello_world();
extern SEXP _wdnet_node_strength_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_binary_directed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_binary_undirected_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_naive_directed_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_naive_undirected_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_nodelist_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_sample_node_cpp(SEXP);

static const R_CMethodDef CEntries[] = {
    {"netSim", (DL_FUNC) &netSim, 12},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_wdnet_dprewire_directed_cpp",        (DL_FUNC) &_wdnet_dprewire_directed_cpp,        11},
    {"_wdnet_dprewire_undirected_cpp",      (DL_FUNC) &_wdnet_dprewire_undirected_cpp,      10},
    {"_wdnet_fill_weight_cpp",              (DL_FUNC) &_wdnet_fill_weight_cpp,               3},
    {"_wdnet_find_node_cpp",                (DL_FUNC) &_wdnet_find_node_cpp,                 2},
    {"_wdnet_find_node_undirected_cpp",     (DL_FUNC) &_wdnet_find_node_undirected_cpp,      4},
    {"_wdnet_fx",                           (DL_FUNC) &_wdnet_fx,                            3},
    {"_wdnet_hello_world",                  (DL_FUNC) &_wdnet_hello_world,                   0},
    {"_wdnet_node_strength_cpp",            (DL_FUNC) &_wdnet_node_strength_cpp,             5},
    {"_wdnet_rpanet_binary_directed",       (DL_FUNC) &_wdnet_rpanet_binary_directed,       15},
    {"_wdnet_rpanet_binary_undirected_cpp", (DL_FUNC) &_wdnet_rpanet_binary_undirected_cpp, 11},
    {"_wdnet_rpanet_naive_directed_cpp",    (DL_FUNC) &_wdnet_rpanet_naive_directed_cpp,    15},
    {"_wdnet_rpanet_naive_undirected_cpp",  (DL_FUNC) &_wdnet_rpanet_naive_undirected_cpp,  11},
    {"_wdnet_rpanet_nodelist_cpp",          (DL_FUNC) &_wdnet_rpanet_nodelist_cpp,           8},
    {"_wdnet_sample_node_cpp",              (DL_FUNC) &_wdnet_sample_node_cpp,               1},
    {NULL, NULL, 0}
};

void R_init_wdnet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
