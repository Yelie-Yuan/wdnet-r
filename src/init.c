#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void netSim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_binary_directed_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_binary_undirected_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_naive_directed_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_naive_undirected_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _wdnet_fillWeight_cpp(SEXP, SEXP, SEXP);
extern SEXP _wdnet_findNode_cpp(SEXP, SEXP);
extern SEXP _wdnet_findNode_undirected_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_fx(SEXP, SEXP, SEXP);
extern SEXP _wdnet_hello_world();
extern SEXP _wdnet_nodeStrength_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rewire_directed_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rewire_undirected_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_sampleNode_cpp(SEXP);

static const R_CMethodDef CEntries[] = {
    {"netSim",                       (DL_FUNC) &netSim,                       12},
    {"rpanet_binary_directed_cpp",   (DL_FUNC) &rpanet_binary_directed_cpp,   29},
    {"rpanet_binary_undirected_cpp", (DL_FUNC) &rpanet_binary_undirected_cpp, 17},
    {"rpanet_naive_directed_cpp",    (DL_FUNC) &rpanet_naive_directed_cpp,    29},
    {"rpanet_naive_undirected_cpp",  (DL_FUNC) &rpanet_naive_undirected_cpp,  17},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_wdnet_fillWeight_cpp",          (DL_FUNC) &_wdnet_fillWeight_cpp,           3},
    {"_wdnet_findNode_cpp",            (DL_FUNC) &_wdnet_findNode_cpp,             2},
    {"_wdnet_findNode_undirected_cpp", (DL_FUNC) &_wdnet_findNode_undirected_cpp,  4},
    {"_wdnet_fx",                      (DL_FUNC) &_wdnet_fx,                       3},
    {"_wdnet_hello_world",             (DL_FUNC) &_wdnet_hello_world,              0},
    {"_wdnet_nodeStrength_cpp",        (DL_FUNC) &_wdnet_nodeStrength_cpp,         5},
    {"_wdnet_rewire_directed_cpp",     (DL_FUNC) &_wdnet_rewire_directed_cpp,     11},
    {"_wdnet_rewire_undirected_cpp",   (DL_FUNC) &_wdnet_rewire_undirected_cpp,   10},
    {"_wdnet_rpanet_cpp",              (DL_FUNC) &_wdnet_rpanet_cpp,               8},
    {"_wdnet_sampleNode_cpp",          (DL_FUNC) &_wdnet_sampleNode_cpp,           1},
    {NULL, NULL, 0}
};

void R_init_wdnet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}