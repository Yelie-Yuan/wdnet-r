#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void netSim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_directed_general_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_directed_general_recip_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rpanet_undirected_general_cpp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _wdnet_directed_rewire_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_findNode_cpp(SEXP, SEXP);
extern SEXP _wdnet_findNode_undirected_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_fx(SEXP, SEXP, SEXP);
extern SEXP _wdnet_hello_world();
extern SEXP _wdnet_nodeStrength_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_simple_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_sampleNode_cpp(SEXP);
extern SEXP _wdnet_sampleNode_simple_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_undirected_rewire_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"netSim",                            (DL_FUNC) &netSim,                            12},
  {"rpanet_directed_general_cpp",       (DL_FUNC) &rpanet_directed_general_cpp,       17},
  {"rpanet_directed_general_recip_cpp", (DL_FUNC) &rpanet_directed_general_recip_cpp, 21},
  {"rpanet_undirected_general_cpp",     (DL_FUNC) &rpanet_undirected_general_cpp,     14},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"_wdnet_directed_rewire_cpp",     (DL_FUNC) &_wdnet_directed_rewire_cpp,     15},
  {"_wdnet_findNode_cpp",            (DL_FUNC) &_wdnet_findNode_cpp,             2},
  {"_wdnet_findNode_undirected_cpp", (DL_FUNC) &_wdnet_findNode_undirected_cpp,  4},
  {"_wdnet_fx",                      (DL_FUNC) &_wdnet_fx,                       3},
  {"_wdnet_hello_world",             (DL_FUNC) &_wdnet_hello_world,              0},
  {"_wdnet_nodeStrength_cpp",        (DL_FUNC) &_wdnet_nodeStrength_cpp,         5},
  {"_wdnet_rpanet_cpp",              (DL_FUNC) &_wdnet_rpanet_cpp,               8},
  {"_wdnet_rpanet_simple_cpp",       (DL_FUNC) &_wdnet_rpanet_simple_cpp,        8},
  {"_wdnet_sampleNode_cpp",          (DL_FUNC) &_wdnet_sampleNode_cpp,           1},
  {"_wdnet_sampleNode_simple_cpp",   (DL_FUNC) &_wdnet_sampleNode_simple_cpp,    4},
  {"_wdnet_undirected_rewire_cpp",   (DL_FUNC) &_wdnet_undirected_rewire_cpp,   10},
  {NULL, NULL, 0}
};

void R_init_wdnet(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}