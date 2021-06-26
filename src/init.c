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
extern SEXP _wdnet_findNode_cpp(SEXP, SEXP, SEXP);
extern SEXP _wdnet_fx(SEXP, SEXP, SEXP);
extern SEXP _wdnet_hello_world();
extern SEXP _wdnet_nodeStrength_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rewire_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_rpanet_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _wdnet_sampleNode_cpp(SEXP);

static const R_CMethodDef CEntries[] = {
    {"netSim", (DL_FUNC) &netSim, 12},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_wdnet_findNode_cpp",     (DL_FUNC) &_wdnet_findNode_cpp,     3},
    {"_wdnet_fx",               (DL_FUNC) &_wdnet_fx,               3},
    {"_wdnet_hello_world",      (DL_FUNC) &_wdnet_hello_world,      0},
    {"_wdnet_nodeStrength_cpp", (DL_FUNC) &_wdnet_nodeStrength_cpp, 5},
    {"_wdnet_rewire_cpp",       (DL_FUNC) &_wdnet_rewire_cpp,       5},
    {"_wdnet_rpanet_cpp",       (DL_FUNC) &_wdnet_rpanet_cpp,       7},
    {"_wdnet_sampleNode_cpp",   (DL_FUNC) &_wdnet_sampleNode_cpp,   1},
    {NULL, NULL, 0}
};

void R_init_wdnet(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}