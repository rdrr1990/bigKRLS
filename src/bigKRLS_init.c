#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bigKRLS_BigCrossProd(SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigDerivMat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigEigen(SEXP, SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigGaussKernel(SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigMultDiag(SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigSolveForc(SEXP, SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigTCrossProd(SEXP, SEXP, SEXP);
extern SEXP bigKRLS_BigTempKernel(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bigKRLS_BigCrossProd",   (DL_FUNC) &bigKRLS_BigCrossProd,   3},
    {"bigKRLS_BigDerivMat",    (DL_FUNC) &bigKRLS_BigDerivMat,    8},
    {"bigKRLS_BigEigen",       (DL_FUNC) &bigKRLS_BigEigen,       4},
    {"bigKRLS_BigGaussKernel", (DL_FUNC) &bigKRLS_BigGaussKernel, 3},
    {"bigKRLS_BigMultDiag",    (DL_FUNC) &bigKRLS_BigMultDiag,    3},
    {"bigKRLS_BigSolveForc",   (DL_FUNC) &bigKRLS_BigSolveForc,   4},
    {"bigKRLS_BigTCrossProd",  (DL_FUNC) &bigKRLS_BigTCrossProd,  3},
    {"bigKRLS_BigTempKernel",  (DL_FUNC) &bigKRLS_BigTempKernel,  4},
    {NULL, NULL, 0}
};

void R_init_bigKRLS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
