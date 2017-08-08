#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ICRanks_PartitioningRankingBlock(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ICRanks_PartitioningRankingLevel(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ICRanks_PartitioningRankingBlock", (DL_FUNC) &_ICRanks_PartitioningRankingBlock, 9},
    {"_ICRanks_PartitioningRankingLevel", (DL_FUNC) &_ICRanks_PartitioningRankingLevel, 5},
    {NULL, NULL, 0}
};

void R_init_ICRanks(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
