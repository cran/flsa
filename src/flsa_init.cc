#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C"
{
    // declaring functions for exporting
    SEXP FLSA(SEXP y); 
    SEXP FLSAexplicitSolution(SEXP solTree, SEXP lambdaVec); 
    SEXP FLSAGeneralMain(SEXP connList, SEXP startValues, SEXP lambdas, SEXP maxSplitSize, SEXP verbose, SEXP thr, SEXP maxGrpNum);
    SEXP FLSAGeneralExplicitSolution(SEXP solObj, SEXP nodes, SEXP lambdas);
    SEXP conn2Dim(SEXP dimensions);
    SEXP FLSATopDown(SEXP, SEXP, SEXP);
    
    void R_init_flsa(DllInfo *info)
    {
        static const R_CallMethodDef CallEntries[] = {
            {"conn2Dim",                    (DL_FUNC) &conn2Dim,                    1},
            {"FLSA",                        (DL_FUNC) &FLSA,                        1},
            {"FLSAexplicitSolution",        (DL_FUNC) &FLSAexplicitSolution,        2},
            {"FLSAGeneralExplicitSolution", (DL_FUNC) &FLSAGeneralExplicitSolution, 3},
            {"FLSAGeneralMain",             (DL_FUNC) &FLSAGeneralMain,             7},
            {"FLSATopDown",                 (DL_FUNC) &FLSATopDown,                 3},
            {NULL, NULL, 0}
        };
        
        R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
        R_useDynamicSymbols(info, FALSE);
    }

    void R_unload_flsa(DllInfo *info)
    {
    
    }
}


