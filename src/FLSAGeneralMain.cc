#include <iostream>
#include <fstream>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "PenaltyGraph.h"
#include "MaxFlowGraph.h"
#include "Groups.h"
#include <set>
#include "FLSAGeneral.h"
#include "helperFuncsSEXP.h"

using namespace std;

// find the maximum of a vector of integers


extern "C" {

// main function that starts all the related FLSAGeneral classes

SEXP FLSAGeneralMain(SEXP connList, SEXP startValues, SEXP lambdas, SEXP maxSplitSize, SEXP verbose, SEXP thr, SEXP maxGrpNum)
{
    // find the highest nodeNumber
    int highNode = LENGTH(connList)-1;
    // cout << "Highest node is" << highNode << endl;
    SEXP sol;
    double highestLambda=infinite;
    
    if(IS_NUMERIC(lambdas))
    {
        highestLambda = maxRDoubleVec(lambdas);
    }
    
    FLSAGeneral FLSAGeneralObj(highNode, connList, startValues, maxSplitSize, verbose, thr, maxGrpNum, highestLambda);

    if(!IS_NUMERIC(lambdas))
    {
        sol=FLSAGeneralObj.solutionGraph();
    }
    else
    {
        SEXP nodes;
        PROTECT(nodes = allocVector(INTSXP, highNode+1));
        for(int i=0; i<=highNode; ++i) {
            INTEGER(nodes)[i]=i;
        }
        sol =FLSAGeneralObj.solution(nodes, lambdas);
        UNPROTECT(1);
    }
    return(sol);
};


SEXP FLSAGeneralExplicitSolution(SEXP solObj, SEXP nodes, SEXP lambdas)
{
    Groups groups(solObj);
    return(groups.solution(nodes, lambdas));
}

}; // end of extern C

