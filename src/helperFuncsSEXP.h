#ifndef _HELPERFUNCSSEXP_
#define _HELPERFUNCSSEXP_

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

// functions for generating a 2 dimensional connection object
int maxRIntVec(SEXP x);
double maxRDoubleVec(SEXP x);

extern "C" {
SEXP conn2Dim(SEXP dimensions);
};

#endif
