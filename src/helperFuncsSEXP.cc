#include "helperFuncsSEXP.h"
#include "GeneralFunctions.h"
#include <list>

using namespace std;

list<int> pointConn(int r, int c, int dimRow, int dimCol, int counter)
{
    list<int> foo;
    if(c>0) // connection to the left
    {
        foo.push_front(counter-dimRow);
    }
    if(c<(dimCol-1))// connection to the right
    {
        foo.push_front(counter+dimRow);
    }
    if(r>0) // connection above
    {
        foo.push_front(counter-1);
    }
    if(r<(dimRow-1))
    {
        foo.push_front(counter+1);
    }
    return(foo);
}



// generates the list of connections for the 2-dimensional fused lasso
// done in C++ to increase speed
// all variables are assumed checked for right format
// dimensions is a vector of two integers
// the nodes will be numbered starting at 0 along the columns (like an R matrix)
SEXP conn2Dim(SEXP dimensions)
{
    SEXP conn, bar;
    int dimRow = INTEGER(dimensions)[0];
    int dimCol = INTEGER(dimensions)[1];
    PROTECT(conn = allocVector(VECSXP, dimRow*dimCol));
    
    list<int> foo; // used to temporarily save the connections
    int counter = 0;
    // go through all gridpoints
    for(int j=0; j<dimCol; ++j)
    {
        for(int i=0; i<dimRow; ++i)
        {
            // get the connections of the current point
            foo=pointConn(i,j,dimRow,dimCol,counter);
            
            // now copy it into an R vector
            PROTECT(bar= allocVector(INTSXP, foo.size()));
            for(int k=0; k<LENGTH(bar); ++k)
            {
                INTEGER(bar)[k] = foo.front();
                foo.pop_front();
            }
            SET_VECTOR_ELT(conn,counter,bar);
            UNPROTECT(1);
            ++counter;
        }
    }
    
    UNPROTECT(1);
    return(conn);
};

int maxRIntVec(SEXP x)
{
//    Rprintf("Started max\n");
    int max=0;
    int len = LENGTH(x);
    int *xVec = INTEGER(x);
    for(int i =0; i< len; ++i)
    {
        if(xVec[i]>max)
        {
            max=xVec[i];
        }
    }
    return(max);
}

double maxRDoubleVec(SEXP x)
{
//    Rprintf("Started max\n");
    double max=-infinite;
    int len = LENGTH(x);
    double *xVec = REAL(x);
    for(int i =0; i< len; ++i)
    {
        if(xVec[i]>max)
        {
            max=xVec[i];
        }
    }
    return(max);
}
