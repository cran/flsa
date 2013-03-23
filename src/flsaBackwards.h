#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <Rinternals.h>
#include <R.h>

using namespace std;

struct GroupInfo {
    int left;
    int right;
    int breakPos;
    double breakTime;
};

// the class that keeps the datastructure
class FLSABackwards {
    int n;
    double* y;
    double* beta;
    double* betaDeriv;
    double* updateLambdaBeta;
    double* tau;
    double* tauDeriv;
    double* updateLambdaTau;
   
	SEXP R_solution;
 
    double* lambdas;
    double curLambda;

	int* isBreakpointVec;   
 
    int numSolutions;
    int curNumSolutions;
    double* solution;
    
    set<int> solGroups;
    set<double> solLambdas;
    multimap<double, GroupInfo> groups;

    bool algRun;
    
    // utility functions that are needed; starts searching at from up to to included
    GroupInfo findBreakpoint(int from, int to);
    
    // calculate the break time for position pos
    double getBreakTime(int pos);
    
    // calculate the betaDerivative of the group
    double calcBetaDeriv(int from, int to);
    
    // update the tau in the group
    void updateTau(int from, int to, double lambda);
    
    // calculate the deriative of tau
    void calcTauDeriv(int from, int to, double betaDeriv);
    
    // saves the current beta as the solution (after adjusting using the derivative)
    void saveCurBetaAsSolution(double lambdaHere, bool isBreakpoint);

    pair<GroupInfo, GroupInfo> splitGroup(GroupInfo group);

public:
    /*
    The constructor for the class; 
    y is the vector with the input data
    n is the number of elements in the vector y
    numGroups is a vector with the number of groups desired in the solution
    numSolutions is the length of the numGroups vector
    */
    FLSABackwards(SEXP R_y, SEXP R_groups, SEXP R_lambdas); // Constructor
    
    // run the algorithm
    void runAlgorithm();
    
    // get the solution that has been saved; returns an SEXP matrix
    void allocateSolutions();
	SEXP returnSolutions();
    
    void printResults(ostream &out);
    
    void printVector(ostream &out, double* x, int len);
    
    void printGroups(ostream &out);
    
    void printSolGroups(ostream &out);
};

