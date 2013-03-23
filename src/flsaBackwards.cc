#include "flsaBackwards.h"
#include<limits>
#include<algorithm>

// utility function for finding the maximum


GroupInfo FLSABackwards::findBreakpoint(int from, int to) {
    GroupInfo res;
    res.left = from;
    res.right = to;
    double foo;
    res.breakPos=-1;
    res.breakTime = 0;
    for(int i=from; i<to; ++i) {
        foo = getBreakTime(i);
        if(foo>=res.breakTime) {
            res.breakTime=foo;
            res.breakPos=i;
        }
    }
    return(res);
}

// get the breaktime for position pos; does not check that pos is within 0 and n-2
double FLSABackwards::getBreakTime(int pos) {
    // when does tau hit the upper boundary
    double a = updateLambdaTau[pos]+(updateLambdaTau[pos]-tau[pos])/(tauDeriv[pos]-1);
    // and the lower boundary
    double b =  updateLambdaTau[pos]+(-updateLambdaTau[pos]-tau[pos])/(1+tauDeriv[pos]);
    // will never hit the boundaries
    if(a > updateLambdaTau[pos]) {a=0;}
    if(b > updateLambdaTau[pos]) {b=0;}
    // return first hit
    return(max(a,b));
}

// calculates the derivative of beta
double FLSABackwards::calcBetaDeriv(int from, int to) {
    int len = to-from+1;
    double boundary = 0;
    if(from>0) { // not the left boundary
        if(tau[from-1] >0) {
            boundary +=1;
        }
        else {
            boundary -=1;
        }
    }
    if(to<n-1) { // also not the right boundary
        if(tau[to] > 0) {
            boundary -= 1;
        }
        else {
            boundary += 1;
        }
    }
    return(boundary/len);
}

// update the tau in the group
void FLSABackwards::updateTau(int from, int to, double lambda) {
    for(int i=from; i<to; ++i) {
        tau[i] = tau[i] + tauDeriv[i] * (lambda-updateLambdaTau[i]);
        updateLambdaTau[i] = lambda;
    }
}

// calculate the deriative of tau
void FLSABackwards::calcTauDeriv(int from, int to, double betaDeriv) {
    double leftDeriv;
	bool violated = false;
    if(from==0) {
        leftDeriv = 0;
    }
    else if(tau[from-1] < 0) {
        leftDeriv = -1;
    }
    else {
        leftDeriv = 1;
    }

    for(int i=from; i<to; ++i) {
        leftDeriv-=betaDeriv;
        tauDeriv[i]=leftDeriv;
   		if(tauDeriv[i]<-1 || tauDeriv[i]>1) {
			violated = true; 
		}
	}
}

void FLSABackwards::saveCurBetaAsSolution(double lambdaHere, bool isBreakpoint) {
    multimap<double, GroupInfo>::iterator mapIt;
    // walk through all groups that are currently saved (are all in the map)
    // and save the beta
    for(mapIt=groups.begin(); mapIt!=groups.end(); ++mapIt) {
        int left = mapIt->second.left;
        int right = mapIt->second.right;
        double betaHere = beta[left] + betaDeriv[left] * (lambdaHere - updateLambdaBeta[left]);
        for(int i=left; i<=right; ++i) {
            solution[curNumSolutions*n + i] = betaHere;
        }
    }
    // also save the lambda
    lambdas[curNumSolutions]=lambdaHere;
	isBreakpointVec[curNumSolutions] = isBreakpoint;
    curNumSolutions++;
}




pair<GroupInfo, GroupInfo> FLSABackwards::splitGroup(GroupInfo group) {
    double betaHere = beta[group.left] + betaDeriv[group.left] * (curLambda-updateLambdaBeta[group.left]);
    updateTau(group.left, group.right, curLambda);
    
    // get the beta movement of the two new groups
    double leftBetaDeriv = calcBetaDeriv(group.left, group.breakPos);
    double rightBetaDeriv = calcBetaDeriv(group.breakPos+1, group.right);
    
    // update the derivative of tau in these two groups
    calcTauDeriv(group.left, group.breakPos, leftBetaDeriv);
    calcTauDeriv(group.breakPos+1, group.right, rightBetaDeriv);
    
    // find the splitpoint of each of the new groups
    GroupInfo left = findBreakpoint(group.left, group.breakPos);
    GroupInfo right = findBreakpoint(group.breakPos+1, group.right);
    
    // set the new betas and the new derivatives
    beta[left.left] = betaHere;
    betaDeriv[left.left] = leftBetaDeriv;
    updateLambdaBeta[left.left] = curLambda;
    beta[right.left] = betaHere;
    betaDeriv[right.left] = rightBetaDeriv;
    updateLambdaBeta[right.left] = curLambda;
    
    // add the new groups to the list and remove the old one
    pair<GroupInfo, GroupInfo> res;
    res.first = left;
    res.second = right;
    return(res);
}


// Constructor
FLSABackwards::FLSABackwards(SEXP R_y, SEXP R_groups, SEXP R_lambdas):solGroups(),groups()
{
	// first initialize the variables that we need
    n = LENGTH(R_y);
    y = REAL(R_y); // only to be used for reading, not writing
    beta = new double[n];
    betaDeriv = new double[n];
    updateLambdaBeta = new double[n];
    tau = new double[n-1];
    tauDeriv = new double[n-1];
    updateLambdaTau = new double[n-1];
    
	// use groups, and/or lambdas
    // now find out how many groups are needed
    for(int i=0; i< LENGTH(R_groups); ++i) {
        int group = INTEGER(R_groups)[i];
        if(group>=1 && group<=n) {
            solGroups.insert(group);
        }
    }

	// now the regular lambdas
	for(int i=0; i< LENGTH(R_lambdas); ++i) {
		double lambdaHere = REAL(R_lambdas)[i];
		if(lambdaHere >=0) {
			solLambdas.insert(lambdaHere);
		}
	}
    numSolutions = solGroups.size() + solLambdas.size();
    curNumSolutions = 0;
    
    // create the memory space for the solutions; these don't need to be initialized
    allocateSolutions();
    algRun = false;
    
    // initialize beta, betaDeriv, tau, tauderiv
    GroupInfo group;
    group.left=0;
    group.right=n-1;
    
    // calculate the mean
    double mean = 0;
    for(int i=0; i<n; ++i)
    {
        mean += y[i];
    }
    mean /= n;
    beta[0] = mean;
    betaDeriv[0] = 0;
    
    // now set tau, tauDeriv
    tau[0]=y[0]-mean;
    group.breakTime = fabs(tau[0]);
    group.breakPos = 0;
    tauDeriv[0]=0;
    for(int i=1; i<n-1; ++i) {
        tau[i] = tau[i-1]+y[i]-mean;
        tauDeriv[i]=0;
        if(group.breakTime < fabs(tau[i])) {
            group.breakTime = fabs(tau[i]);
            group.breakPos = i;
        }
    }
    // insert the group into the map
    groups.insert(pair<double, GroupInfo>(group.breakTime, group));
   
    // now initialize the lambda values to 1 higher than the first break point
	// or at the highest lambda
	double highestLambda;
	if(solLambdas.size()>0) {
		highestLambda = *(solLambdas.rbegin());
	}
	else {
	 	highestLambda = 0;
	}
	curLambda = max(group.breakTime, highestLambda)+1;
    // initialize the vectors that store when beta and tau were last updated
    updateLambdaBeta[0] = curLambda;
    for(int i=0; i<n-1; ++i) {
        updateLambdaTau[i]=curLambda;
    }
}

// run the algorithm
void FLSABackwards::runAlgorithm()
{
 	double nextBreakpoint;
	double nextLambda;
    int counter=1;
    // while there is a next solution to calculate
    while(solGroups.size()>0 || solLambdas.size()>0) {
        // find the current group (the one with the highest value of lambda
        multimap<double, GroupInfo>::reverse_iterator nextGroup = groups.rbegin();
        nextBreakpoint = nextGroup->first;
		
		// first save the solutions for specific lambda values
        // check if a solution should be saved and if yes, erase it from solGroups afterwards
        
		// now save the solutions for specific lambda values
		while(solLambdas.size()>0 && *(solLambdas.rbegin())>=nextBreakpoint) {
		 	nextLambda = *(solLambdas.rbegin());
			saveCurBetaAsSolution(nextLambda, false);
			solLambdas.erase(nextLambda);
		}
		
		if(solGroups.find(groups.size())!=solGroups.end()) {
			saveCurBetaAsSolution(nextBreakpoint, true);
			solGroups.erase(groups.size());
	    }
		
		// set the current lambda value
		curLambda = nextBreakpoint;
		// no more work to do
		if((int) groups.size() == n) {
			break;
		}

        // now split the current group at the right position
        // also adds the new groups to the groups map
        pair<GroupInfo, GroupInfo> newGroups = splitGroup(nextGroup->second);
        // delete the current group and add two new ones
        multimap<double, GroupInfo>::iterator eraseElement;
        if(groups.size()==1) {
            eraseElement = groups.begin(); // only one element there
        }
        else {
            eraseElement = (++nextGroup).base(); // gives back the last element as a regular iterator
        }
        groups.erase(eraseElement);
        groups.insert(pair<double, GroupInfo>(newGroups.first.breakTime, newGroups.first));
        groups.insert(pair<double, GroupInfo>(newGroups.second.breakTime, newGroups.second));
        
        counter++;
    }
    algRun = true;
}
    
// get the solution that has been saved; returns an SEXP matrix
void FLSABackwards::allocateSolutions()
{
    PROTECT(R_solution = allocVector(VECSXP, 3));

    // set the names of the components
    SEXP names = getAttrib(R_solution, R_NamesSymbol);
    names = allocVector(STRSXP,3);
    SET_STRING_ELT(names, 0, mkChar("Solution"));
    SET_STRING_ELT(names, 1, mkChar("Lambdas"));
	SET_STRING_ELT(names, 2, mkChar("isBreakpoint"));
    setAttrib(R_solution, R_NamesSymbol, names);
    
    // set the matrix with the solutions 
    // and the vector with the values of lambda at the breakpoints
    SET_VECTOR_ELT(R_solution, 0, allocMatrix(REALSXP,n,numSolutions));
    SET_VECTOR_ELT(R_solution, 1, allocVector(REALSXP,numSolutions));
    SET_VECTOR_ELT(R_solution, 2, allocVector(LGLSXP,numSolutions));

    // copy the data over
    solution = REAL(VECTOR_ELT(R_solution,0));
    lambdas = REAL(VECTOR_ELT(R_solution,1));
    isBreakpointVec = LOGICAL(VECTOR_ELT(R_solution,2));
}


SEXP FLSABackwards::returnSolutions() {
	UNPROTECT(1);
	return(R_solution);
}


void FLSABackwards::printResults(ostream &out) {
    out << "--------------------------------------------------" << endl;
    out << "y:";
    printVector(out,y,n);
    
    out << "Beta:";
    printVector(out,beta,n);
    
    out << "BetaDeriv:";
    printVector(out, betaDeriv, n);
    
    out << "UpdateLambdaBeta:";
    printVector(out, updateLambdaBeta, n);
    
    out << "Tau:";
    printVector(out, tau,n);
    
    out << "TauDeriv:";
    printVector(out, tauDeriv, n);
    
    out << "UpdateLambdaTau:";
    printVector(out, updateLambdaTau,n);
    out << "----------------------------------------------------" << endl;
}

void FLSABackwards::printVector(ostream &out, double* x, int len) {
    for(int i=0; i< len; ++i) {
        out << x[i] << " ";
    }
    out << endl;
}


void FLSABackwards::printGroups(ostream &out) {
    multimap<double, GroupInfo>::iterator mapIt;
    for(mapIt=groups.begin(); mapIt!=groups.end(); ++mapIt) {
        out << "Left: " << mapIt->second.left;
        out << "Right: " << mapIt->second.right;
        out << "BreakPos: " << mapIt->second.breakPos;
        out << "BreakTime: " << mapIt->second.breakTime << endl;
    }
}


void FLSABackwards::printSolGroups(ostream &out) {
    set<int>::iterator setIt=solGroups.begin();
    
    for(;setIt!=solGroups.end(); ++setIt) {
        out << *setIt << " ";
    }
    out << endl;
}


// interface function to R

extern "C" {
SEXP FLSATopDown(SEXP R_y, SEXP R_groups, SEXP R_lambdas) {
    FLSABackwards flsa(R_y, R_groups, R_lambdas);
    
    flsa.runAlgorithm();
    
    SEXP R_return = flsa.returnSolutions();
    
    return(R_return);
}
}

