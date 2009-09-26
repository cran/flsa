\name{flsa}
\alias{flsa}
\alias{FLSA}
\alias{flsaGetSolution}
\alias{flsaTopDown}
\title{Fused Lasso Signal Approximator}
\description{
  These functions are the main interface functions for calculating FLSA solutions
}
\usage{
flsa(y, lambda1=0, lambda2=NULL, connListObj = NULL, splitCheckSize=1e+09, verbose=FALSE, thr = 1e-09, maxGrpNum=4*length(y))
flsaTopDown(y, lambda1=0, groups=1:length(y), lambda2=NULL)
flsaGetSolution(solObj, lambda1=0, lambda2=NULL, dim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ response variable; numeric }
  \item{lambda1}{ penalty parameter vector (non-negative) for the absolute values uf the coefficients; numeric }
  \item{lambda2}{ penalty parameter vector (non-negative) for the difference of certain coefficients; numeric }
  \item{groups}{ Return solutions for which the given number of groups is present - solutions found exactly at the breakpoint}	
  \item{connListObj}{an object specifying which differences are to be penalized by lambda2. If NULL, then the dimensionalty of y is being used. If y is a vector, the differences of neighbouring coefficients are penalized. If y is a matrix, differences of neighbouring coefficients in 2 dimensions are being penalized. For more information see \code{\link{connListObj}}}
  \item{splitCheckSize}{a parameter specifying from which size on, groups of variables are not being checked for breaking up; can be used to reduce computation time; may lead to inaccurate results}
  \item{solObj}{Solution object as returned by FLSA if lambda2=NULL}
  \item{dim}{dimensions how the result should be formatted for a specific lambda. Used to format the 2-dimensional FLSA as a matrix in the response. For this, just include the dimensions of \code{y} as \code{dim}}
  \item{verbose}{print status messages during fitting}
  \item{thr}{the error threshold used in the algorithm}
  \item{maxGrpNum}{if every step of the algorithm, a group with a higher number is generated; this limits the number of steps the algorithm can take}
}
\details{
\code{flsa} is the main function for calculate a \code{flsa} fit. If lambda2=NULL, then it returns an object that encodes the whole solution path. Solutions for specific values of lambda1 and lambda2 can then be obtained by using \code{flsaGetSolution}.

\code{flsaTopDown} calculates the solution of the 1-dimensional FLSA, starting at large values of lambda2. If only solutions for large values of lambda2 are needed, this is more efficient.
}
\seealso{\code{\link{connListObj}}}
\examples{
library(flsa)
# generate some artificial data, 1 and 2 dimensional
y <- rnorm(100)
y2Dim = matrix(rnorm(100), ncol=10)

### apply function flsa and get solution directly
lambda2= 0:10/10
res <- flsa(y, lambda2=lambda2)
res2Dim <- flsa(y2Dim, lambda2=lambda2)

### apply the function and get the solution later
resSolObj <- flsa(y, lambda2=NULL)
resSolObjTopDown <- flsaTopDown(y)
resSolObj2Dim <- flsa(y2Dim, lambda2=NULL)

res2 <- flsaGetSolution(resSolObj, lambda2=lambda2)
### here note that the solution object does not store that the input was 2 dimensional
### therefore, res2Dim does not give out the solution as a 2 dimensional matrix (unlike the direct version above)
res2Dim2 <- flsaGetSolution(resSolObj2Dim, lambda2=lambda2)

}
\author{Holger Hoefling}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{multivariate}





