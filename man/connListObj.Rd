\name{is.connListObj}
\alias{ConnListObj}
\alias{is.connListObj}
\alias{connListObj}
\title{Connection List Objects}
\description{
  Describes the makeup of a connection list object
}
\usage{
is.connListObj(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{the object to be tested}
}
\details{
A connection list object can be used to specifiy which differences in
\code{fusedlasso} or \code{flsa} functions are to be penalized. Here, it
is assumed that the n coefficients in the model are numbered from 0 to
n-1. The connection list object is a list of length n with each element
corresponding to one of the coefficients. The i-th element of the list
here corresponds to coefficient with number i-1. Each element of the
list is a vector of integers, naming the numbers of the coefficients to
which the coefficient corresponding to the current list element is
linked (i.e. the difference of the two coefficients is being penalized).
I.e., assume that value $j$ is a member of the list under list element
$i$. Then this means that coeffient $i-1$ and coefficient $j$ are being
penalized. To understand this, consider that R-lists when viewed in
C-code are being numbered starting with 0, not 1 and note that all
computation is being done in C-code.

Furthermore, the connection list object has class \code{connListObj}.

Also note that the vectors in the list are of type \code{integer} not \code{numeric}. An empty vector should be set to \code{NULL}.
}
\seealso{\code{\link{connListObj}}}
\author{Holger Hoefling}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{multivariate}
\examples{

connList <- vector("list", 4)
y <- 1:4

class(connList) = "connListObj"
connList[[1]] = as.integer(c(1,2))
connList[[2]] = as.integer(c(0,3))
connList[[3]] = as.integer(c(3,0))
connList[[4]] = as.integer(c(2,1))
names(connList) <- as.character(0:3) ## not necessary, just for illustration

res <- flsa(y, connListObj=connList)
res2 <- flsa(matrix(y, nrow=2))

res$BeginLambda
res2$BeginLambda
flsaGetSolution(res, lambda2=c(0, 0.5, 1))
flsaGetSolution(res2, lambda2=c(0, 0.5, 1))
}