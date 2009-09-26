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
A connection list object can be used to specifiy which differences in \code{fusedlasso} or \code{flsa} functions are to be penalized. Here, it is assumed that the n coefficients in the model are numbered from 1 to n. The connection list object is a list of length n with each element corresponding to one of the coefficients. The i-th element of the list is a vector of integers. If j is an element of the i-th element of the connection list, then the difference between the i-th and j-th coefficients are being penalized. It is assumed that this is symmetric, i.e. that if j is an element of the i-th vector, then also i is an element of the j-th vector. Furthermore, the connection list object has class \code{connListObj}.

Also note that the vectors in the list are of type \code{integer} not \code{numeric}. An empty vector should be set to \code{NULL}.
}
\seealso{\code{\link{connListObj}}}
\author{Holger Hoefling}
\keyword{regression}% at least one, from doc/KEYWORDS
\keyword{multivariate}