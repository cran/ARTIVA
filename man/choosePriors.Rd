\name{choosePriors}
\alias{choosePriors}
\title{
Function to plot an overview of possible priors for the number of changepoints and/or edges. 
}
\description{
   Plots an overview of some possible priors used in the
   \code{\link{ARTIVAsubnet}} function  for the number of changepoints
   (resp. incoming edges) according to a given number of maximum
   changepoints \code{maxCP} (resp. incoming Edges \code{maxPred}) when
   parameters (\code{alphaCP}, \code{betaCP} for the CPs or
   \code{alphaEdges},  \code{betaEdges} for the edges) in function
   \code{\link{ARTIVAsubnet}} are set to default (alpha=1, beta=0.5).
   In the \code{\link{ARTIVAsubnet}} procedure, the number of CPs
   (respectively the number of incoming edges) is sampled from a
   truncated Poisson with mean lambda, where lambda is  drawn from an
   Inverse Gamma distribution (alpha, beta), see Lebre et al. (2010) for
   more details.
}
\usage{
choosePriors(kmax)
}
\arguments{
\item{kmax}{Maximum number of changepoints or incomins Edges (parents)}
}
\value{
NULL, a graph is plotted.
}
\references{
  S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais (2010)
  Statistical inference of the time-varying structure of gene-regulation
  networks, BMC Systems Biology, 4:130.
}
\author{
S. Lebre and G. Lelandais.
}

\seealso{\code{\link{ARTIVAsubnet}}, \code{\link{ARTIVAnet}} 
}
\examples{
# See some prior probabbility density when the maximal number of
# CPs/parents is equal to 5
choosePriors(kmax=5)
}
\keyword{util}