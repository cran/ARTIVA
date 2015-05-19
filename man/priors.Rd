\name{priors}
\alias{priors}
\docType{data}
\title{
 Set of possible priors for the number of changepoints (CPs) or incoming edges. 
}
\description{
   Set of possible priors used in the \code{\link{ARTIVAsubnet}} function for the
   number of changepoints (resp. incoming edges) according to a given number
   of maximum changepoints \code{maxCP} (resp. incoming Edges \code{maxPred}) 
   when parameters (\code{alphaCP}, \code{betaCP} for the CPs or
   \code{alphaEdges},  \code{betaEdges} for the edges) in function
   \code{\link{ARTIVAsubnet}} are set to default (alpha=1, beta=0.5).
   In the \code{\link{ARTIVAsubnet}} procedure, the number of CPs
   (respectively the number of incoming edges) is sampled from a
   truncated Poisson with mean lambda, where lambda is  drawn from an
   Inverse Gamma distribution (alpha, beta), see Lebre et al. (2010) for
   more details.
}
\usage{data(priors)}
\format{
  A matrix of 96 rows by 44 columns (kmax, alpha, beta and the probability for k=0 to 40
  according to the chosen values of kmax, alpha and beta).
}
\references{
 Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 2010, 4:130.
}
\seealso{\code{\link{choosePriors}},\code{\link{ARTIVAsubnet}}, \code{\link{ARTIVAnet}} 
}
\examples{
# See some prior probability density when the maximal number of
# CPs/parents is equal to 5
data(priors)
choosePriors(kmax=5,priors)
}
\keyword{datasets}
