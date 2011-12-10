\name{ARTIVA-package}
\alias{ARTIVA-package}
\alias{ARTIVA}
\docType{package}
\title{
 Time-varying DBN inference with the ARTIVA (Auto Regressive TIme VArying) model.
}
\description{
This package generates Reversible Jump MCMC (RJ-MCMC) 
sampling for approximating the posterior distribution of a time varying 
regulatory network, under the Auto Regressive TIme VArying (ARTIVA) model 
(for a detailed description of the algorithm, see Lebre et al. 
BMC Systems Biology, 2010). 

Starting from time-course gene expression measurements for a gene of
interest (referred to as "target gene") and a set of genes (referred to
as "parent genes") which may explain the expression of the target gene, the
ARTIVA procedure identifies temporal segments for which a set of
interactions occur between the "parent genes" and the "target gene".
The time points that delimit the different temporal segments are referred to as changepoints (CP).}
\details{
\tabular{ll}{
Package: \tab ARTIVA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-12-10\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
S. Lebre and G. Lelandais.

Maintainer:  S. Lebre <sophie.lebre@lsiit-cnrs.unistra.fr>. 
}
\references{
Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 2010, 4:130.
}
\keyword{DBN}
\keyword{inference}
\keyword{network inference}
\keyword{graphical model}
\keyword{time series}
