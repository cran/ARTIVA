\name{drosophila}
\alias{drosophila}
\docType{data}
\title{
Drosphila life cycyle time series by Arbeitman et al 2002.
}
\description{
Gene expression time series for 4028 genes involved in the life cycle
of Drosophila melanogaster. The microarray data measured
gene expression levels during all four major stages of morphogenesis:
embryo, larva, pupa and adult (67 time points).
Data published in Arbeitman et al. (2002) and used in Lebre et al. (2010)
}
\usage{data(drosophila)}
\format{
  A matrix of 4028 row (genes) by 67 columns (timepoints) containing the
  gene expression time series.
}
\source{
The data has been published in Arbeitman et al (2002).
}
\references{
 Gene expression during the life cycle of Drosophila melanogaster.
 M.N Arbeitman, E.E.M. Furlong, F. Imam, E. Johson, B.H. Null, B.S. Baker, M.A. Krasnow, M.P. Scott, R.W. Davis and K.P. White.  Science, 297(5590):2270-2275, 2002.

Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 2010, 4:130.
}
\examples{
#load the data
data(drosophila)

#plot time serie for a chosen gene 
gene=1
plot(1:67,drosophila[gene,],type="l",xlab="Timepoints", ylab="LogRatio",
main=row.names(drosophila)[gene])
}
\keyword{datasets}
