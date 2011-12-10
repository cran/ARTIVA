\name{simulatedProfiles}
\alias{simulatedProfiles}
\docType{data}
\title{
 Simulated gene expression profiles dataset. 
}
\description{Simulated gene expression time series for 55 genes and 30
  timepoints generated as in Lebre et al (2010).
  In the simulation model, every gene may be target gene but only genes
  51 to 55 are parent genes.
}
\usage{data(simulatedProfiles)}
\format{
  A matrix of 55 row (genes) by 30 columns (timepoints).
}
\references{
 Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 2010, 4:130.
}
\examples{

# Load the ARTIVA R package
library(ARTIVA)

# Load the simulated profiles dataset
data(simulatedProfiles)

# Name of the target gene to be analyzed with ARTIVA 
targetGene = 1

# Names of the parent genes (typically transcription factors) 
parentGenes = c("TF1", "TF2", "TF3", "TF4", "TF5")

# Plot of the gene expression profiles for target gene and parent genes
traceGeneProfiles(targetData= simulatedProfiles[targetGene,], 
		  parentData= simulatedProfiles[parentGenes,])

###
# ARTIVA analysis searching for potential interactions between the target 
# genes and a predefined list of parent genes. 
###

# Note that the number of iterations in the RJ-MCMC sampling is reduced 
# to in this example to 'niter=20000' in order obtain a quick overview of
# the ARTIVAnet fonction, but it should be increased (e.g. up to 50000)
# for a better parameter estimation.
\dontrun{
ARTIVAtest = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  niter = 5000,
  savePictures = FALSE)
}
}
\keyword{datasets}
