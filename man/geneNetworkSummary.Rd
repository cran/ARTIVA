\name{geneNetworkSummary}
\alias{geneNetworkSummary}
\title{
Function to 
}
\description{
  This function is used for printing a summary of the gene network
  estimated with the ARTIVA procedure (\code{\link{ARTIVAnet}}, \code{\link{ARTIVAsubnet}}) for Auto Regressive TIme-VArying network inference.
}
\usage{
geneNetworkSummary(ARTIVAnet, edgesThreshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ARTIVAnet}{
Table containing the information to plot a time-varying regulatory
network. In particular, this table can be obtained with function
\code{\link{ARTIVAsubnet}},

\code{\link{ARTIVAsubnetAnalysis}} (output value \code{network}) or
\code{\link{ARTIVAnet}} (unique output value). Each row of the table
describes one edge. The columns, entitled Target, CPini, CPfinal,
Parent, PostProb, describe the name of the target gene, the changepoints
defining the start and the end of the regulation, the parent name and
the estimated posterior probability of the edge.
  }
  \item{edgesThreshold}{Probability threshold for the selection of the
    edges to be plotted.
  }

}

\value{
NULL
}

\references{
Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 4:130, 2010.
}

\author{
Original version by S. Lebre and G. Lelandais, contribution of
D. Servillo to the final version.
}

\seealso{
  \code{\link{ARTIVAnet}},\code{\link{ARTIVAsubnet}},
  \code{\link{ARTIVAsubnetAnalysis}}, \code{\link{CP.postDist}},
  
  \code{\link{segmentModel.postDist}}, \code{\link{plotCP.postDist}}
}

\examples{
# Load the ARTIVA R package
library(ARTIVA)

# Load the dataset with simulated gene expression profiles
data(simulatedProfiles)

# Name of the target gene to be analyzed with ARTIVA 
targetGene = 1

# Names of the parent genes (typically transcription factors) 
parentGenes = c("TF1", "TF2", "TF3", "TF4", "TF5")


# Run the ARTIVAsubnet function
# Note that the number of iterations in the RJ-MCMC sampling is reduced 
# to 'niter=20000' in this example, but it should be increased (e.g. up to
# 50000) for a better estimation.

\dontrun{
ARTIVAtest = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.6, 
  niter= 2000,
  savePictures=FALSE)

# Print a summary of the obtained network 
geneNetworkSummary(ARTIVAtest$network, edgesThreshold = 0.3)

# List of target genes to be analyzed independantly with ARTIVA
targetGenes = c("TF3", 45, 50)
ARTIVAtest2 = ARTIVAnet(targetData = simulatedProfiles[targetGenes,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGenes,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.6, 
  niter= 2000,
  savePictures=FALSE)

# Print a summary of the obtained network 
geneNetworkSummary(ARTIVAtest2, edgesThreshold = 0.3)

# Re-compute a time-varying network from the output of function 
# ARTIVAsubnet with new analysis parameters
analysis2 = ARTIVAsubnetAnalysis(ARTIVAsubnet=ARTIVAtest,
  segMinLength = 3,
  edgesThreshold = 0.5,
  outputPath="ARTIVAsubnet2",
  savePictures=FALSE)

# Print a summary of the network obtained with the 2nd analysis.
geneNetworkSummary(analysis2$network, edgesThreshold = 0.3)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}
\keyword{util}