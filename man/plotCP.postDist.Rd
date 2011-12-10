\name{plotCP.postDist}
\alias{plotCP.postDist}
\title{
Function to plot the estimated posterior distribution for the changepoints (CPs)
  number and position
}
\description{
This function is used for plotting the estimated changepoint number and
position posterior distribution after running the ARTIVA procedure
(function \code{\link{ARTIVAsubnet}}) for Auto Regressive TIme-VArying network inference.
}
\usage{
plotCP.postDist(CPpostDist, targetName = NULL, onepage = TRUE,
color1 = "green", color2 = "black", estimatedCPpos=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CPpostDist}{
 A list of 2 tables :
  1)\code{CPpostDist$CPnumberPostDist}: A table containing the distribution for
  the number of CPs approximated with \code{\link{ARTIVAsubnet}}. 
  2)\code{CPpostDist$CPpositionPostDist}: A table containing the
  distribution for the position of the CPs approximated with function
  \code{\link{ARTIVAsubnet}} or \code{\link{CP.postDist}} 
}
  \item{targetName}{
Name of the target gene (optional, default: \code{targetName=NULL}).
}
  \item{onepage}{
Boolean, if TRUE the two estimated posterior distributions are
plotted in one window next to each other (optional, default: \code{mfrow=TRUE}).
}
  \item{color1}{
Color for plotting the estimated posterior distribution for the
changepoints (CPs) number (default: \code{color1="green"}).
  }
  \item{color2}{
Color for plotting the estimated posterior distribution for the
changepoints (CPs) position (default: \code{color2="black"}).
}
\item{estimatedCPpos}{CP positions to be highlighted as most
  significant, e.g. CP positions estimated with function
  \code{\link{CP.postDist}}  (optional, default:
  \code{estimatedCPpos=NULL}, if  \code{estimatedCPpos=NULL} then the
  number of highlighted CPs is the maximum of \code{CPpostDist$CPnumberPostDist} and the positions are the top best of \code{CPpostDist$CPpositionPostDist}).
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
S. Lebre and G. Lelandais.
}

\seealso{
\code{\link{ARTIVAnet}}, \code{\link{ARTIVAsubnet}},
\code{\link{CP.postDist}}, \code{\link{segmentModel.postDist}},
\code{\link{ARTIVAsubnetAnalysis}}
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

# run ARTIVAsubnet
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
  niter= 20000,
  savePictures=FALSE)

# compute the PC posterior distribution with other parameters
outCPpostDist = CP.postDist(ARTIVAtest$Samples$CP, burn_in=500, 
			    segMinLength=3)

# plot the CP posterior distribution
plotCP.postDist(outCPpostDist, targetName=paste("Target", targetGene), 
		  estimatedCPpos=outCPpostDist$estimatedCPpos)

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}
\keyword{util}