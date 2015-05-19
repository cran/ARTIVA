\name{CP.postDist}
\alias{CP.postDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute the CPs posterior distribution for the ARTIVA
  network model from the the ouput samples of function ARTIVAsubnet.
}
\description{
  Using  the ouput RJ-MCMC samples of functions \code{\link{ARTIVAsubnet}}, this function estimates posterior distributions for the number of CPs and their position.
}
\usage{
CP.postDist(CPsamples, burn_in=NULL, segMinLength=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CPsamples}{A matrix with the different iterations (in row)
performed with the \code{\link{ARTIVAsubnet}} function and in column the identified positions for CPs.
  }
  \item{burn_in}{Number of initial iterations that are discarded for the
    estimation of the model distribution (posterior
    distribution). The \code{\link{ARTIVAsubnet}} function is a RJ-MCMC
    algorithm which, at each iteration, randomly samples a new
    configuration of the time-varying regulatory network from
    probability distributions based on constructing a Markov chain that
    has the network model distribution as its equilibrium distribution
    (The equilibrium distribution is obtained when the Markov Chain
    converges, which requires a large number of iterations).
    Typically, initial iterations are notconfident because the Markov
    Chain has not stabilized. The burn-in samples allow to not consider
    these initial iterations in the final analysis (optional, default:
    \code{burn_in=NULL}, if \code{burn_in=NULL} then the first 25\% of
    the iterations is left for \code{burn_in}).
  }
  \item{segMinLength}{Minimal length (number of time points) to define a
    temporal segment. Must be - strictly -  greater than 1 if there is
    no repeated measurements for each time point in arguments
    \code{targetData} and \code{parentData} (optional, default: \code{segMinLength=2}).
  }
}

\value{
  A list of 4 elements:
  
  1) \code{CPnumber}: a table containing the approximate posterior distribution for the number of CPs.

  2) \code{CPposition}: a table containing the approximate posterior
  distribution for the CPs position.

  3) \code{estimatedCPnumber}: number of CP position with the greatest
  posterior probability according to the approximate posterior
  distribution for the number of CPs \code{CPnumber}.

  4) \code{estimatedCPpos}: a table containing the
  \code{estimatedCPnumber} most significant CP positions according to
  \code{CPnumber}, \code{CPposition} and \code{segMinLength} (if parameter \code{dyn=1}, first CP is \code{2} and final CP is \code{n+1}, where \code{n} is the number of time points).

}
\references{
Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology, 4:130, 2010.
}
\author{
S. Lebre and G. Lelandais.
}

\seealso{
\code{\link{ARTIVAsubnet}}, \code{\link{ARTIVAnet}},
\code{\link{plotCP.postDist}}, \code{\link{ARTIVAsubnetAnalysis}}
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
outCPpostDist = CP.postDist(ARTIVAtest$Samples$CP, burn_in=10000, 
			    segMinLength=3)

# plot the CP posterior distribution
plotCP.postDist(outCPpostDist, targetName=paste("Target", targetGene), 
		  estimatedCPpos=outCPpostDist$estimatedCPpos)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{util}