\name{ARTIVAsubnetAnalysis}
\alias{ARTIVAsubnetAnalysis}
\title{
  Function to estimate a regulatory time-varying network from the
  output of function ARTIVAsubnet.
}
\description{
  This function estimates a regulatory time-varying network from the
  output of function \code{\link{ARTIVAsubnet}}.
  A graphical representation in a pdf file and estimated values are
  provided in text files.
  This function is used in function \code{\link{ARTIVAsubnet}}
  when parameter \code{segmentAnalysis=TRUE}. This function can be used
  separately for re-computing a time-varying network from the output of function
  \code{\link{ARTIVAsubnet}}  with new analysis parameters
  \code{segMinLength},  \code{edgesThreshold}, \code{CPpos},
  \code{layout}, ... see detail below.
}

\usage{
ARTIVAsubnetAnalysis(ARTIVAsubnet=NULL, CPpostDist=NULL, CPsamples=NULL,
coefSamples=NULL, TFnumber=NULL, segMinLength=2, edgesThreshold=0.5,
burn_in=NULL, CPpos=NULL,targetData=NULL, parentData=NULL,
targetName=NULL,parentNames=NULL, savePictures=TRUE,saveEstimations=TRUE,
outputPath=NULL,layout="fruchterman.reingold", silent=FALSE,
inARTIVAsubnet=FALSE , onepage= FALSE)
}

\arguments{
  \item{ARTIVAsubnet}{
Ouput of function \code{\link{ARTIVAsubnet}}, a list containing \code{Samples},
  \code{Counters}, \code{CPpostDist}, \code{nbSegs},
  \code{SegmentPostDist}, \code{network}, ...
  (optional, default: \code{ARTIVAsubnet=NULL}, if
  \code{ARTIVAsubnet=NULL} then parameters \code{CPpostDist},
  \code{CPsamples}, \code{coefSamples}, \code{TFnumber} must be not null.
  }


  \item{CPpostDist}{ A list of 2 tables :
  1)\code{CPpostDist$CPnumberPostDist}: A table containing the distribution for
  the number of CPs approximated with \code{\link{ARTIVAsubnet}}. 
  2)\code{CPpostDist$CPpositionPostDist}: A table containing the distribution for the position of the CPs approximated with \code{\link{ARTIVAsubnet}}.
    (optional, default: \code{CPpostDist=NULL}, but
  \code{CPpostDist} must be given when parameter \code{ARTIVAsubnet=NULL})}

  \item{CPsamples}{A matrix with the different iterations (in row)
performed with the \code{\link{ARTIVAsubnet}} function and in column the identified positions for CPs.  (optional, default: \code{CPsamples=NULL}, but
  \code{CPsamples} must be given when parameter \code{ARTIVAsubnet=NULL})
  }

  \item{coefSamples}{A matrix with the different (in row)performed with
    the \code{\link{ARTIVAsubnet}} function and in column the 
coefficient values corresponding to the identified regulatory interactions.
  (optional, default: \code{coefSamples=NULL}, but
  \code{coefSamples} must be given when parameter \code{ARTIVAsubnet=NULL}) }

  \item{TFnumber}{Number of parent genes in the data
  \code{parentData} used in the \code{\link{ARTIVAsubnet}} function.
   (optional, default: \code{TFnumber=NULL}, but
  \code{TFnumber} must be given when parameter \code{ARTIVAsubnet=NULL})}
  
  \item{segMinLength}{Minimal length (number of time points) to define a
    temporal segment. Must be - strictly -  greater than 1 if there is
    no repeated measurements for each time point in arguments
    \code{targetData} and \code{parentData} (optional, default: \code{segMinLength=2}).}
  
  \item{edgesThreshold}{Probability threshold for the selection of the
    edges of the time-varying regulatory network when
    \code{segmentAnalysis=TRUE} (optional, default:  \code{edgesThreshold=0.5}).
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
    
 \item{CPpos}{A table containing the  desired most significant CP positions (optional, default:
    \code{CPpos=NULL}, if \code{CPpos=NULL} then \code{CPpos} is
  evaluated as in \code{\link{CP.postDist}} function. )
  }

  \item{targetData}{
      A vector with the temporal gene expression measurements for the
      target gene (i.e. the gene whose regulation factors are looked
      for). (optional, default: \code{targetData=NULL}, if not null then
      the target data is plotted).
  }
  \item{parentData}{  A matrix (or a vector if only 1 parent gene) with the temporal gene expression measurements for the proposed parent genes (i.e. potential
    regulation factors). Parent genes are shown in row and expression values
    in column. (optional, default: \code{parentData=NULL}, if not null then
      the parent data is plotted).
  }

  
  \item{targetName}{
     Name of the target gene (optional, default: targetName="Target").
  }

  \item{parentNames}{
     A vector with the names for parent gene(s) (optional, default: \code{parentNames=NULL}).
  }

  \item{savePictures}{
    Boolean, if \code{TRUE} all estimated posterior distributions and
    networks are plotted in a pdf file either in a new sub folder named
    "ARTIVA_Results" created by default in the current folder or in a
    folder specified with argument \code{outputPath} (see below)
    (optional, default: \code{savePictures=TRUE}).}

  \item{saveEstimations}{
   Boolean, if \code{TRUE} all estimated posterior distributions are saved
    as text files either in a new sub folder named "ARTIVA_Results" created
    by default in the current folder or in a folder specified with argument
    \code{outputPath} (see below) (optional, default:
    \code{saveEstimations=TRUE}).}

  \item{outputPath}{
   File path to a folder in which the output results have to be saved,
    either a complete path or the name of a folder to be created in the
    current directory (optional, default: \code{outputPath=NULL}).}

  \item{layout}{
  Name of the function determining the placement of the vertices for
  drawing a graph, possible values among others: \code{"random",
    "circle", "sphere", "fruchterman.reingold", "kamada.kawai",
    "spring",}

  \code{"reingold.tilford", "fruchterman.reingold.grid"},
  see package \code{igraph} for more details (default:
  \code{layout="fruchterman.reingold"}).}
 
\item{silent}{
  Boolean, if \code{TRUE} messages are printed along the ARTIVA procedure (optional, default: \code{silent=FALSE}).
  }

  \item{inARTIVAsubnet}{
  Boolean, if \code{TRUE}, general information already printed in
  function  \code{\link{ARTIVAsubnet}} are not printed a second time (optional, default: \code{inARTIVAsubnet=FALSE}).
  }
  
  \item{onepage}{
 Boolean, if \code{TRUE}, all output pictures are plotted on one page
 only (optional, default: \code{onepage=FALSE}.
  }
 
 
   
}
\value{

  \item{nbSegs}{An integer equal to the number of temporal segments with
  the largest value observed in the CP number posterior distribution
  (from \code{CPpostDist$CPnumberPostDist}). 
}

  \item{CPposition}{A table containing the most significant CP positions
    that delimit \code{nbSegs} temporal segments, according to
  \code{CPnumber}, \code{CPposition} and \code{segMinLength} (if
  parameter \code{dyn=1}, first CP is \code{2} and final CP is
  \code{n+1}, where \code{n} is the number of time points).
}

  \item{SegmentPostDist}{Output of function
    \code{\link{segmentModel.postDist}}. A list of tables:
    
  1) \code{SegmentPostDist$CPpos}: A table containing the most
  significant CP positions that delimit \code{nbSegs} temporal segments,
  according to \code{CPpostDist$CPnumber}, \code{CPpostDist$CPposition} and
  \code{segMinLength} (if parameter \code{dyn=1}, first CP is \code{2}
  and final CP is \code{n+1}, where \code{n} is the number of time points).

  2) \code{SegmentPostDist$edgesPostDist}: A table containing the approximate posterior distribution for the incoming edges (regulatory
interaction between parent and target genes) for each temporal segment
delimited by the CP given in \code{SegmentPostDist$CPpos} (see previously). Each raw corresponds to a segment, ordered by time.

  3) \code{SegmentPostDist$edgesCoeff} A table containing the estimated coefficients for the incoming edges (regulatory interaction between parent and target
  genes) for each temporal segment delimited by the CP given in \code{SegmentPostDist$CPpos} (see previously). Each raw corresponds to a segment, ordered by time.
}
   \item{network}{
    A table containing the information to plot (see function
\code{\link{traceNetworks}}) the network estimated with the \code{\link{ARTIVAsubnet}} procedure.
  }
}

\references{
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais (2010)
Statistical inference of the time-varying structure of gene-regulation
networks
BMC Systems Biology, 4:130.
}
\author{
S. Lebre and G. Lelandais
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ARTIVAsubnet}}, \code{\link{ARTIVAnet}},
  \code{\link{traceNetworks}}, \code{\link{traceGeneProfiles}},
  \code{\link{CP.postDist}}, \code{\link{plotCP.postDist}}.
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


# Note that the number of iterations in the RJ-MCMC sampling is reduced 
# to 'niter=20000' in this example, but it should be increased (e.g. up to
# 50000) for a better estimation.

# Run the ARTIVAsubnet function
\dontrun{
ARTIVAtest = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.6, 
  niter= 20000,
  savePictures=FALSE)

# Re-compute a time-varying network from the output of function 
# ARTIVAsubnet with new analysis parameters
analysis2 = ARTIVAsubnetAnalysis(ARTIVAsubnet=ARTIVAtest,
  segMinLength = 3,
  edgesThreshold = 0.5,
  outputPath="ARTIVAsubnet2",
  savePictures=FALSE)

# Trace the obtained network.
traceNetworks(analysis2$network, edgesThreshold = 0.3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DBN}
\keyword{inference}
\keyword{network inference}
\keyword{graphical model}
\keyword{time series}
