\name{ARTIVAnet}
\alias{ARTIVAnet}
\title{
Function to run the ARTIVA procedure for Auto Regressive TIme-VArying
network inference on several target genes
}
\description{
This function runs function \code{\link{ARTIVAsubnet}} for all target
genes in \code{targetData} successively. This function generates
Reversible Jump MCMC (RJ-MCMC) sampling for approximating the posterior
distribution of a time varying regulatory network, under the Auto
Regressive TIme VArying (ARTIVA) model (for a detailed description of
the algorithm, see \code{\link{ARTIVAsubnet}} and see Lebre et al. 
BMC Systems Biology, 2010). A network representing the interactions
between the factor genes and the target genes is estimated and plotted.
}
\usage{
ARTIVAnet(targetData, parentData, targetNames = NULL, parentNames =
NULL, dataDescription=NULL, saveEstimations=TRUE, saveIterations=FALSE,
savePictures = TRUE, outputPath=NULL, dyn=1, segMinLength=2, maxCP=NULL,
maxPred=NULL, nbCPinit=NULL, CPinit=NULL, niter=50000, burn_in=NULL,
PSRFactor=FALSE,PSRF_thres=1.05, segmentAnalysis=TRUE,
edgesThreshold=0.5, layout = "fruchterman.reingold", cCP= 0.5,
cEdges=0.5, alphaCP=1, betaCP=0.5, alphaEdges=1, betaEdges=0.5, v0=1,
gamma0=0.1, alphad2=2, betad2=0.2, silent=FALSE)
}
\arguments{
  \item{targetData}{
    A vector with the temporal gene expression measurements for the target gene (i.e. the gene whose regulation factors are looked for).
  }
  \item{parentData}{
    A matrix (or a vector if only 1 parent gene) with the temporal gene
    expression measurements for the proposed parent genes (i.e. potential
    regulation factors). Parent genes are shown in row and expression values
    in column. For computational reasons, we advise not to test
    simultaneously more than 10 parent genes.
  }
  \item{targetNames}{
      A vector with the names for target gene(s) (optional, default:
    \code{targetNames=NULL}
   }
  \item{parentNames}{
    A vector with the names for parent gene(s) (optional, default:
    \code{parentNames=NULL}).
  }
  \item{dataDescription}{
    (Required only when the gene expression measurements contain repeated
    values for the same time points). A vector indicating the ordering of
    the time measurements in the data. For example
    dataDescription=rep(1:n, each=2), if there are two measurements for
    each time point AND the repetitions for each time point are next to
    each other. Note that temporal gene expression measurements have to be
    organized identically in arguments \code{targetData} and
    \code{parentData} (optional, default: \code{dataDescription=NULL}).
  }
  \item{saveEstimations}{
    Boolean, if \code{TRUE} all estimated posterior distributions are saved
    as text files either in a new sub folder named "ARTIVA_Results" created
    by default in the current folder or in a folder specified with argument
    \code{outputPath} (see below) (optional, default:
    \code{saveEstimations=TRUE}).
  }
\item{saveIterations}{
    Boolean, if \code{TRUE} the configuration for all iterations is saved as
    text files either in a new sub folder named "ARTIVA_Results" created
    by default in the current folder or in a folder specified with argument
    \code{outputPath} (see below) (optional, default:
    \code{saveIterations=FALSE}).
  }
  \item{savePictures}{
    Boolean, if \code{TRUE} all estimated posterior distributions and
    networks are plotted in a pdf file either in a new sub folder named
    "ARTIVA_Results" created by default in the current folder or in a
    folder specified with argument \code{outputPath} (see below)
    (optional, default: \code{savePictures=TRUE}).
  }
  \item{outputPath}{
    File path to a folder in which the output results have to be saved,
    either a complete path or the name of a folder to be created in the
    current directory (optional, default: \code{outputPath=NULL}).
  }
  \item{dyn}{
    Time delay to be considered in the auto-regressive process, between the
    temporal expression measurements of the analyzed target gene and the
    ones of the parent genes (optional, default: \code{dyn=1}).
  }
  \item{segMinLength}{
    Minimal length (number of time points) to define a temporal
    segment. Must be - strictly -  greater than 1 if there is no
    repeated measurements for each time point in arguments
    \code{targetData} and \code{parentData} (optional, default: \code{segMinLength=2}).
}
  \item{maxCP}{
    Maximal number of CPs to be considered by the ARTIVA inference
    procedure. Note that for a temporal course with \code{n} time points
  the maximal number of CPs is \code{(n-1-dyn)} (see before for a
  description of argument \code{dyn}). For long temporal courses (more
  than 20 time points), we advise - for computational reasons - to limit
  the maximal number of CP to 15 (optional, default: \code{maxCP=NULL},
  if \code{maxCP=NULL} then the maximal number of CPs is set to
  \code{maxCP=min(round((n-1-dyn)/segMinLength)-1, 15)} where \code{n} is the
  number of time points).
}
  \item{maxPred}{
    Maximal number of simultaneous incoming edges for each segment of the
    regulatory network estimated for the target gene (default:
    \code{maxPred=NULL}, if \code{maxPred=NULL} then the maximal number of
    incoming edges is set to
    
    \code{maxPred=min(dim(parentData)[1],15)})
}
  \item{nbCPinit}{
Number of CPs to be considered at the algorithm initialization
  (optional, default: \code{nbCPinit=NULL}, if \code{nbCPinit=NULL} then the
  initial number of CPs is randomly set to a value between 0 and \code{maxCP/2}.).  
}
  \item{CPinit}{
A vector with the initial positions for potential CPs.
(optional, default: \code{CPinit = NULL}, if \code{CPinit = NULL}
then the initial vector is chosen randomly according to priors, with
  \code{nbCPinit} changepoints (see previous argument \code{nbCPinit}) ). 
}
  \item{niter}{
Number of iterations to be performed in the RJ-MCMC sampling 
(optional, default: \code{niter = 50000}).
}
  \item{burn_in}{
Number of initial iterations that are discarded for the estimation of
the model distribution (posterior distribution). The
  \code{\link{ARTIVAnet}} function is a RJ-MCMC algorithm which, at
  each iteration, randomly samples a new configuration of the
  time-varying regulatory network from probability distributions based
  on constructing a Markov chain that has the network model distribution
  as its equilibrium distribution (The equilibrium distribution is
  obtained when the Markov Chain converges, which requires a large
  number of iterations).
  Typically, initial iterations are notconfident because the Markov Chain
  has not stabilized.  The burn-in samples allow to not consider these
  initial iterations in the final analysis (optional, default:
  \code{burn_in=NULL}, if \code{burn_in=NULL} then the first 25\% of the
  iterations is left for \code{burn_in}).

}
  \item{PSRFactor}{
    Boolean, if \code{TRUE} the RJ-MCMC procedure is stopped when the Potential 
Scale Ratio Factor or PSRF (which is a usual convergence criterion) is
  below a specified threshold (see   documentation for argument
  \code{PSRF_thres} below) or when the maximal number of iterations is
  reached (see documentation for argument   \code{niter} previously) (optional, default: \code{PSRFactor=FALSE}). 
For more details about the PSRF criterion, see Gelman and Rubin (1992).
}  
  \item{PSRF_thres}{
(Only when \code{PSRFactor=TRUE}) RJ-MCMC stopping threshold: the
  RJ-MCMC procedure is stopped when the Potential Scale Ratio Factor
  (PSRF) is below this threshold (see documentation for argument
  \code{PSRFactor} previously). Values of the PSRF below 1.1 are usually
  taken as indication of sufficient convergence (optional, default:
  \code{PSRF_threshold=1.05}).  
} 
  \item{segmentAnalysis}{
Boolean, if \code{TRUE} the posterior distribution for the edges in each
  temporal segment delimited by the estimated changepoints (CPs) is
  computed. (The estimated CPs are the \code{k} CPs with maximal posterior
  probability, where \code{k} is the number of CPs having the maximal
  porsterior probability (such that each temporal segment is larger than
  \code{segMinLength}) (see  documentation below for output values
  \code{CPpostDist}, \code{nbSegs}, \code{CPposition}. If \code{segmentAnalysis=FALSE} only the posterior  distributions of the CPs number and CPs positions are computed  (optional, default: \code{segmentAnalysis=TRUE}).
}
  \item{edgesThreshold}{
Probability threshold for the selection of the edges of the time-varying
regulatory network when \code{segmentAnalysis=TRUE}
(optional, default:  \code{edgesThreshold=0.5}).
  }
  \item{layout}{
Name of the function determining the placement of the vertices for
  drawing a graph, possible values among others: \code{"random",
    "circle", "sphere", "fruchterman.reingold", "kamada.kawai",
    "spring",}

  \code{"reingold.tilford", "fruchterman.reingold.grid"},
  see package \code{igraph0} for more details (default:
  \code{layout="fruchterman.reingold"}).
  }
  
  \item{cCP}{
Maximal probability to propose the birth (resp. death) of a changepoint
(CP) during the RJ-MCMC iterations (optional, default: \code{cCP=0.5}).
}
  \item{cEdges}{
Maximal probability - when a move update of the network topology is
chosen - to propose ecah edge move (birth or death of an edge) within a
temporal segment (optional, default: \code{cEdges=0.5}).
}
  \item{alphaCP}{
Hyperparameter for sampling the number \code{k} of CPs. \code{k} follows
a Gamma distribution: Gamma(alphaCP, betaCP (see below)) (optional,
default: \code{alphaCP=1}). Note that function
\code{\link{choosePriors}} can be used to choose  \code{alphaCP} (or \code{betaCP}) according to the desired dimension penalisation.
}
  \item{betaCP}{
Hyperparameter for sampling the number \code{k} of CPs. \code{k} follows
a Gamma distribution: Gamma(alphaCP (see before), \code{betaCP})
(optional, default:  \code{betaCP=0.5}). Note that function
\code{\link{choosePriors}} can be used to choose \code{betaCP} (or
\code{alphaCP}) according to the desired dimension penalisation.
}
  \item{alphaEdges}{
Hyperparameter for sampling the number \code{l} of regulatory
interactions between the target gene and the parent genes. \code{l}
follows a Gamma distribution: \code{rgamma(1, shape=alphaEdges,
rate=betaEdges)} (see \code{betaEdges} below), (default:
\code{alphaEdges=1}). Function \code{\link{choosePriors}} can be used to
choose \code{alphaEdges} (or \code{betaEdges}) according to
the desired dimension penalisation.
}
  \item{betaEdges}{
Hyperparameter for sampling the number \code{l} of regulatory
interactions between the target gene and the parent genes. \code{l}
follows a \code{Gamma} distribution: \code{rgamma(1, shape=alphaEdges,
rate=betaEdges)} (see \code{alphaEdges} upper) (optional, default:
\code{betaEdges=0.5}). Function \code{\link{choosePriors}} can be used
to choose \code{betaEdges} (or \code{alphaEdges}) according to the desired
dimension penalisation.
}
  \item{v0}{
Hyperparameter for sampling the noise variance (denoted by \code{Sig2})
in the auto-regressive model defining the regulatory time vrarying network. The prior
distribution of the noise variance is an Inverse Gamma distribution with
shape parameter \code{v0/2} and scale parameter \code{gamma0/2}: \code{rinvgamma(1,shape=v0/2, scale = gamma0)}, (optional, default: \code{v0=1}).

}
  \item{gamma0}{
Hyperparameter for sampling the noise variance (denoted by \code{Sig2}) in the ARTIVA package) in the auto-regressive
model defining the regulatory time vrarying network. The prior
distribution of the noise variance is an Inverse Gamma distribution with
shape parameter \code{v0/2} and scale parameter \code{gamma0/2}: \code{rinvgamma(1, shape=v0/2, scale = gamma0)}, (optional, default: \code{gamma0=0.1}).
}
  \item{alphad2}{
Hyperparameter for sampling a parameter  that represents the expected
  signal-to-noise ratio (denoted by \code{delta2} in the ARTIVA
  package). It is sampled according to an Inverse Gamma distribution:
  \code{rinvgamma(1, shape=alphad2, scale=betad2)}, (optional, default:
  \code{alphad2=2}).
}
  \item{betad2}{
Hyperparameter for sampling a parameter  that represents the expected
  signal-to-noise ratio (denoted by \code{delta2} in the ARTIVA
  package). It is sampled according to an Inverse Gamma distribution:
  \code{rinvgamma(1, shape=alphad2, scale=betad2)}, (optional, default:
  \code{betad2=0.2}).
}
  \item{silent}{
   Boolean, if \code{TRUE} messages are printed along the ARTIVA procedure (optional, default: \code{silent=FALSE})
  }
}
\value{
A table containing the information to plot (see function
\code{traceNetworks}) the global network estimated by the
  \code{ARTIVAnet} procedure.
All results are plotted in a pdf file (when choosing \code{savePictures =
TRUE}) in folder \code{outputPath}.
All numerical results (see \code{\link{ARTIVAsubnet}} output values
  documentation)) are saved in text files (when choosing \code{saveEstimations=TRUE} and/or \code{saveIterations=TRUE}) in folder \code{outputPath}.
}
\references{
Statistical inference of the time-varying structure of gene-regulation networks
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais, BMC Systems Biology 2010, 4:130.

Inference from iterative simulation using multiple sequences.
Gelman, A. and D. Rubin, Statistical science 7 (4), 457-472, 1992.
}
\author{
S. Lebre and G. Lelandais
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ARTIVAsubnet}}, \code{\link{choosePriors}},
\code{\link{ARTIVAsubnetAnalysis}},\code{\link{CP.postDist}},

\code{\link{segmentModel.postDist}}
}
\examples{
library(ARTIVA)

# Load the dataset with simulated gene expression profiles
data(simulatedProfiles)

# List of target genes to be analyzed independantly with ARTIVA 
targetGenes = c(1, 10, 20, "TF3", 45, 50)

# Names of the parent genes (typically transcription factors) 
parentGenes = c("TF1", "TF2", "TF3", "TF4", "TF5")

###
# ARTIVA analysis searching for potential interactions between each target 
# genes and a predefined list of parent genes. 
###

# Note that the number of iterations in the RJ-MCMC sampling is reduced 
# to 'niter=20000' in this example, but it should be increased (e.g. up to
# 50000) for a better estimation.

\dontrun{
ARTIVAtest1=ARTIVAnet(targetData = simulatedProfiles[targetGenes,],
  parentData = simulatedProfiles[parentGenes,], 
  targetNames= targetGenes,
  parentNames = parentGenes,
  niter = 20000,
  savePictures=FALSE)

# Note that function ARTIVAnet calls fonction ARTIVAsubnet for each 
# target gene successively and provides a global estimated regulatory 
# network entitled "ARTIVA_FinalNetwork.pdf" in addition to the output 
# results given by function ARTIVAsubnet.

## Gene names for the target and the parent genes, minimum segment length, 
## threshold for the edges selection can be specified as follow:
ARTIVAtest2=ARTIVAnet(targetData = simulatedProfiles[targetGenes,],
  parentData = simulatedProfiles[parentGenes,], 
  targetNames= targetGenes,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.5, 
  niter = 20000,
  outputPath = "ARTIVA-test")
}

# By default, the output results (pictures and estimation values) are
# saved in a folder named "ARTIVAnet" created in the current directory
# In order to save the results in a specific folder, for example a
# folder entitled "ARTIVA-test" in the current directory:

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DBN}
\keyword{inference}
\keyword{network inference}
\keyword{graphical model}
\keyword{time series}
