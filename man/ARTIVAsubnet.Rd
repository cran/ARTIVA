\name{ARTIVAsubnet}
\alias{ARTIVAsubnet}
\title{
Function to recover Auto Regressive TIme-VArying interactions between a gene of
interest (referred to as "target gene") and a set of genes which may explain the expression of the target gene.
}
\description{
This function generates Reversible Jump MCMC (RJ-MCMC) 
sampling for approximating the posterior distribution of a time varying 
regulatory network, under the Auto Regressive TIme VArying (ARTIVA) model 
(for a detailed description of the algorithm, see Lebre et al. 
BMC Systems Biology, 2010). 

Starting from time-course gene expression measurements for a gene of
interest (referred to as "target gene") and a set of genes (referred to
as "parent genes") which may explain the expression of the target gene, the
ARTIVA procedure identifies temporal segments for which a set of
interactions occur between the "parent genes" and the "target gene".
The time points that delimit the different temporal segments are referred to as changepoints (CP).

If the measurements time delay is short enough so that the expression of
the target gene depends more likely on the expression of the parent
genes at some previous time points, then the time delay for the
interactions can be chosen with argument \code{dyn}. In that case the set of
parent genes may contain the target gene (auto-regulation). Otherwise,
contemporaneous measurements of the parent genes are used to explain the
expression of the target gene and argument \code{dyn} is set to 0. In the
latter case, the target gene must be kept out of the set of possible
parent genes.

The ARTIVA algorithm uses a combination of efficient and robust methods:
(1) dynamical Bayesian networks (DBN) to model directed regulatory
interactions between the parent genes and the analyzed target gene
(2) RJ-MCMC sampling for inferring - simultaneously - (a) the time
position (changepoints) at which the regulatory interactions between
the parent genes and the target gene changes and (b) the regulatory
network topologies (interactions between parent and target genes)
associated with the temporal segments delimited by the changepoints.

If available, repeated measurements can be used, the design of
experiments must be specified with argument \code{dataDescription}.
}
\usage{
ARTIVAsubnet(targetData, parentData, targetName="Target",parentNames=NULL, 
dataDescription=NULL, saveEstimations=TRUE, saveIterations=FALSE, 
savePictures = TRUE, outputPath=NULL, dyn=1,  segMinLength=2, maxCP=NULL, 
maxPred=NULL, nbCPinit=NULL, CPinit=NULL, niter=50000, burn_in=NULL, 
PSRFactor=FALSE, PSRF_thres=1.05, segmentAnalysis=TRUE,  edgesThreshold=0.5,
layout = "fruchterman.reingold", cCP= 0.5, cEdges=0.5, alphaCP=1,
betaCP=0.5, alphaEdges=1, betaEdges=0.5, v0=1, gamma0=0.1, alphad2=2,
betad2=0.2, silent=FALSE)
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
  \item{targetName}{
    Name of the target gene (optional, default: targetName="Target").
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
  initial number of CPs is set to \code{min(floor(n/2), 3)}, where
  \code{n} is the number of time points in thedata).
  
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
  \code{\link{ARTIVAsubnet}} function is a RJ-MCMC algorithm which, at
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
  see package \code{igraph} for more details (default:
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
   Boolean, if \code{TRUE} messages are printed along the ARTIVA procedure (optional, default: \code{silent=FALSE}).
  }

}
\value{
\item{Samples}{
Results obtained at each iteration of the RJ-MCMC procedure. \code{Samples} is list composed of the following elements: 

1) \code{Samples$CP}: a matrix with in row the different iterations 
performed with the \code{\link{ARTIVAsubnet}} function and in column the identified positions for CPs (if the parameter dyn=1, first CP=2 and final CP=n+1, with n the number of time points);

2) \code{Samples$Edges}: a matrix with in row the different ARTIVA iterations and in column the 
number of regulatory interactions identified in each temporal phases (-1 values are default 
values, when no temporal phase exist); 

3) \code{Samples$coeff}: a matrix with in row the different ARTIVA iterations and in column the 
coefficient values corresponding to the identified regulatory interactions.

4) \code{Samples$variance}: a matrix with in row the different ARTIVA iterations and in column the 
variance values modelling data noise for each identified temporal phase. 
}
\item{Counters}{
Results obtained at each iteration of the RJ-MCMC procedure. \code{Counters} is list composed of the following elements: 

1) \code{Counters$CPMovesCount}: Number of modifications PROPOSED during ARTIVA iterations in term
of CPs (i.e. birth of a new CP, death of an existing CP, move of an existing
CP or upate of regulatory models in each temporal phases).
  
2) \code{Counters$CPMovesAcceptationPrct}: Percentage of modifications ACCEPTED during ARTIVA iterations in term
of CPs (i.e. birth of a new CP, death of an existing CP, move of an existing
CP or upate of regulatory models in each temporal phases).

3) \code{Counters$EdgesMoveCount}: Number of modifications PROPOSED during
ARTIVA iterations in term of regulatory models (i.e. birth of a new edge between parent and target genes, death of an
existing edge or update of the regression coefficient for existing edges).

4) \code{Counters$EdgesMovesAcceptationPrct}: Percentage of modifications ACCEPTED
during ARTIVA iterations in term of regulatory models (i.e. birth of a
new edge between parent and target genes, death of an existing edge or update of the regression coefficient for existing edges).

5) \code{Counters$iterations}: Total number of iterations generated by the
ARTIVA procedure.

6) \code{Counters$rcvgce}: (only if \code{PSRFactor=TRUE}) Number of iterations before
the stopping criterion is reached, i.e. before the PSRF factor is below
the threshold specified with argument \code{PSRF_thres}. The returned value
is NULL if the stopping criterion is not reached.

}
\item{CPpostDist}{
  A list of 2 tables :
  
  1) \code{CPpostDist$CPnumberPostDist}: A table containing the approximated distribution for the number of CPs. 

  2)\code{CPpostDist$CPpositionPostDist}: A table containing the approximated distribution for the position of the CPs.
}
\item{nbSegs}{
An integer equal to the number of temporal segments with the largest value observed in the posterior distribution (see previously \code{CPnumberPostDist}).  
}

\item{SegmentPostDist}{(only when parameter
  \code{segmentAnalysis=TRUE}) A list of tables:
  
  1) \code{SegmentPostDist$CPpos}: A table containing the most
  significant CP positions that delimit \code{nbSegs} temporal segments,
  according to \code{CPpostDist$CPnumber}, \code{CPpostDist$CPposition} and
  \code{segMinLength} (if parameter \code{dyn=1}, first CP is \code{2}
  and final CP is \code{n+1}, where \code{n} is the number of time
  points).

  2) \code{SegmentPostDist$edgesPostDist}: A table containing the approximate posterior distribution for the incoming edges (regulatory
interaction between parent and target genes) for each temporal segment
delimited by the CP given in \code{SegmentPostDist$CPpos} (see previously). Each raw corresponds to a segment, ordered by time.


  3) \code{SegmentPostDist$edgesCoeff} A table containing the estimated
  coefficients for the incoming edges (regulatory interaction between
  parent and target genes) for each temporal segment delimited by the CP given in \code{SegmentPostDist$CPpos} (see previously). Each raw corresponds to a segment, ordered by time.
}

\item{network}{
A table containing the information to plot (see function
\code{\link{traceNetworks}}) the network estimated with the \code{\link{ARTIVAsubnet}} procedure.
}

\item{GLOBvar}{
  A list of parameters used in the \code{\link{ARTIVAsubnet}} procedure
}

\item{HYPERvar}{
  A list of hyperparameters used in the \code{\link{ARTIVAsubnet}} procedure
}
\item{OUTvar}{
  A list of output parameters used in the \code{\link{ARTIVAsubnet}} procedure.

}
\item{targetData}{
  \code{targetData} vector given as input of the \code{\link{ARTIVAsubnet}} procedure.
}

\item{parentData}{
  \code{parentData} matrix given as input of the \code{\link{ARTIVAsubnet}} procedure.
}

}
\references{
S. Lebre, J. Becq, F. Devaux, M. P. H. Stumpf, G. Lelandais (2010)
Statistical inference of the time-varying structure of gene-regulation
networks
BMC Systems Biology, 4:130.

Gelman, A. and D. Rubin (1992)
Inference from iterative simulation using multiple sequences.
Statistical science 7 (4), 457-472.
}
\author{
S. Lebre and G. Lelandais
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ARTIVAnet}},\code{\link{ARTIVAsubnetAnalysis}},
\code{\link{choosePriors}}, \code{\link{CP.postDist}}, \code{\link{plotCP.postDist}},

\code{\link{segmentModel.postDist}}, \code{\link{traceNetworks}}
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

###
# ARTIVA analysis searching for potential interactions between the target 
# genes and a predefined list of parent genes. 
###

# Note that the number of iterations in the RJ-MCMC sampling is reduced 
# to in this example to 'niter=20000' (in order obtain a quick overview
# of the ARTIVAnet fonction, but it should be increased (e.g. up to
# 50000) for a better parameter estimation.

# Without saving the output pictures "savePictures=FALSE"
\dontrun{
ARTIVAtest = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.5, 
  niter = 20000,
  savePictures=FALSE)

# By default, the output results (pictures and estimation values) are
# saved in a folder named "ARTIVAsubnet" created in the current directory
ARTIVAtest = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.5, 
  niter = 20000)

# In order to save the results in a specific folder, for example a
# folder entitled "ARTIVA-test" in the current directory:

ARTIVAtest2 = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
  parentData = simulatedProfiles[parentGenes,],
  targetName = targetGene,
  parentNames = parentGenes,
  segMinLength = 2,
  edgesThreshold = 0.5, 
  niter = 20000,
  outputPath = "ARTIVA-test")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DBN}
\keyword{inference}
\keyword{network inference}
\keyword{graphical model}
\keyword{time series}
