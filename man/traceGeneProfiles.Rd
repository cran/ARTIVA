\name{traceGeneProfiles}
\alias{traceGeneProfiles}
\title{
Function to plot the gene expression profiles
}
\description{
  This function is used for plotting the target and parent gene
  expression profiles.
}
\usage{
traceGeneProfiles(targetData, parentData,  dataDescription=NULL,
targetColor = "#046380", parentColor = "#8E3557", onepage=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targetData}{
A matrix (or a vector if only one target gene) with the temporal gene expression measurements for
the target genes (i.e. the genes whose regulation factors are looked for).
Target genes are shown in row and expression values in column. 
}
  \item{parentData}{
A matrix (or a vector if only 1 parent gene) with the temporal gene
expression measurements for the proposed parent genes (i.e. potential
regulation factors). Parent genes are shown in row and expression values in column.
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

  \item{targetColor}{
    Color for plotting the target genes expression profiles (optional, default \code{targetColor= "#046380"})
  }

  \item{parentColor}{
 Color for plotting the parent genes expression profiles (optional,
 default \code{parentColor= "#8E3557"})
  }

  \item{onepage}{
 Boolean, if \code{TRUE}, all output pictures are plotted on one page
 only (optional, default: \code{onepage=TRUE}.
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
\code{\link{ARTIVAsubnet}}, \code{\link{ARTIVAnet}} 
}

\examples{
# Load the R package ARTIVA
library(ARTIVA)

# Load the dataset with simulated gene expression profiles
data(simulatedProfiles)

# Name of the target gene to be analyzed with ARTIVA 
targetGene = 1

# Names of the parent genes (typically transcription factors) 
parentGenes = c("TF1", "TF2", "TF3", "TF4", "TF5")

# Plot of the gene expression profiles for target gene and parent genes
traceGeneProfiles(targetData= simulatedProfiles[targetGene,], 
		  parentData= simulatedProfiles[parentGenes,])

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}
\keyword{util}