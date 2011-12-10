\name{yeast}
\alias{yeast}
\docType{data}
\title{
  Yeast stress response data
}
\description{
  This dataset was produced by Lucau-Danila et al. (2005).
In this study the authors measured the changes in the mRNA concentrations for
each gene at successive times after addition of an antimitotic drug (benomyl),
in the growth media of Saccharomyces cerevisiae cells. Parallel experiments were 
conducted in different genetic contexts: the wild type strain and knock out (KO)
strains in which the genes coding for different transcription factors connected 
to drug response (YAP1, PDR1, PDR3 and YRR1) were deleted. 
For 78 genes, the measured expression values for 5 different time 
points (at 30s, 2min, 4min, 10min and 20min) are available, in each yeast strains 
(WT, DeltaYAP1, DeltaPDR1, DeltaPDR3 and DeltaYRR1). These genes are divided 
into 18 clusters of genes with condordant transcription profiles. 
In this context, regulatory associations between parent and target genes are 
proposed if the deletion of a parent gene significantly alters the expression measurements of the target genes.
}

\usage{data(yeast)}
\format{A list that comprises information for the 18 clusters of genes whose expression is identically modified in strains deleted for 
YAP1, PDR1, PDR3 and YRR1 transcription factors, compared to the wild type strain. Detailed description and more information 
concerning these clusters can be found in Lebre et al. (2010).
  \describe{
 \item{\code{Cluster 1}}{
Experimental data for genes that belong to the Cluster 1.
\code{yeast$Cluster1} is a list composed of the following elements:
1) \code{yeast$Cluster1$InitialCluster}: A matrix with in row the genes that belong to 
Cluster 1 and in column the different experimental measurements related to the 
wild type strain (WT), and strains deleted for the transcription factor YAP1 
(DeltaYAP1), PDR1 (DeltaPDR1), PDR3 (DeltaPDR3) and YRR1 (DeltaYRR1). For each
genetic context, 5 time points are available T1 = 30s, T2 = 2min, T3 = 4min,
T4 = 10min and T5 = 20min.
2) \code{yeast$Cluster1$targetData}: A vector with all gene expression measurements 
correctly formatted to be analyzed with the \code{\link{ARTIVAsubnet}} function.
3) \code{yeast$Cluster1$targetName}: The name of the analyzed cluster.
4) \code{yeast$Cluster1$GeneList}: A vector with the names of the genes that belong to 
the analyzed cluster.
5) \code{yeast$Cluster1$parentData}: A matrix with in row the four transcription factors 
for which the corresponding genes were deleted (independantly) in knock out 
strains. 1 = the gene coding for the transcription factor is present and 0 = 
the gene coding for the transcription factor is deleted. This matrix is 
correctly formatted to be analyzed with the \code{\link{ARTIVAsubnet}} function, therefore
searching for regulatory interactions between YAP1, PDR1, PDR3 and YRR1 
transcription factors and the genes that belong to the analyzed cluster.
6) \code{yeast$Cluster1$dataDescription} : A vector indicating the
ordering of the time measurements in the \code{targetData} and
\code{parentData} variables.
}
\item{\code{Cluster2 to 18}}{Experimental data for genes that belong to the Cluster 2 to 18 (see the documentation for Cluster 1).}
}
}
\source{
  Expression data were obtained from the website:

  \code{http://www.biologie.ens.fr/lgmgml/publication/benomyl/}
}
\references{
  The Early Expression of Yeast Genes Affected by Chemical Stress
A. Lucau-Danila, G. Lelandais, Z. Kozovska, V. Tanty, T. Delaveau, F. Devaux 
and C. Jacq., Mol Cell Biol., 25(5):1860-8, 2005.

Statistical inference of the time-varying structure of gene regulation networks
S. Lebre, J. Becq, F. Devaux, MP Stumpf, G. Lelandais., BMC Systems Biology 4:130, 2010.
}
\examples{
####
# Datasets related to the analysis of the genomic response of the yeast 
# Saccharomyces cerevisiae to an environmental stress induced by 
# benomyl (a toxic compound).
# Analysis of the yeast data is presented in the original article of 
# ARTIVA (Lebre et al. BMC Syst. Biol, 2010)
####

# Load the yeast dataset
data(yeast)
# This is a a list that comprises information for the 18 clusters of genes
# whose expression is identically modified in strains deleted for
# YAP1, PDR1, PDR3 and YRR1 transcription factors,
# compared to the wild type strain.

# As an illustration : analysis of one cluster
cluster=4

# Different genes in a cluster is considered as repeated measurements.
# Organisation of the different time point measurements is described in
# variable : yeast[[cluster]]$dataDescription
# Beacause of repeated measurements, the minimum segment length is set to
# segMinLength = 1.
# The parentdata is the experiment design (YAP1, PDR1, PDR3 and YRR1
# deletion) described in variable: yeast[[cluster]]$parentData
# Time delay between parent and target genes is fixed to dyn=0.
ARTIVAtest = ARTIVAsubnet(targetData = yeast[[cluster]]$targetData,
  targetName = yeast[[cluster]]$targetName,
  parentData = yeast[[cluster]]$parentData,
  parentNames = row.names(yeast[[cluster]]$parentData),
  dataDescription = yeast[[cluster]]$dataDescription,
  outputPath = paste("ARTIVA_Results_Cluster", cluster, sep = ""),
  dyn = 0,
  segMinLength = 1,
  edgesThreshold = 0.7,
  niter = 20000)
# Detailed results can be found in the folder named
# "ARTIVA_Results_Cluster4" (with the subfolders "Estimations" for
# detailled results of the estimated parameters and "Pictures" for
# graphical representations).
}
\keyword{datasets}
