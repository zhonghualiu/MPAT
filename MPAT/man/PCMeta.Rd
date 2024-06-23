\name{PCMeta}
\alias{PCMeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use PCMeta method to combine Z testing statistics.
}
\description{
Use PCMeta method to combine Z testing statistics.
}
\usage{
PCMeta(Z.vec, Sigma,SigmaMeta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Z.vec}{ Z.vec is a vector representing the Z testing statistics for multiple phenotypes.
}
  \item{Sigma}{
  Sigma is the correlation matrix among the Z testing statistics. }
    \item{SigmaMeta}{
 The correlation among X which can be estimated by simulation using SigmaMetaEstimate.}
}
\value{
\item{p-value}{An overall p-value assessing the association between a genetic variant and multiple phenotypes is returned.}
}
\references{
MPAT: an R package for multivariate phenotype association studies.
}
\author{
Zhonghua Liu  <zliu@mail.harvard.edu>
}

\examples{
data(lipids)
lipids = as.data.frame(lipids)
lipids_zscore = as.matrix(lipids[,c("Zscore.HDL","Zscore.LDL","Zscore.TG","Zscore.TC")])
Sigma = cor(lipids_zscore )
SigmaMeta = SigmaMetaEstimate(Sigma,simNum=1000) 
PCMeta(Z.vec=lipids_zscore[1,],Sigma=Sigma,SigmaMeta=SigmaMeta) ## p-value for the first SNP
}
