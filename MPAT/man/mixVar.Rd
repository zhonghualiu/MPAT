\name{mixVar}
\alias{mixVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Z-statistics for multiple phenotypes using inverse variance weighting method in a linear mixed effects model.
}
\description{
Z-statistics for multiple phenotypes are modeled using a linear mixed effects model. Two independent score testing statistics
are combined using inverse variance weighting method. 
}
\usage{
mixVar(Z.vec, Sigma,method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Z.vec}{ Z.vec is a vector representing the Z testing statistics for multiple phenotypes.
}
  \item{Sigma}{
  Sigma is the correlation matrix among the Z testing statistics. }
    \item{method}{
  The method to compute the p-value which takes three values: davies,liu and liumod .}
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
mixVar(Z.vec=lipids_zscore[1,],Sigma) ## p-value for the first SNP
}
