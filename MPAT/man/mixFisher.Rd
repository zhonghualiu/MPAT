\name{mixFisher}
\alias{mixFisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Z-statistics for multiple phenotypes using Fisher's method in a linear mixed effects model.  
}
\description{
Z-statistics for multiple phenotypes are modeled using a linear mixed effects model. Two independent score testing statistics
are combined using Fisher's method. 
}
\usage{
mixFisher(Z.vec, Sigma,method)
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
\item{p_group}{p-value for the group effects of a genetic variant on the multiple phenotypes}
\item{p_individual}{p-value for the individual effects of a genetic variant on the multiple phenotypes}
\item{p_overall}{An overall p-value assessing the association between a genetic variant and multiple phenotypes is returned.}
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
mixFisher(Z.vec=lipids_zscore[1,],Sigma) ## p-value for the first SNP
}