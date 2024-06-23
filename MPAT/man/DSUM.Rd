\name{DSUM}
\alias{DSUM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Direct summation of Z-statistics
}
\description{
For SNP, the Z testing statistics for multiple phenotypes are summed together as an overall testing statistic.  
}
\usage{
DSUM(Z.vec, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z.vec}{ Z.vec is a vector representing the Z testing statistics for multiple phenotypes.}
  \item{Sigma}{ Sigma is the correlation matrix among the Z testing statistics that can be estimated by its sample version over the SNPs.}
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
DSUM(Z.vec=lipids_zscore[1,],Sigma) ## p-value for the first SNP
}
