\name{PCFisher}
\alias{PCFisher}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fisher's method to combine principal components based on Z testing statistics
}
\description{
Z testing statistics for multiple phenotype are orthogonally transformed onto principal component axes and combined using Fisher's method. 
}
\usage{
PCFisher(Z.vec, Sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Z.vec}{Z.vec is a vector representing the Z testing statistics for multiple phenotypes.
}
  \item{Sigma}{
  Sigma is the correlation matrix among the Z testing statistics. }
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
PCFisher(Z.vec=lipids_zscore[1,],Sigma) ## p-value for the first SNP
}