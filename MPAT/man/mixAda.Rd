\name{mixAda}
\alias{mixAda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adaptively select optimal combination coefficient of two independentscore statistics based on linear mixed effects model. 
}
\description{
mixAda combines the Z testing statistics for multiple phenotypes at a genetic variant in an adaptive fashion to maximize the power. }
\usage{mixAda(Z.vec, Sigma)}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Z.vec}{ Z.vec is a vector representing the Z testing statistics for multiple phenotypes.
}
  \item{Sigma}{
  Sigma is the correlation matrix among the Z testing statistics. }
}

\value{
 \item{p-value }{An overall p-value assessing the association between a genetic variant and multiple phenotypes is returned.}
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
mixAda(Z.vec=lipids_zscore[1,],Sigma) ## p-value for the first SNP
}

