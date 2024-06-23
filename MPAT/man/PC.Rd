\name{PC}
\alias{PC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PC is a function used to combine correlated GWAS summary statistics for multiple phenotypes at a genetic locus. 
}
\description{
By combining univariate summary statistics for multiple phenotypes, we can obtain an overall association p-value that assess the association between 
a genetic variant and multiple phenotypes as a whole. PC allows users to conduct dimension reduction by specifying which principal component to use to represent 
the original multiple summary statistics. 
}
\usage{
PC(Z.vec, Sigma, PCorder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z.vec}{ Z.vec is a column vector representing the Wald-type Z testing statistics for multiple phenotypes.
}
  \item{Sigma}{
       Sigma is the correlation matrix of the multiple Z testing statistics. 
}
  \item{PCorder}{
PCorder specifys which principal component to use, it takes integer values. For example, PCorder=1 means we choose to use the first principal component. 
}
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
PC(Z.vec=lipids_zscore[1,],Sigma,PCorder=1) ## p-value for the first SNP using PC1
}


