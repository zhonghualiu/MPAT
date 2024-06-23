\name{SigmaMetaEstimate}
\alias{SigmaMetaEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To estimate the correlation matrix of among X used for PCMeta
}
\description{
Compute SigmaX which is required for PCMeta function. 
}
\usage{
SigmaMetaEstimate(Sigma,simNum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Sigma}{
  Sigma is the correlation matrix among the Z testing statistics. }
  
 \item{simNum}{ The number of simulations to be performed to estimate SigmaX.
}

   
}
\value{
\item{SigmaMeta}{A correlation matrix among X.}
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
SigmaMetaEstimate(Sigma,simNum=1000) 
}
