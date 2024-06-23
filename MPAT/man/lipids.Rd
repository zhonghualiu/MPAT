\name{lipids}
\alias{lipids}
\docType{data}
\title{
An example data set that contains summary statistics from GWAS studies of four lipids levels. 
}
\description{
This is an example data set that contains summary statistics (Z-scores) for 2000 genetic variants and four lipids levels. 
}
\usage{data(lipids)}
\format{
  A data frame with 2000 observations on the following 15 variables.
  \describe{
    \item{\code{MarkerName}}{SNP names}
    \item{\code{Allele1}}{This the allele used as the effect allele}
    \item{\code{Allele2}}{This is the "other" allele}
    \item{\code{W.HDL}}{The sum of the individual study weights (typically, N) for this marker}
    \item{\code{Zscore.HDL}}{Z scores for HDL}
    \item{\code{Pvalue.HDL}}{P-values for HDL}
    \item{\code{W.LDL}}{The sum of the individual study weights (typically, N) for this marker}
    \item{\code{Zscore.LDL}}{Z scores for LDL}
    \item{\code{Pvalue.LDL}}{P-values for LDL}
    \item{\code{W.TC}}{The sum of the individual study weights (typically, N) for this marker}
    \item{\code{Zscore.TC}}{Z scores for TC}
    \item{\code{Pvalue.TC}}{P-values for TC}
    \item{\code{W.TG}}{The sum of the individual study weights (typically, N) for this marker}
    \item{\code{Zscore.TG}}{Z scores for TG}
    \item{\code{Pvalue.TG}}{P-values for TG}
  }
}

\source{
http://www.sph.umich.edu/csg/abecasis/public/lipids2010/
}
\references{
http://www.nature.com/nature/journal/v466/n7307/full/nature09270.html
}
\examples{
data(lipids)
head(lipids)
}
\keyword{datasets}
