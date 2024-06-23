PCMinP <-
function(Z.vec,Sigma)
{
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  z.tmp = Z.vec
  K = length(z.tmp)
  PC.vec = t(eigen.vec)%*%z.tmp
  PC.std = (PC.vec)^2/lambdas
  PC.p = pchisq(PC.std,df=1,lower.tail = FALSE)
  p.min = min(PC.p)
  PCMinP.p = 1- (1-p.min)^K
  return(PCMinP.p)
}
