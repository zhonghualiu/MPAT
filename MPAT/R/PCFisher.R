PCFisher <-
function(Z.vec,Sigma)
{
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  z.tmp = as.matrix(Z.vec)
  K = length(z.tmp)
  PC.vec = t(eigen.vec)%*%z.tmp
  PC.std = (PC.vec)^2/lambdas
  PC.p = pchisq(PC.std,df=1,lower.tail=FALSE)
  PC.Fisher.stat = -2*sum(log(PC.p))
  PC.Fisher.p = pchisq(PC.Fisher.stat,df=2*K,lower.tail = FALSE)
  return(PC.Fisher.p)
}
