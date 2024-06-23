PCLogit <-
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
  
  PC.Logit.stat = sum(log(PC.p/(1-PC.p)))
  PC.Logit.stat.std = abs(PC.Logit.stat)*sqrt(3*(5*K+4)/(K*pi^2*(5*K+2)))
  PC.Logit.p = 2*pt(PC.Logit.stat.std,df=5*K+4,lower.tail = FALSE)
  return(PC.Logit.p)
}
