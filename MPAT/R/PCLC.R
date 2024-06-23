PCLC <-
function(Z.vec,Sigma)
{
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  z.tmp = as.matrix(Z.vec)
  K = length(z.tmp)
  PC.vec = t(eigen.vec)%*%z.tmp
  w.vec = 1/lambdas
  PCLC.stat = sum(w.vec*PC.vec)
  PCLC.stat.std = PCLC.stat^2/sum(w.vec)
  PCLC.p = pchisq(PCLC.stat.std,df=1,lower.tail = FALSE)
  return(PCLC.p)
}
