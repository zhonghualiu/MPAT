PCSS <-
function(Z.vec,Sigma)
{ 
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  z.tmp = as.matrix(Z.vec)
  K = length(z.tmp)
  PC.vec = t(eigen.vec)%*%z.tmp
  w.vec = 1/lambdas
  PCSS.stat = sum(w.vec*PC.vec^2)
  PCSS.p = pchisq(PCSS.stat,df=K,lower.tail = FALSE)
  return(PCSS.p)
}
