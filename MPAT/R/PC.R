PC <-
function(Z.vec,Sigma,PCorder) # Z.vec must be a column vector 
{
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  z.tmp = as.matrix(Z.vec) ## must be a column vector D
  
  lambda = lambdas[PCorder]
  u = eigen.vec[,PCorder]
  PC = t(u)%*%z.tmp
  PC.std = (PC/sqrt(lambda))^2
  PC.p = pchisq(PC.std,df=1,lower.tail = FALSE) ## use lower.tail can have better accuray than 1 - pchisq
  return(PC.p)
}
