WI <-function(Z.vec,Sigma,method = "davies")
{  
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  
  z.tmp = as.matrix(Z.vec)  # column vector
  T.WI = t(z.tmp) %*%z.tmp # Quadratic form
  if(method=="liu"){
  pval.WI = liu(T.WI,lambdas)
  }
  else if(method=="liumod"){
  pval.WI = liumod(T.WI,lambdas)

  }
  else if(method=="davies"){
  pval.WI = davies(T.WI,lambdas)$Qq
  if(pval.WI ==0 | pval.WI < 0){
  pval.WI = liumod(T.WI,lambdas)
  }
  }
  else {stop("Please specify a valid method: davies, liu, liumod.\n")}
  return(pval.WI)
}


