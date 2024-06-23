Wald <-
function(Z.vec,Sigma)
{
    K = length(Z.vec)
    R.inv = solve(Sigma) ## inverse
    z.tmp =as.matrix(Z.vec) # column vector
    T.Wald = t(z.tmp)%*%R.inv%*%z.tmp # Quadratic form
    pval.Wald = pchisq(T.Wald,df=K,lower.tail=FALSE) # 
    return(pval.Wald) 
}
