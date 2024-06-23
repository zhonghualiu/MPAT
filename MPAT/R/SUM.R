SUM <-
function(Z.vec,Sigma)
{
    z.tmp =as.matrix(Z.vec) # column vector
    K = length(Z.vec)
    J = as.matrix(rep(1,K))
    R.inv = solve(Sigma)
    T.SUM = abs(t(J)%*%R.inv%*%z.tmp)
    pval.SUM= pnorm(T.SUM,mean=0,sd=sqrt(sum(R.inv)),lower.tail=FALSE)*2 # two-sided tests
    return(pval.SUM)
}
