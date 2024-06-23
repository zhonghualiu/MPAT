DSUM <-
function(Z.vec,Sigma)
{
    z.tmp =as.matrix(Z.vec) # column vector
    T.SUM = abs(sum(z.tmp))
    pval.SUM= pnorm(T.SUM,mean=0,sd=sqrt(sum(Sigma)),lower.tail=FALSE)*2 # two-sided tests
    return(pval.SUM)
}
