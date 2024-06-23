mixAda <-
function(Z.vec, Sigma)
{
    Sigma.inv = solve(Sigma)
    K = length(Z.vec)
    J = rep(1,K)
    I = diag(K)
    H = 1/K * J%*%t(J)
    z.tmp = as.matrix(Z.vec) # column vector

    pval.mixAda = getMixAdapPvalue(z.tmp,Sigma=Sigma,B=11) ## B is the number of grid between 0 and 1
    return(pval.mixAda)
}
