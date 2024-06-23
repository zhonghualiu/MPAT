MinP = function (Z.vec, Sigma)
{
    z.tmp = as.matrix(Z.vec)
    chisq.tmp = z.tmp^2
    pvalue.tmp = pchisq(chisq.tmp, df = 1, lower.tail = FALSE)
    T.minp = min(pvalue.tmp)
    K = length(z.tmp)
    upper_bound = rep(qnorm(1 - T.minp/2), K)
    lower_bound = -rep(qnorm(1 - T.minp/2), K)
    a = pmvnorm(lower = lower_bound , upper = upper_bound, mean = rep(0,
        K), corr = Sigma, algorithm = GenzBretz(maxpts = 50000,
        abseps = 1e-07))
    pval.minp = 1 - a[1]
    return(pval.minp)
}