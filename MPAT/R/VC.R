VC <-
function(Z.vec,Sigma,method="davies")
{
    R.inv = solve(Sigma)
    lambdas = eigen(R.inv)$values

    z.tmp =as.matrix(Z.vec) # column vector
    T.VC = t(z.tmp)%*%R.inv%*%R.inv%*%z.tmp # Quadratic form
	if(method=="davies"){
	pval.VC = davies(T.VC,lambdas)$Qq
	if(pval.VC == 0 | pval.VC < 0){
	pval.VC = liumod(T.VC,lambdas)
	}
	}
	else if(method=="liu"){
    pval.VC = liu(T.VC,lambdas) # Liu method to compute p-value
	}
	else if(method=="liumod"){
	pval.VC = liumod(T.VC,lambdas)
	}
	else {stop("Please specify a valid method: davies, liu, liumod.\n")}
    return(pval.VC)
}
