mixSD <-
function(Z.vec,Sigma,method="davies")
{
    z.tmp = as.matrix(Z.vec)
    R.inv = solve(Sigma)
    l = eigen(R.inv)$values
    
    K = length(z.tmp)
    J = rep(1,K)
    I = diag(K)
    H = 1/K * J%*%t(J)
    # compute inverse SD  weight
    varUsq = 2*sum(diag(R.inv%*%J%*%t(J)%*%R.inv%*%J%*%t(J)))
    Lambda.tau = (I-H)%*%R.inv%*%R.inv%*%(I-H)
    varTau = 2*sum(diag(Lambda.tau%*%Sigma%*%Lambda.tau%*%Sigma))
    phi = sqrt(varTau)/(sqrt(varUsq)+sqrt(varTau))

    Lambda.phi.mat = phi*R.inv%*%J%*%t(J)%*%R.inv + (1-phi)*(I -H)%*%R.inv%*%R.inv%*%(I-H)
    l1= eigen(Lambda.phi.mat)$values

    
    mat.tmp = 1/2*( Lambda.phi.mat + t(Lambda.phi.mat))
    T.phi = t(z.tmp)%*%mat.tmp%*%z.tmp
    
    lamb = eigen(Sigma)$values
    P = eigen(Sigma)$vectors
    Sigma.sqrt = P%*%sqrt(diag(lamb))%*%t(P)
    mat.tmp2 = Sigma.sqrt%*%mat.tmp%*%Sigma.sqrt
    mat.tmp2 = 1/2*(mat.tmp2 + t(mat.tmp2))
    lambda.tmp = eigen(mat.tmp2)$values
	if(method=="liu"){
    pval.mixSD = liu(T.phi,lambda.tmp)
	}
	else if(method=="davies"){
	pval.mixSD = davies(T.phi,lambda.tmp)$Qq
	if(pval.mixSD == 0 | pval.mixSD < 0){
     pval.mixSD = liumod(T.phi,lambda.tmp)
	}
	}
	else if(method=="liumod"){
	pval.mixSD = liumod(T.phi,lambda.tmp)
	}
	else {stop("Please specify a valid method: davies, liu, liumod.\n")}
	return(pval.mixSD)
}
