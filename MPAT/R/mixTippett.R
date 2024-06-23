mixTippett <-
function(Z.vec,Sigma,method="davies")
{   
    z.tmp = as.matrix(Z.vec) ## transform to be a column vector 
    R.inv = solve(Sigma)
    K = length(Z.vec)
    J = rep(1,K)
    I = diag(K)
    H = 1/K * J%*%t(J)
  
    ######## compute U.mu
    U.mu = (t(J)%*%R.inv%*%z.tmp)^2/(t(J)%*%R.inv%*%J) ## follow chisq df =1
    pval.mu = pchisq(U.mu,df=1,lower.tail=FALSE) # compute p-value for U.mu0
    
    #  compute U.tau
    mu.hat = mean(z.tmp)
    U.tau = t(z.tmp-mu.hat)%*%R.inv%*%R.inv%*%(z.tmp-mu.hat)
    
    lamb = eigen(Sigma)$values
    P = eigen(Sigma)$vectors
    Sigma.sqrt = P%*%sqrt(diag(lamb))%*%t(P)
    mat.tmp = Sigma.sqrt%*%(I-H)%*%R.inv%*%R.inv%*%(I-H)%*%Sigma.sqrt
    mat.tmp = 1/2*(mat.tmp + t(mat.tmp))
    lambda.tmp = eigen(mat.tmp)$values
    
	
	if(method=="davies"){
    pval.tau = davies(U.tau,lambda.tmp)$Qq
	  if(pval.tau ==0 | pval.tau <0){
      	pval.tau = liumod(U.tau,lambda.tmp) 
	  }  
	}
	else if(method=="liu"){
	pval.tau = liu(U.tau,lambda.tmp)
	}
	else if(method=="liumod"){
	pval.tau = liumod(U.tau,lambda.tmp)
	}
	else {stop("Please specify a valid method: davies, liu, liumod.\n")}
	
    # Tippet method to combine pval.mu and pval.tau
    pval.mixTippett = 1 - (1-min(pval.mu,pval.tau))^2
    
    p_group = pval.mu
    p_individual = pval.tau
    p_overall = pval.mixTippett
    out = c(p_overall,p_group,p_individual)
    names(out) = c("p_overall","p_group","p_individual")
    return(out)    
}
