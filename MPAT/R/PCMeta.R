## input: Z-score: Z.vec; Sigma: the correlation matrix among Z-scores
## SigmaMeta: is the correlation matrix among X, which is inverse normal transformatoin of MinP 
## SigmaMeta: should be pre-computed using simulation 
PCMeta = function(Z.vec,Sigma,SigmaMeta){
 ## avoid Z.vec is a row vector which won't work later
 ## Z.vec should be a column vector object, not a 1*K matrix 
 if(!is.vector(Z.vec)){
 Z.vec = as.vector(Z.vec)
 }

##########
p.PCMinP = PCMinP(Z.vec,Sigma)
p.PCFisher = PCFisher(Z.vec,Sigma)
p.PCLC = PCLC(Z.vec,Sigma)
p.WI = WI(Z.vec,Sigma)
p.Wald = Wald(Z.vec,Sigma)
p.VC = VC(Z.vec,Sigma)
T.minp = min(p.PCMinP,p.PCFisher,p.PCLC, p.WI,p.Wald,p.VC)

K.X = dim(SigmaMeta)[1]
if(K.X!=6){
stop("The dimension of SigmaMeta is not 6!!")
}

upper_bound = rep(Inf,K.X)
lower_bound = rep(qnorm(T.minp), K.X)

a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
        K.X), sigma = SigmaMeta, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
pval = 1 - a[1]
return(pval)
}