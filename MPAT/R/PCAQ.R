## input: Z-score: Z.vec; Sigma: the correlation matrix among Z-scores
## SigmaX: is the correlation matrix among X, which is inverse normal transformatoin of MinP 
## SigmaX: should be pre-computed using simulation 
PCAQ = function(Z.vec,Sigma,SigmaX){
 ## avoid Z.vec is a row vector which won't work later
 ## Z.vec should be a column vector object, not a 1*K matrix 
 if(!is.vector(Z.vec)){
 Z.vec = as.vector(Z.vec)
 }


p.WI = WI(Z.vec,Sigma)
p.Wald = Wald(Z.vec,Sigma)
p.VC = VC(Z.vec,Sigma)
##take minimum
T.minp = min(c(p.WI,p.Wald,p.VC))
## use multivariate normal CDF to compupte p-value

K.X = dim(SigmaX)[1]
if(K.X!=3){
stop("The dimension of SigmaX is not 3!!")
}

upper_bound = rep(Inf,K.X)
lower_bound = rep(qnorm(T.minp), K.X)

a = pmvnorm(lower = lower_bound, upper = upper_bound, mean = rep(0, 
        K.X), sigma = SigmaX, algorithm = GenzBretz(maxpts = 250000, abseps = 1e-13))
p.PCAQ = 1 - a[1]
return(p.PCAQ)
}