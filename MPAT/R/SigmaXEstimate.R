# Estimate the SigmaX for PCAQ p-value computation
# This is a pre-computation step for PCAQ 
# input: the correlation matrix among Z-scores
# output: the SigmaX matrix 
SigmaXEstimate = function(Sigma,simNum=1000){

X.WI = rep(NA,simNum)
X.Wald =rep(NA,simNum)
X.VC =rep(NA,simNum)
for(i in 1:simNum){
K = dim(Sigma)[1]
Z.vec = t(rmvnorm(1,rep(0,K),Sigma))
X.WI[i] = qnorm(WI(Z.vec,Sigma))
X.Wald[i] = qnorm(Wald(Z.vec,Sigma))
X.VC[i] = qnorm(VC(Z.vec,Sigma))
}
## if VC gives pvalue=0, then X.VC is inf, then cor function can't compute it 
# remove all rows with non-finite values, this could happen when X.WI is Inf...
X.mat = cbind(X.WI,X.Wald,X.VC)
X.mat = X.mat[!rowSums(!is.finite(X.mat)),]
SigmaX = cor(X.mat)
return(SigmaX)
}