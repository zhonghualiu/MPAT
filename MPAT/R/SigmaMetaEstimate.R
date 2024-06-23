# Estimate the SigmaMeta for PCMeta p-value computation
# This is a pre-computation step for PCMeta
# input: the correlation matrix among Z-scores
# output: the SigmaMeta matrix 
SigmaMetaEstimate = function(Sigma,simNum=1000){

X.PCMinP = rep(NA,simNum)
X.PCFisher = rep(NA,simNum)
X.PCLC = rep(NA,simNum)
X.WI = rep(NA,simNum)
X.Wald =rep(NA,simNum)
X.VC =rep(NA,simNum)

for(i in 1:simNum){
K = dim(Sigma)[1]
Z.vec = t(rmvnorm(1,rep(0,K),Sigma))

X.PCMinP[i] = qnorm(PCMinP(Z.vec,Sigma))
X.PCFisher[i] = qnorm(PCFisher(Z.vec,Sigma))
X.PCLC[i] = qnorm(PCLC(Z.vec,Sigma))
X.WI[i] = qnorm(WI(Z.vec,Sigma))
X.Wald[i] = qnorm(Wald(Z.vec,Sigma))
X.VC[i] = qnorm(VC(Z.vec,Sigma))
}

X.mat = cbind(X.PCMinP,X.PCFisher,X.PCLC,X.WI,X.Wald,X.VC)
X.mat = X.mat[!rowSums(!is.finite(X.mat)),] ## remove rows containing NA, NaN,Inf
SigmaMeta = cor(X.mat)
return(SigmaMeta)

}
  