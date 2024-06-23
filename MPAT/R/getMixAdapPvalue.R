getMixAdapPvalue <-
function(Z,Sigma,B=11){

### Step 1:
B = B
phi.vec = seq(0.01,0.99,length.out=B)

#### Step 2 and 3:
K = length(Z) ## number of traits
J = rep(1,K)
I = diag(K)
H = 1/K * J%*%t(J)
Sigma = Sigma
Sigma.inv = solve(Sigma)

###########################
p.phi.vec=rep(NA,B)
param.mat<-NULL
c1=rep(NA,4)

for(i in 1:B){

    Lambda.phi.mat = phi.vec[i]*Sigma.inv%*%J%*%t(J)%*%Sigma.inv + (1-phi.vec[i])*(I -H)%*%Sigma.inv%*%Sigma.inv%*%(I-H)
    Lambda.phi.mat = 1/2*(Lambda.phi.mat+t(Lambda.phi.mat))
    T.phi = t(Z)%*%Lambda.phi.mat%*%Z
    lamb = eigen(Sigma)$values
    P = eigen(Sigma)$vectors
    Sigma.sqrt = P%*%sqrt(diag(lamb))%*%t(P)
    mat.temp = Sigma.sqrt%*%Lambda.phi.mat%*%Sigma.sqrt
    mat.temp = 1/2*( mat.temp + t(mat.temp ))
    lambda.temp = eigen(mat.temp)$values
    #lambda.temp = eigen(Sigma%*%Lambda.phi.mat)$values
    #p.phi.vec[i] =liu.p.mod(T.phi,lambda=lambda.temp) # could be other methods
    #print(Lambda.phi.mat)
    # print("=================")
    # print(Sigma)
     #print("=================")
     #print(mat.temp)
     #print("=================")
     ## could be solved by rounding to enforce real symmetric positive definite matrix
     #print(lambda.temp) ## could be complex numbers simply due to numerical rounding errors!! very serious in high dimensions!
     #print("==================")
  

    c1[1]<-sum(lambda.temp)
    c1[2]<-sum(lambda.temp^2)
    c1[3]<-sum(lambda.temp^3)
    c1[4]<-sum(lambda.temp^4)
    param<-Get_Liu_Params_Mod(c1)

    DF = param$l
    muQ = param$muQ
    sigmaQ = param$sigmaQ  ### Found this mistake!!! should  be square!!!
    ##note that in central case, muX = l, sigmaX=sqrt(2*l), l=DF
    muX = param$muX
    sigmaX = param$sigmaX

    tstar = (T.phi- muQ)/sigmaQ
    q_approx = tstar*sigmaX + muX
    p.phi.vec[i] =pchisq(q_approx,df=DF,ncp=param$d,lower.tail=F)
    param.mat<-rbind(param.mat,c(muQ,sigmaQ,DF))

}


#### Step 4: Compute related statistics
P.ada = min(p.phi.vec)
q.phi.vec = rep(NA,B)
tau.phi.vec = rep(NA,B)

for(i in 1:B){
    muQ<-param.mat[i,1]
    sdQ<-param.mat[i,2] 
    DF<-param.mat[i,3]

    q.temp = qchisq(1-P.ada,df=DF) ## quantitles for chisq_df=DF
    q.phi.vec[i] = (q.temp - DF)/sqrt(2*DF) *sdQ + muQ ## transform back to the quantiles of T.phi
    tau.phi.vec[i] = phi.vec[i]*t(J)%*%Sigma.inv%*%J

}

#### Step 5: Numerical integration
### be careful about this integrand
## x must a vector and output must be a vector too!!!
integrand.func = function(x){
    #print(x)
    #print(length(x))
    n.x<-length(x)
    # operator %x% represents Kronecker product
    ## temp1 is matrix with n.x columns, and B rows.
	temp1<-tau.phi.vec%x%t(x)
    #print(tau.phi.vec)
    #print(temp1)
    delta.x.temp = (q.phi.vec - temp1)/(1-phi.vec)# matrix of dim, B*n.r
    #print(delta.x.temp)
    delta.x <- apply(delta.x.temp,2,min) # vector of length n.x
    
    # lamb = eigen(Sigma)$values
    # P = eigen(Sigma)$vectors
    # Sigma.sqrt = P%*%sqrt(diag(lamb))%*%t(P)
    # mat.tmp=Sigma.sqrt%*%(I-H)%*%Sigma.inv%*%Sigma.inv%*%(I-H)%*%Sigma.sqrt
    lamb = eigen(Sigma)$values
    P = eigen(Sigma)$vectors
    Sigma.sqrt = P%*%sqrt(diag(lamb))%*%t(P)

    mat.tmp=Sigma.sqrt%*%(I-H)%*%Sigma.inv%*%Sigma.inv%*%(I-H)%*%Sigma.sqrt
    #print(mat.tmp)
    mat.tmp = 1/2*(mat.tmp + t(mat.tmp)) ## to force symmetric to avoid complex eigenvalues. 
    lambda.temp.tau = eigen(mat.tmp)$values
    # print("======================")
    #print(lambda.temp.tau)
    mu.kappa = sum(lambda.temp.tau )
    sigma.kappa = sqrt(2*sum(lambda.temp.tau^2)) ###!!!!
    Ker.kappa <-sum(lambda.temp.tau^4)/(sum(lambda.temp.tau^2))^2 * 12
	Df<-12/Ker.kappa
    #print(Df)
    delta.x.star = (delta.x - mu.kappa)/sigma.kappa*sqrt(2*Df) + Df
    #print(delta.x.star)
	re<-pchisq(delta.x.star ,df=Df) * dchisq(x,df=1)
    #print(re)
    return(re)
}
## minimum value for p-value: 2.151763e-10
# xx = seq(0,40,length=10000)
# plot(xx,integrand.func(xx))
 pvalue= 1 - integrate(integrand.func,lower=0, upper=40, subdivisions=2000,abs.tol = 10^-25)$value

pvalue

}
