library(unmarked)
library(jagsUI)
library(MASS)
library(ggmcmc)
##########################################
##########################################
### Single-season occupancy simulation ###
##########################################
########################################## 
set.seed(567)
M<-500 # Number of sites
J<-3 # Number of repeat surveys per site
a0<-0.5 # True occupancy intercept (logit scale)
a1<-1 # True occupancy slope (logit scale)
a2<--0.5 # True occupancy slope (logit scale)
X1<-rnorm(M) # site covariate 1 (mean 0, sd 1)
X2<-rbinom(M,1,0.35) # site covariate 3 (binary-coded)
psi<-plogis(a0 + a1*X1 + a2*X2) # true occupancy probability
z<-rbinom(M,1,psi) # true (latent) occurrence (z = 1 is present; z = 0 is absent)
y<-array(NA,dim=c(M,J))
X3<-array(NA, dim=c(M,J)) # empty detection covariate (X4) matrix
p<-array(NA, dim = c(M,J)) # empty detection probability matrix
b0<-1 # true detection intercept (logit scale)
b1<-0.5 # true slope associated with X1 (logit scale)
b2<--2 # true slope associated with X4 (logit scale)
for (i in 1:M){
  for (j in 1:J){
    X3[i,j]<-runif(1,-1,1)
    p[i,j]<-plogis(b0 + b1*X1[i] + b2*X3[i,j])
    y[i,j]<-rbinom(1,1,p[i,j]*z[i]) # survey data matrix
    }
  }
occ.umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(X1 = X1, X2 = X2), obsCovs = list(X3 = X3))
SSoccunmarked<-occu(~X1 + X3 ~X1 + X2, data = occ.umf) #Order of formula is ~Detection ~Occupancy
summary(SSoccunmarked)

# Sum up number of occupied sites
OccSites<-ranef(SSoccunmarked)
sum(ranef(SSoccunmarked)@post[,2,])

####################################################################
### Bayesian implementation in JAGS (just another gibbs sampler) ###
####################################################################  

zst<-apply(y,1,max) # initial values for latent occupancy state z
jags.SS.data<-list(y=y, nSites=M, nOccasions = J, X1=X1, X2=X2, X3=X3)

cat(file="jagsSSocc.txt", "
model{

  	for (i in 1:nSites){
		z[i]~dbern(psi[i])
		logit(psi[i])<- alpha[1] + alpha[2]*X1[i] + alpha[3]*X2[i]
		}
	for (i in 1:nSites){
		for (j in 1:nOccasions){
			logit(p[i,j]) <- beta[1] + beta[2]*X1[i] + beta[3]*X3[i,j]
			muy[i,j]<- z[i]*p[i,j]
			y[i,j]~dbern(muy[i,j]) # observed data
				}
			}
for (j in 1:3){
	alpha[j]~dnorm(0,0.37)
}
for (j in 1:3){
	beta[j]~dnorm(0,0.37)
}
occupiedSites<-sum(z[]) # sum up number of occupied sites
}")

inits<-function()list(alpha=c(runif(1,-3,3),runif(1,-3,3),runif(1,-3,3)), beta=c(runif(1,-3,3),runif(1,-3,3),runif(1,-3,3)), z=zst)
params=c("alpha", "beta", "occupiedSites")
ni<-1000
na<-500
nb<-500
nc<-3
nt<-1

# run a model with parallel processing
SSoccjag<-jags(jags.SS.data, inits, params, "jagsSSocc.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.adapt=na, n.burnin = nb, parallel=T)
SSoccjag
summary(SSoccunmarked)

# diagnostics and plots (prints a pdf to your working directory)
ggmcmc(ggs(SSoccjag$samples))