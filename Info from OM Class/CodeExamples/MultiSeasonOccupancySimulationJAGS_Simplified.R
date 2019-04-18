library(unmarked)
library(jagsUI)
library(MASS)
library(ggmcmc)
##############################################################
### First set up some required arrays and parameter values ###
##############################################################

set.seed(567)
# Specify the number of Sites, Years, and Occasions ###
R<-500 # number of sites
J<-5   # number of years (primary occasions)
K<-3   # occasions (secondary occasions)

# Initial occupancy
a0<--0.25 # initial occupancy intercept
a1<-0.50  # slope for site-level occupancy covariate
X1<-runif(R,-3,3) # generate some random values for the site-level covariate, X1 (standardized)
psi<-plogis(a0 + a1*X1) # calculate initial site-specific occupancy probabilities

# Detection
b0<-1.00 # detection intercept
b1<-0.50 # slope for site-level detection covariate
b2<--0.20 # slope for observation-level covariate
X2<-array(dim = c(R,K,J)) # array for site-level detection covariate, X2 (standardized)
### Create occasion-specific detection covariate ###
for (i in 1:R){
    for (k in 1:K){
      for (j in 1:J){
        X2[i,k,j]<-runif(1,-3,3)
      }
   }
}
# Year-specific, site-level covariates X3 and X4, which are slightly negatively correlated with mean 0, variance = 1, covariance = -0.3. Only X3 affects colonization and extinction
yearcovs <- as.data.frame(mvrnorm(R*J, mu = c(0,0), 
                     Sigma = matrix(c(1,-0.6,-0.6,1), ncol = 2), 
                     empirical = TRUE))
names(yearcovs)<-c("X3","X4")
X3<-array(yearcovs$X3,dim=c(R,J)) # year-specific, site-level covariate X3
X4<-array(yearcovs$X4,dim=c(R,J)) # year-specific, site-level covariate X4

# plot correlated covariates
plot(X3,X4)

# Colonization
g0<--1.00 # colonization intercept
g1<- 0.50 # slope for site-level colonization covariate
g2<--1.50 # slope for year-specific site-level colonization covariate
gamma<-array(dim = c(R,J-1)) # array for colonization probabilties

# Extinction
e0<--1.25 # extinction intercept
e1<- 0.50 # slope for site-level extinction intercept
e2<--0.50 # slope for year-specific site-level covariate
epsilon<-array(dim = c(R,J-1)) # array for extinction probabilities

# Additional arrays/matrices, used below
z<-array(dim=c(R,J)) # latent occupancy state, indexed by Site and Year
muZ<-array(dim=c(R,J)) # latent occupancy state, indexed by Site and Year
y<-array(dim = c(R,K,J)) # array for detection/non-detection data, indexed by Site, Occasion, Year
p<-array(dim = c(R,K,J)) # detection probability array (indexed by site, occasion, year)muZ<-array(dim=c(R,J))
prob<-array(dim = c(R,K,J)) # array for effective detection probability (z*p)

####################################################
####################################################
### Now simulate observation and state processes ###
####################################################
#################################################### 

# First calculate occasion-specific p
for (i in 1:R){
    for (k in 1:K){
      for (j in 1:J){
        p[i,k,j]<-plogis(b0 + b1*X1[i] + b2*X2[i,k,j])
    }
  }
}
# Calculate site and year-specific colonization (gamma) and extinction (epsilon) probabilities. Note that although there are only 4 transitions (4 intervals over which colonization and extinction events occur). Unmarked requires a site by year matrix (R x K), which in this case is a 250 x 5 matrix. Also note that when you fit the model, unmarked will assume that the first column corresponds with the first interval, and so on, and that the final column is ignored....just something to keep in mind when fitting a model with real data and formatting things.   
for (i in 1:R){
  for (j in 1:(J-1)){
      gamma[i,j]<-plogis(g0 + g1*X1[i] + g2*X3[i,j])
      epsilon[i,j]<-plogis(e0 + e1*X1[i] + e2*X4[i,j])
    }
  }
# Determine intitial occupancy state (latent occupancy z)
z[,1]<-rbinom(R,1,psi)

# Simulate state process: temporal changes in occupancy conditional upon intital state
for (i in 1:R){
  for (j in 2:J){
  muZ[i,j]<-z[i,j-1]*(1-epsilon[i,j-1]) + (1-z[i,j-1])*gamma[i,j-1]
  z[i,j]<-rbinom(1,1,muZ[i,j])
   }
}
# Simulate detection/nondetection data (i.e., fake survey data)
for (i in 1:R){
    for (k in 1:K){
      for (j in 1:J){
        prob[i,k,j]<-z[i,j]*p[i,k,j]
        y[i,k,j]<-rbinom(1,1,prob[i,k,j])
    }
  }
}
# sitecovs data frame (we only have one site-level covariate, but we could have more)
sitecovs<-data.frame(X1)

# Reformat observation data and convert to data frame for unmarked; there's probably a more elegant way to convert a 3D array to a data frame ###
yfin<-matrix(y,nrow(y),J*K, byrow = F)

# Reformat observation-level covariates (influencing p) as an R X T*J matrix
X2fin<-matrix(X2,nrow(X2),J*K, byrow = F)

# Creat unmarkedMultFrame 
SimDat<-unmarkedMultFrame(y = yfin,numPrimary=J,siteCovs=sitecovs,obsCovs=list(X2=X2fin),yearlySiteCovs=list(X3=X3,X4=X4))

# save(SimDat, file = "SimDat.RData")
# load(file = "SimDat.RData")

#########################################
#########################################
### Run model using "colext" function ###
#########################################
######################################### 
occmultunm<- colext(psiformula=~X1,gammaformula=~X1+X3,epsilonformula=~X1+X4,pformula=~X1+X2,data=SimDat)
summary(occmultunm)

####################################################################
####################################################################
### Bayesian implementation in JAGS (just another gibbs sampler) ###
####################################################################  
####################################################################

# To run the mdoel in JAGS we need to make sure the arrays are in the proper format, because
# the JAGS model is indexed slightly differently from the R simulation (I could have just re-written the JAGS model). 
# Also, the aperm function enables you to re-arrange a multi-dimensional array.
# Here, I'm taking the 3rd dimension of y and changing it to the first dimension, and so on. What I want is a T x R x K
# array for y and bX in jags.
yjags<-aperm(y, c(3,1,2)) 
X2jags<-aperm(X2, c(3,1,2))
# Here I want a T x R array for X3 and X4
X3jags<-aperm(X3, c(2,1))
X4jags<-aperm(X4, c(2,1))

# JAGS is picky about initial values for the latent occupancy state, so you have to supply reasonable initial values; here, we're just supplying the max observation across sites and occasions for each year, so we preserve dimentions 1 and 3 (i.e. taking max across occasions)
zst<-apply(y,c(1,3),max) # have to supply good initial values
zst<-aperm(zst,c(2,1)) # then re-arrange it a bit to fit our model as it's written

# now package all of the data up. Notice that we need to lag the site x time covariates by 1; otherwise, they'd start at the second column and read to the 5th column. 
jags.mult.data<-list(y=yjags, nSites=R, nYears=J, nOccasions = K, X1=sitecovs$X1, X2=X2jags, X3=X3jags, X4=X4jags)

cat(file="jagsdynocc.txt", "
model{
for (t in 1:1){
	for (i in 1:nSites){
		z[1,i]~dbern(psi[1,i])
		logit(psi[1,i])<- alpha[1] + alpha[2]*X1[i] 
		}
}
for (t in 2:nYears){
	for (i in 1:nSites){
		z[t,i]~dbern(muz[t,i])
		muz[t,i] <- z[t-1,i]*(1-epsilon[t-1,i]) + (1-z[t-1,i])*gamma[t-1,i]
		logit(epsilon[t-1,i])<- eps[1] + eps[2]*X1[i] + eps[3]*X4[t-1,i] 
		logit(gamma[t-1,i])<- gam[1] + gam[2]*X1[i] + gam[3]*X3[t-1,i]
	}
}
for (t in 1:nYears){
	for (i in 1:nSites){
		for (j in 1:nOccasions){
			logit(p[t,i,j]) <- beta[1] + beta[2]*X1[i] + beta[3]*X2[t,i,j]
			muy[t,i,j]<- z[t,i]*p[t,i,j]
			y[t,i,j]~dbern(muy[t,i,j]) # observed data
				}
			}
	   }
for (j in 1:2){
	alpha[j]~dnorm(0,0.37)
}
for (j in 1:3){
	eps[j]~dnorm(0,0.37)
}
for (j in 1:3){
	gam[j]~dnorm(0,0.37)
}
for (j in 1:3){
	beta[j]~dnorm(0,0.37)
}
for (t in 1:nYears){
OccupiedSites[t]<- sum(z[t,])
}
}")

inits<-function()list(alpha=c(runif(1,-3,3),runif(1,-3,3)), beta=c(runif(1,-3,3),runif(1,-3,3),runif(1,-3,3)), gam = c(runif(1,-3,3),runif(1,-3,3),runif(1,-3,3)), eps=c(runif(1,-3,3),runif(1,-3,3),runif(1,-3,3)), z=zst)

params=c("alpha", "beta", "eps", "gam", "OccupiedSites")
ni<-1000
nt=1
na = 500
nb=500
nc=1

occmultjags<-jags(jags.mult.data, inits, params, "jagsdynocc.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=F)

occmultjags
summary(occmultunm)
ggmcmc(ggs(occmultjags$samples))