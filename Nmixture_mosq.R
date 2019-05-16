source("./From Madsen/gen.nmix.R")
library(readr)
dat_ready__5_9_2019 <- read_csv("./dat_ready__5_9_2019.csv")
ls()
#Need to define these functions first
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=dat_ready__5_9_2019


#counts
n.it<-data[,paste(rep("y",33), c(1:33), sep="")]
R<-nrow(n.it)
T<-ncol(n.it)

##Site Level covariates
lon<-data[,"lon"]
lat<-data[,"lat"]
total<-data[,"TOTAL"]
#how to deal with habitat..?

#observation-level predictors
temp.obs<-data[,paste(rep("temp.obs",33),c(1:33),sep="")]
temp.max<-data[,paste(rep("temp.max",33),c(1:33),sep="")]
DATE<-data[,paste(rep("t",33),c(1:33),sep="")]
DATE2<-DATE^2


#model: Kery, Royle, and Schmid (2005)####
X.const = rep(1,R) 
X.lam = cbind(X.const,lon,lat,total)

p.date.lin = as.vector(t(DATE))
Z.const = rep(1,R*T)
Z.p = cbind(Z.const,p.date.lin)

#Note: we got to skip a lot of code here because our data
#is not centered and standardized

DATE.2 = DATE

#Now set the assumed primary period length. Our data is only sampled
#every two weeks, so this should probably be small
pim.period.length <- 11
DATE.3 = ceiling(DATE.2/pim.period.length)
DATE.3

#Optimization
#from toutorial: "original N-mixture model with the negative binomial prior and no covariates (ie, the null model) 
#Note:had to use as.matrix(n.it) because it was reading n.it as a list
Method = "BFGS"
K.lim = 50
start.vals = c(.1,.1,.1) #also tried 0,0,0. Both times converged on initial values 
model.null = optim(start.vals, nmix.mig, method=Method, hessian=TRUE, n=as.matrix(n.it),
                   X=X.const, Z=Z.const, migration="none", prior="NB", Date=DATE.3, K=K.lim)

model.null$conv #should come back "0" if we found max

#Evaluate the stability
ev.null = eigen(model.null$hessian)$values
ev.null
cn.null = max(ev.null)/min(ev.null)
cn.null #the condition number

model.null$par #this gives the MLE for log(lambda), logit(p), and log(dispersion parameter alpha)
lambda.est = exp(model.null$par[1])
p.est = expit(model.null$par[2])
c(lambda.est, p.est)

se = sqrt(diag(solve(model.null$hess)))
se #the standard errors

nll = model.null$val
aic = nll + 2*length(model.null$par)

#In fitting the N-mixture model with covariates for lambda and p, the MLEs from the 
#null model are used as initial values.

model.closed = optim(c(model.null$par[1],0,0,model.null$par[2],0,model.null$par[3]),
                     nmix.mig, method=Method, hessian=TRUE, n=as.matrix(n.it), X=X.lam, Z=Z.p, migration="none",
                     prior="NB",Date=DATE.3,K=K.lim)
model.closed$conv
ev.closed = eigen(model.closed$hessian)$values
cn.closed = max(ev.closed)/min(ev.closed)
cn.closed #check condition number


#Open Model ####
model.open = optim(c(model.closed$par[1:5],-2,2,model.closed$par[6]), nmix.mig, 
                   method=Method, hessian=TRUE, n=as.matrix(n.it), X=X.lam, Z=Z.p, migration="constant", 
                   prior="NB", Date=DATE.3, K=K.lim)

model.open$conv
ev.open = eigen(model.open$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others


#Closure Test
t.stat= 2*(model.closed$val - model.open$val)
obs.inf = -1*model.open$hess
obs.inf.nn = obs.inf[1:5,1:5]
obs.inf.np = obs.inf[1:5,6:7]
obs.inf.pn = obs.inf[6:7,1:5]
obs.inf.pp = obs.inf[6:7,6:7]
I.tilda = obs.inf.pp - obs.inf.pn%*%solve(obs.inf.nn)%*%obs.inf.np
prop = acos(I.tilda[1,2]/(sqrt(I.tilda[1,1]*I.tilda[2,2])))/(2*pi)
prop.0 = 0.5 - prop
prop.1 = 0.5
prop.2 = prop
p.value = prop.0*(0) + prop.1*(1-pchisq(t.stat,1)) + prop.2*(1-pchisq(t.stat,2))
p.value

ests.closed = ests(model.closed$par,model.closed$hess, migration="none", n=n.it, X=X.lam, Z=Z.p, T=3, prior="NB")
ests.closed

ests.open = ests(model.open$par, model.open$hess, migration="constant", n=n.it, X=X.lam, Z=Z.p, T=3, prior="NB")
ests.open

gamma = ests.open[9]
gamma
omega = ests.open[10]
omega


se = sqrt(diag(solve(model.open$hess)))
se.gamma = se[6]
se.omega = se[7]

gamma.ci = exp( c( model.open$par[6]- 1.96*se.gamma , model.open$par[6] +1.96*se.gamma))
gamma.ci
omega.ci = expit( c( model.open$par[7] - 1.96*se.omega, model.open$par[7] + 1.96*se.omega))
omega.ci
