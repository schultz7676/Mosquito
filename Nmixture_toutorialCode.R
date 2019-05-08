source("C:/Users/Jacob/Dropbox/Grad School/2018-2019/Spring/Consulting Class/Git/Mosquito/From Madsen/gen.nmix.R")
ls()

expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=mallard.data
elev<-data[,"elev"]
length<-data[,"length"]
forest<-data[,"forest"]
n.it<-data[,c("count1","count2","count3")]
R<-nrow(n.it)
T<-ncol(n.it)
DATE<-data[,c("date1","date2","date3")]
IVEL<-data[,c("ivel1","ivel2","ivel3")]
DATE2<-DATE^2

#model: Kery, Royle, and Schmid (2005)####
X.const = rep(1,R)
X.lam = cbind(X.const,elev,forest)

p.date.lin = as.vector(t(DATE))
Z.const = rep(1,R*T)
Z.p = cbind(Z.const,p.date.lin)


DATE.vec = c(DATE[,1],DATE[,2],DATE[,3])
stds = sort(na.exclude(unique(DATE.vec)))
diffs = numeric(length(stds)-1)
for(i in 1:(length(stds)-1)){ diffs[i] = stds[i+1] - stds[i] }

diffs.days = round(diffs/min(diffs))
day.unique = numeric(length(diffs.days)+1)

day.unique[1]=1
for(i in 1:(length(diffs.days))){ day.unique[i+1] = day.unique[i] + diffs.days[i] }
day.unique = c(day.unique,NA)
DATE.mat = cbind(DATE.vec,seq(1:length(DATE.vec)))
DATE.sort = DATE.mat[order(DATE.mat[,1]),]
dup = duplicated(DATE.sort[,1])

j=1
for(i in 1:length(DATE.vec)){
  if(dup[i]==TRUE) DATE.sort[i,1]=DATE.sort[i-1,1]
  if(dup[i]!=TRUE) { 
      DATE.sort[i,1]= day.unique[j]
      j=j+1
    }
  }

DATE.mat2 = DATE.sort[order(DATE.sort[,2]),]
DATE.2 = matrix(DATE.mat2,nrow=239,ncol=3,byrow=FALSE)

DATE.3 = ceiling(DATE.2/30)


Method = "BFGS"
K.lim = 40
model.null = optim(c(0,0,0), nmix.mig, method=Method, hessian=TRUE, n=n.it,
                   X=X.const, Z=Z.const, migration="none", prior="NB", Date=DATE.3, K=K.lim)

model.null$conv

ev.null = eigen(model.null$hessian)$values
cn.null = max(ev.null)/min(ev.null)
cn.null #the condition number
model.null$par #this gives the MLE for log(lambda), logit(p), and log(dispersion parameter)
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
                       nmix.mig, method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, migration="none",
                       prior="NB",Date=DATE.3,K=K.lim)
model.closed$conv
ev.closed = eigen(model.closed$hessian)$values
cn.closed = max(ev.closed)/min(ev.closed)
cn.closed #check condition number


#Open Model ####
model.open = optim(c(model.closed$par[1:5],-2,2,model.closed$par[6]), nmix.mig, 
                   method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, migration="constant", 
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
