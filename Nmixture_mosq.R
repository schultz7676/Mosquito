source("./From Madsen/gen.nmix.R")
library(readr)
library("unmarked")
library("tictoc") #So we can see our code's run time

dat_ready__5_9_2019 <- read_csv("./dat_ready__5_9_2019.csv")

#Need to define these functions first
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=dat_ready__5_9_2019
<<<<<<< HEAD
View(data)
data=data[-28,] #eliminates the empty row
#Subset
#timespan of study: # of weeks to include in the dataset
#First year: T in [1,14], Second year: T in [24:33]
=======
data = data[-28,]
>>>>>>> 37701953bb18cf3909195030fd58459392cbe4f4

tspan <- 33

#counts
n.it<-data[,paste(rep("y",tspan), c(1:tspan), sep="")]
R<-nrow(n.it)
T<-ncol(n.it)

##Site Level covariates
lon<-data[,"lon"]
lat<-data[,"lat"]
total<-data[,"TOTAL"]
#how to deal with habitat..?

#observation-level predictors
temp.obs<-data[,paste(rep("temp.obs",tspan),c(1:tspan),sep="")]
temp.max<-data[,paste(rep("temp.max",tspan),c(1:tspan),sep="")]
DATE<-data[,paste(rep("t",tspan),c(1:tspan),sep="")]
DATE2<-DATE^2

#model: Kery, Royle, and Schmid (2005)####
X = cbind(lon,lat,total)

p.date.lin = as.vector(t(DATE))
tobs.lin = as.vector(t(as.matrix(temp.obs)))
Z = cbind(tobs.lin,p.date.lin)

#Note: we got to skip a lot of code here because our data
#is not centered and standardized

DATE.2 = DATE

#Now set the assumed primary period length. Our data is only sampled
#every two weeks, so this should probably be small
pim.period.length <- 1
DATE.3 = ceiling(DATE.2/pim.period.length)
DATE.4 = as.matrix(DATE.3)
DATE.4[is.na(DATE.4)] = max(as.vector(DATE.3),na.rm=TRUE) 
mode(DATE.4) = "integer"

#Optimization
#from toutorial: "original N-mixture model with the negative binomial prior and no covariates (ie, the null model) 
#Note:had to use as.matrix(n.it) because it was reading n.it as a list
<<<<<<< HEAD
Method = "SANN"
K.lim = 100
start.vals = c(.5,0) 
#Testing Likelihood
nmix.mig(start.vals,ceiling(as.matrix(n.it)/100),X.const,Z.const,Date=DATE.3,K=K.lim)

model.null = optim(start.vals, nmix.mig, method=Method, hessian=TRUE, n=as.matrix(n.it), X=X.const, Z=Z.const, migration="none", prior="poisson", Date=DATE.3, K=K.lim,control=list(trace=2,reltol=1e-2,ndeps=c(50,.05)))

model.null$conv #should come back "0" if we found max
=======
Method = "BFGS"
K.lim = 200
start.vals = c(0,0)
y = unmarkedFramePCount(ceiling(as.matrix(n.it)/100),
						siteCovs=NULL,
						obsCovs=NULL,
						mapInfo=NULL)
model.null = pcount(~1~1,
					y,
					K.lim,
					mixture="P",
					start.vals,
					method=Method,
					se=TRUE,
					engine="C")

summary(model.null)
>>>>>>> 37701953bb18cf3909195030fd58459392cbe4f4

#Evaluate the stability
ev.null = eigen(model.null@opt$hessian)$values
cn.null = max(ev.null)/min(ev.null)
cn.null #the condition number. Should not be negative or close to 0.


#In fitting the N-mixture model with covariates for lambda and p, the MLEs from the 
#null model are used as initial values.
y = unmarkedFramePCount(ceiling(as.matrix(n.it)/100),
						siteCovs=as.data.frame(X),
						obsCovs=as.data.frame(Z),
						mapInfo=NULL)
model.closed = pcount(~1 ~1,
					  y,
					  K.lim,
					  mixture="P",
					  model.null@opt$par,
					  method=Method,
					  se=TRUE,
					  engine="C")

summary(model.closed)

ev.closed = eigen(model.closed@opt$hessian)$values
cn.closed = max(ev.closed)/min(ev.closed)
cn.closed #check condition number


#Open Model #### This takes a long time to run, but it does seem to converge.
y = unmarkedFramePCO(ceiling(as.matrix(n.it)/100),
					 siteCovs=as.data.frame(X),
					 obsCovs=as.data.frame(Z),
					 yearlySiteCovs=NULL,
					 mapInfo=NULL,
					 numPrimary=33,
					 primaryPeriod=DATE.4)
tic()
model.open = pcountOpen(~1,~1,~1,~1,
						y,
						mixture="P",
						K.lim,
						dynamics="constant",
						fix="none",
						starts=c(model.closed@opt$par[1],0,0,model.closed@opt$par[2]),
						method=Method,
						se=TRUE,
						immigration=FALSE,
						iotaformula=~1)
toc() #714.39 sec on Jacob's machine (12 mins)

summary(model.open)

lam <- exp(coef(model.open, type="lambda"))
gam <- exp(coef(model.open, type="gamma"))
om <- plogis(coef(model.open, type="omega"))
p <- plogis(coef(model.open, type="det"))
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open@opt$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others


#Closure Test
if(FALSE){
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
}
