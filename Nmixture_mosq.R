source("./From Madsen/gen.nmix.R")
library(readr)
library("unmarked")
dat_ready__5_9_2019 <- read_csv("./dat_ready__5_9_2019.csv")

#Need to define these functions first
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=dat_ready__5_9_2019
data = data[-28,]


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

#Evaluate the stability
ev.null = eigen(model.null@opt$hessian)$values
cn.null = max(ev.null)/min(ev.null)
cn.null #the condition number


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

summary(model.open)

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
