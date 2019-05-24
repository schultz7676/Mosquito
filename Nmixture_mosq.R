source("./From Madsen/gen.nmix.R")
library(readr)
library("unmarked")
library("tictoc") #So we can see our code's run time

dat_ready__5_9_2019 <- read_csv("./dat_ready__5_9_2019.csv")

#Need to define these functions first
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=dat_ready__5_9_2019
#View(data)

data=data[-28,] #eliminates the empty row
#Subset

tspan <- 33
#timespan of study: # of weeks to include in the dataset
#First year: T in [1,14], Second year: T in [24:33]

#counts
n.it<-data[,paste(rep("y",tspan), c(1:tspan), sep="")]
R<-nrow(n.it)
T<-ncol(n.it)

##Site Level covariates
lon<-data[,"lon"]
lat<-data[,"lat"]
total<-data[,"TOTAL"]
trap_num<-data[,"trap_num"]
#how to deal with habitat..?

#observation-level predictors
temp.obs<-data[,paste(rep("temp.obs",tspan),c(1:tspan),sep="")]
temp.max<-data[,paste(rep("temp.max",tspan),c(1:tspan),sep="")]
DATE<-data[,paste(rep("t",tspan),c(1:tspan),sep="")]
DATE2<-DATE^2

#model: Kery, Royle, and Schmid (2005)####
X = cbind(lon,lat,total,trap_num)
#X = cbind(lon,lat,total)


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

#Closed model, binned responses ####
Method = "BFGS"
K.lim = 200
start.vals = c(0,0)
y = unmarkedFramePCount(ceiling(as.matrix(n.it)/100),
						siteCovs=as.data.frame(X),
						obsCovs=as.data.frame(Z),
						mapInfo=NULL)
model.closed = pcount(~1 ~1,
					  y,
					  K.lim,
					  mixture="P",
					  start.vals,
					  method=Method,
					  se=TRUE,
					  engine="C")

summary(model.closed)

ev.closed = eigen(model.closed@opt$hessian)$values
cn.closed = max(ev.closed)/min(ev.closed)
cn.closed #check condition number

(lam <- 100*exp(model.closed@opt$par[1]))
(p <- plogis(coef(model.closed, type="det")))
c(lam,p)
#Open intercept Model, binned responses + closure test #### 

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
						starts=c(model.closed@opt$par[1],log(.4),logit(.4),model.closed@opt$par[2]),
						method=Method,
						se=TRUE,
						immigration=FALSE,
						iotaformula=~1)
toc() #714.39 sec on Jacob's machine (12 mins)

#note: starting at log(.4), logit(.4) instead of 0, 0 
#saved one optim interation 59->58. Not a huge improvement

summary(model.open)

lam <- 100*exp(coef(model.open, type="lambda")) # *100 because of binned responses
gam <- exp(coef(model.open, type="gamma"))
om <- plogis(coef(model.open, type="omega"))
p <- plogis(coef(model.open, type="det"))
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open@opt$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others

#closure test
#method 1: extract log-likelihoods from fitted models
-2*(model.open@opt$value-model.closed@opt$value)
model.open2@opt$value
#method 2: use unmarked built in LRT
LRT(model.open,model.closed)
#Result in test stat of 2765.518 in both cases

#Open Model covariates, binned responses ####

y = unmarkedFramePCO(ceiling(as.matrix(n.it)/100),
                     siteCovs=as.data.frame(X),
                     obsCovs=as.data.frame(Z),
                     yearlySiteCovs=NULL,
                     mapInfo=NULL,
                     numPrimary=33,
                     primaryPeriod=DATE.4)
tic()
model.open.siteLam = pcountOpen(~as.factor(trap_num),~1,~1,~1,
                        y,
                        mixture="P",
                        K.lim,
                        dynamics="constant",
                        fix="none",
                        starts=c(model.closed@opt$par[1],rep(log(.4),63),logit(.4),model.closed@opt$par[2]),
                        method=Method,
                        se=TRUE,
                        immigration=FALSE,
                        iotaformula=~1)
toc()
summary(model.open.siteLam)

lam <- 100*exp(coef(model.open.siteLam, type="lambda"))
gam <- exp(coef(model.open.siteLam, type="gamma"))
om <- plogis(coef(model.open.siteLam, type="omega"))
p <- plogis(coef(model.open.siteLam, type="det"))
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open.siteLam@opt$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others

#Open intercept model without binning responses####
K.lim.raw = max(as.matrix(n.it),na.rm=TRUE)+1

y = unmarkedFramePCO(as.matrix(n.it),
                     siteCovs=as.data.frame(X),
                     obsCovs=as.data.frame(Z),
                     yearlySiteCovs=NULL,
                     mapInfo=NULL,
                     numPrimary=33,
                     primaryPeriod=DATE.4)

tic()
model.open.raw = pcountOpen(~1,~1,~1,~1,
                        y,
                        mixture="P",
                        K.lim2,
                        dynamics="constant",
                        fix="none",
                        starts=c(model.closed@opt$par[1],log(.4),logit(.4),model.closed@opt$par[2]),
                        method=Method,
                        se=TRUE,
                        immigration=FALSE,
                        iotaformula=~1)
toc()

summary(model.openraw)


summary(model.open.raw)

lam <- 100*exp(coef(model.open.raw, type="lambda"))
gam <- exp(coef(model.open.raw, type="gamma"))
om <- plogis(coef(model.open.raw, type="omega"))
p <- plogis(coef(model.open.raw, type="det"))
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open.raw@opt$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others