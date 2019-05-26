source("./From Madsen/gen.nmix.R")
library(readr)
library("unmarked")
library("tictoc") #So we can see our code's run time

dat_ready__5_24_2019 <- read_csv("./dat_ready__5_24_2019.csv")
#Need to define these functions first
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }

data=dat_ready__5_24_2019
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

# Site Level covariates####
lon<-data[,"lon"]
lat<-data[,"lat"]
total<-data[,"TOTAL"]
desert<-data[,"DESERT"]
sltmrsh<-data[,"SLTMRSH"]
dkpnd<-data[,"DKPND"]
rcrp<-data[,"RCRP"]
grp<-data[,"GRP"]
cit<-data[,"CIT"]
datefarm<-data[,"DAT"]
pst<-data[,"PST"]
fsh<-data[,"FSH"]
trap_num<-data[,"trap_num"]
d.to.sea<-data[,"d.to.sea"]

X = cbind(lon,
          lat,
          total,
          desert,
          sltmrsh,
          dkpnd,
          rcrp,
          grp,
          cit,
          datefarm,
          pst,
          fsh,
          trap_num,
          d.to.sea)

#observation-level predictors####
temp.obs<-data[,paste(rep("temp.obs",tspan),c(1:tspan),sep="")]
temp.max<-data[,paste(rep("temp.max",tspan),c(1:tspan),sep="")]
mo<-data[,paste(rep("mo",tspan),c(1:tspan),sep="")]
wk<-data[,paste(rep("wk",tspan),c(1:tspan),sep="")]
yr<-data[,paste(rep("yr",tspan),c(1:tspan),sep="")]
season<-data[,paste(rep("season",tspan),c(1:tspan),sep="")]


DATE<-data[,paste(rep("t",tspan),c(1:tspan),sep="")]
DATE2<-DATE^2
DATE.2 = DATE

p.date = as.vector(t(DATE))
temp.obs = as.vector(t(as.matrix(temp.obs)))
temp.max = as.vector(t(as.matrix(temp.max)))
mo = as.vector(t(as.matrix(mo)))
wk = as.vector(t(as.matrix(wk)))
yr = as.vector(t(as.matrix(yr)))
season= as.vector(t(as.matrix(season)))
Z = cbind(p.date,
          temp.obs,
          temp.max,
          mo,
          wk,
          yr,
          season)

#Now set the assumed primary period length. Our data is only sampled
#every two weeks, so this should probably be small
pim.period.length <- 1
DATE.3 = ceiling(DATE.2/pim.period.length)
DATE.4 = as.matrix(DATE.3)
DATE.4[is.na(DATE.4)] = max(as.vector(DATE.3),na.rm=TRUE) 
mode(DATE.4) = "integer"

Method = "BFGS"
K.lim = max(ceiling(as.matrix(n.it)/100),na.rm=TRUE)+1
#Closed model, binned responses ####

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
K.lim = max(ceiling(as.matrix(n.it)/100),na.rm=TRUE)+1
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

lam <- 100*backTransform(model.open,type="lambda")# *100 because of binned responses
gam <- backTransform(model.open,type="gamma")
om <- backTransform(model.open,type="omega")
p<- backTransform(model.open,type="det")
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open@opt$hessian)$values
cn.open = max(ev.open)/min(ev.open)
cn.open #check condition number: this condition number is much larger than the others

#closure test

LRT(model.open,model.closed)
#Result in test stat of 2765.518

#Open Model covariates, binned responses ####

K.lim = max(ceiling(as.matrix(n.it)/100),na.rm=TRUE)+1
y = unmarkedFramePCO(ceiling(as.matrix(n.it)/100),
                     siteCovs=as.data.frame(X),
                     obsCovs=as.data.frame(Z),
                     yearlySiteCovs=NULL,
                     mapInfo=NULL,
                     numPrimary=33,
                     primaryPeriod=DATE.4)

tic()
model.open.covariate = pcountOpen(~d.to.sea,~d.to.sea,~1,~1,
                        y,
                        mixture="P",
                        K.lim,
                        dynamics="constant",
                        fix="none",
                        starts=c(model.closed@opt$par[1],rep(0,4),model.closed@opt$par[2]),
                        method=Method,
                        se=TRUE,
                        immigration=FALSE,
                        iotaformula=~1)
toc()
summary(model.open.covariate)

lam <- 100*exp(coef(model.open.covariate, type="lambda"))
gam <- exp(coef(model.open.covariate, type="gamma"))
om <- plogis(coef(model.open.covariate, type="omega"))
p <- plogis(coef(model.open.covariate, type="det"))
c(lam,gam,om,p) #back transformed fitted values

ev.open = eigen(model.open.covariate@opt$hessian)$values
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
                        K.lim.raw,
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