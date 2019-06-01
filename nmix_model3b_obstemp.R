source("./From Madsen/gen.nmix.R")
library(readr)
library("unmarked")
library("tictoc") #So we can see our code's run time

source("./nmix_model_setup.R")

Method = "BFGS"
K.lim = max(ceiling(as.matrix(n.it)/100),na.rm=TRUE)+20

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

#Open Model covariates, binned responses ####
y = unmarkedFramePCO(ceiling(as.matrix(n.it)/100),
                     siteCovs=as.data.frame(X),
                     obsCovs=as.data.frame(Z),
                     yearlySiteCovs=NULL,
                     mapInfo=NULL,
                     numPrimary=33,
                     primaryPeriod=DATE.4)

tic()
model.open.covariate = pcountOpen(~1,~1,~1,~temp.obs,
                        y,
                        mixture="P",
                        K.lim,
                        dynamics="constant",
                        fix="none",
                        starts=c(model.closed@opt$par[1],rep(0,34),model.closed@opt$par[2]),
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
cn.open
