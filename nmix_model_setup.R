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

