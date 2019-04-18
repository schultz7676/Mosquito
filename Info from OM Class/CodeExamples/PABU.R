### PAINTED BUNTING MULTI-SEASON DATA
library(unmarked)
PaBu<-read.csv("MultiSeasonOccupancy_PaBu.csv")
PaBuObs <- PaBu[,c(2:40)] # but we don't need to reformat the data this time

# site-level covariates
PaBuSiteCovs <- PaBu[,41:44] # subset for site-level covariates; we only need the first 338 observations
dim(PaBuObs)


PaBu.umf <- unmarkedMultFrame(y = PaBuObs, numPrimary = 13, siteCovs = PaBuSiteCovs) # there are no occasion-specific covariates


PaBumod_1 <- colext(psiformula = ~Lat + MeanLHeight, gammaformula = ~Lat+ MeanLHeight, epsilonformula = ~Lat, pformula = ~Shrub + SaltMarsh, data = PaBu.umf)
PaBumod_2 <-colext(psiformula= ~Lat, gammaformula= ~Lat, epsilonformula= ~ Lat, pformula= ~Shrub + SaltMarsh, data=PaBu.umf)
PaBumod_3<-colext(psiformula= ~Lat+ Shrub + SaltMarsh, gammaformula= ~Lat+Shrub, epsilonformula = ~Lat+Shrub, pformula = ~MeanLHeight, data=PaBu.umf)

PaBu_0<-colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1, pformula = ~1, data = PaBu.umf)

summary(PaBumod_1)
summary(PaBumod_2)
summary(PaBumod_3)

models<- fitList(m0=PaBu_0, m1=PaBumod_1, m2=PaBumod_2, m3=PaBumod_3)
modSel(models)

##compute observed chi-square
obs <- mb.chisq(PaBumod_3)
obs
##round to 4 digits after decimal point
print(obs, digits.vals = 4)
##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- mb.gof.test(PaBumod_3, nsim = 100)
obs.boot
print(obs.boot, digits.vals = 4, digits.chisq = 4)

