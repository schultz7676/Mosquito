library(unmarked) 
library(dplyr) # need this for data manipulation
library(tidyr) # need this for data manipulation
library(reshape2)
options(scipen=999) # suppress printing scientific notation
walkeri<-read.csv("SingleSeasonOccupancy_Mwalkeri.csv")
walkeriObs<-walkeri[,17:24] # subset columns for observation data
walkeriObsCovs<-walkeri[,c(1,25:40)] # subset columns for observation-level covariates
walkeriSiteCovs<-walkeri[,1:16] # subset columns for site-level covariates
Effort<-walkeriObsCovs[,1:9]
stEffort<-walkeriObsCovs[,c(1,10:17)]
# melt and sort data using rehape2
Effort2 <- melt(Effort, id.vars = "Site")
Effort2 <- Effort2[with(Effort2, order(Site, variable)),]
stEffort2 <- melt(stEffort, id.vars = "Site")
stEffort2 <- stEffort2[with(stEffort2, order(Site, variable)),]
EffortFin2 <- data.frame(Site = Effort2[,1], Effort = Effort2[,3], stEffort = stEffort2[,3])
print(EffortFin2[1:24,])

# melt and sort data using dplyr/tidyr
Effort3 <- Effort %>% gather(variable, value, -Site) %>% arrange(Site, variable)
stEffort3 <- stEffort %>% gather(variable, value, -Site) %>% arrange(Site, variable)
EffortFin3 <- data.frame(Site = Effort3$Site, Occasion = Effort3$variable, Effort = Effort3$value, stEffort = stEffort3$value)
print(EffortFin3[1:24,])
# create unmarkedFrame - note that if we supply unmarked with a melted (long format) obscovs data frame, we can read it in directly...
# I think another benefit to entering the obscovs this way is that you've got a single dataframe that contains all of the observation-
# level covariates you want to use
walkeri.umf1 <- unmarkedFrameOccu(y = walkeriObs, siteCovs = walkeriSiteCovs, obsCovs = EffortFin3)

# but if we supply unmarked with a Site x Occasion dataframe (or matrix), then we have to enter it as so....
walkeri.umf2 <- unmarkedFrameOccu(y = walkeriObs, siteCovs = walkeriSiteCovs, obsCovs = list(stEffort = stEffort[,2:9]))
Effort
(walk_mod1 <- occu(~stEffort ~ stSpringQcms + I(stSpringQcms^2), walkeri.umf1))
(walk_mod2 <- occu(~stEffort ~ stSpringQcms + I(stSpringQcms^2), walkeri.umf2))