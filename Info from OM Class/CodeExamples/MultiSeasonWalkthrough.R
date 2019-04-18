library(unmarked)
library(reshape2)

######################################
### Read in each of the 4 datasets ###
###################################### 
skink<-read.csv("MultiSeasonOccupancy_GrandSkinks.csv")
skinkObs <- skink[,c(1:5)] # but we need to reformat the data
skinkObs <- melt(skinkObs, id.vars = c("Site", "Time")) # first melt the data
skinkObs <- skinkObs[with(skinkObs, order(Time, Site, variable)),] # then re-order
names(skinkObs)<- c("Site", "Time", "Occasion", "y") # rename the columns
skinkObsFin <-dcast(skinkObs, Site ~ Time + Occasion, value.var = "y") # then transpose the data using dcast
head(skinkObsFin) # this is what we want - a row for each site, with detection/non-detection data in columns representing seasons and occasions.
dim(skinkObsFin) # we need to know this for below (how many sites there are)

# site-level covariates
skinkSiteCovs <- skink[1:338,6:8] # subset for site-level covariates; we only need the first 338 observations

# observation-level covariates: we didn't measure any occasion-specific covariates, but we can make one representing each occasion, which we'll call "Day".
skinkObsCovs <- matrix(c("Day1", "Day2", "Day3"), nrow = nrow(skink), ncol = 3, byrow = TRUE) # note that we need to specify that the matrix is filled by rows and not columns, which is the default. And we also need to make the a matrix with dimensions nSites*nYears x nOccasions (which we'll need to reformat further, below)
head(skinkObsCovs) # take a look
skinkObsCovs <- as.data.frame(skinkObsCovs) # convert to a data.frame
skinkObsCovs$Site <- rep(1:338,5) # add a site number
skinkObsCovs$Year <- rep(1:5, each = 338) # add a year number
skinkObsCovs<-melt(skinkObsCovs, id.vars = c("Site", "Year"))
skinkObsCovsFin <- skinkObsCovs[with(skinkObsCovs, order(Site, Year, variable)),] # then re-order
names(skinkObsCovsFin)<-c("Site", "Year", "Occasion", "Day")
skinkObsCovsFin<-skinkObsCovsFin[,c(1,4)] # just want to keep site and day

# We also didn't measure any year-specific covariates, but we can make one representing each year (season) in the same way as above.
skinkYearCovs <- matrix(c("01", "02", "03", "04", "05"), nrow = nrow(skinkObsFin), ncol = 5, byrow = TRUE)
skink.umf <- unmarkedMultFrame(y = skinkObsFin[,2:16], numPrimary = 5, siteCovs = skinkSiteCovs, yearlySiteCovs = list(Year=skinkYearCovs), obsCovs = skinkObsCovsFin)
(skinkmod_1 <- colext(psiformula = ~HabType, gammaformula = ~Year + HabType, epsilonformula = ~Year + HabType, pformula = ~Year, data = skink.umf))
skinkmod_2 <- colext(psiformula = ~HabType, gammaformula = ~Year + HabType, epsilonformula = ~Year + HabType, pformula = ~Year + HabType, data = skink.umf)
skinkmod_3 <- colext(psiformula = ~HabType, gammaformula = ~Year + HabType, epsilonformula = ~Year + HabType, pformula = ~Year + HabType + Day + Year*Day, data = skink.umf)
summary(skinkmod_1)
library(AICcmodavg)
aictab(list(skinkmod_1,skinkmod_2,skinkmod_3))