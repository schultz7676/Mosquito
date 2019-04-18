

library(reshape2)
nSites <- 100  # number of site
nReps <- 3  # number of surveys
occasion <- rep(1:nReps, nSites)  # survey number
site <- rep(1:nSites, each = 3)  # site number
cov1 <- rep(runif(nSites, 0, 10), each = 3)  # covariate, constant among surveys
cov2 <- runif(length(site), 0, 100)  # covariate, varies among surveys
y <- rbinom(length(site), 1, 0.4)  # survey data with p = 0.4
dater <- data.frame(site, occasion, y, cov1, cov2)  # resulting data
head(dater)


dethist <- dcast(dater, site + cov1 ~ occasion, value.var = "y")
head(dethist)


covhist <- dcast(dater, site ~ occasion, value.var = "cov2")
head(covhist)


Day <- matrix(c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5"), nrow = 200, ncol = 5, byrow = TRUE)
