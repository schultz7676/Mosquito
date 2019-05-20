library(readxl)
NOAAsubset <- read_excel("NOAAsubset.xlsx") #added dataset to git 5/20/19

timeperiod <- NOAAsubset$T
mtemps <- NOAAsubset$TMAX
otemps <- NOAAsubset$TOBS
day <- c(1:length(timeperiod))

#let's look at the scatterplots first
plot(timeperiod,mtemps)
plot(day,mtemps)

#fit max temps over all 42 time periods (sinusoidal)
timeperiod.c<-cos(2*pi*timeperiod/(365.25/14))
timeperiod.s<-sin(2*pi*timeperiod/(365.25/14))
fit.lm<-lm(mtemps~timeperiod.c+timeperiod.s)
summary(fit.lm)

pred <- predict(fit.lm, newdata=data.frame(timeperiod=timeperiod))    

plot(mtemps~timeperiod, data=NOAAsubset, xlim=c(0, 43))
lines(timeperiod,pred,col="blue")


#fit max temps over all 588 days (sinusoidal)
timeperiod.c<-cos(2*pi*day/(365.25))
timeperiod.s<-sin(2*pi*day/(365.25))
fit.lm<-lm(mtemps~timeperiod.c+timeperiod.s)
summary(fit.lm)

pred <- predict(fit.lm, newdata=data.frame(day=day))    

plot(mtemps~day, data=NOAAsubset, xlim=c(0, 588))
lines(day,pred,col="blue")
