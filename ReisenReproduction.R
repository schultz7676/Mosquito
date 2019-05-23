library("readxl")
library("fBasics")
library("dplyr")
#setwd("C:/Users/Joe/Documents/PSU/Mosquito")

# Read the data from Microsoft Excel file, clean up extraneous data
filename = "C:/Users/Joe/Documents/PSU/Mosquito/9495co2t.xlsx"
cxt = read_excel(filename,skip=1)
rows = dim(cxt)[1]
cxt = cxt[1:(rows-5),]
cxt = subset(cxt,select=-c(X__1))
rows = dim(cxt)[1]

# Define some functions
source("ReisenHelperFuncs.R")

# Add column of grand total counts including all species
cxt = cxt %>% mutate(All_Species = select(.,CXT:CSINC) %>% rowSums(na.rm=TRUE))

# Add site locations
library(readr)
TrapSites = read_csv("TrapSites.csv")
TrapSites$description <- NULL
colnames(TrapSites)[colnames(TrapSites)=="X"] <- "LON"
colnames(TrapSites)[colnames(TrapSites)=="Y"] <- "LAT"
TrapSites = TrapSites %>% mutate(Name = tochr(Name))
cxt = cxt %>% inner_join(TrapSites,by=c("TRAP_NUM" = "Name"))

# Add habitat data
habitats = read_excel("Habitats from Hugh L.xls")
habitats$X__1 <- NULL
habitats$X__2 <- NULL
habitats <- habitats[-c(63:72),]
habitats[63,] <- list("0064",156.46,0,0,0,0,0,0,0,0,156.46)
habitats = habitats %>% mutate(TRAP_NUM = truncstr(TRAP_NUM))
cxt = cxt %>% inner_join(habitats,by="TRAP_NUM")

# Add temperature data
noaa = read_excel("NOAAsubset.xlsx")
noaa = noaa %>% mutate(T2 = date_to_T(DATE))
noaaavg = noaa %>%
          select(T,TMAX,TOBS) %>%
		  group_by(T) %>%
		  summarize("TMAX" = mean(TMAX,na.rm=TRUE),
					"TOBS" = mean(TOBS,na.rm=TRUE))
cxt = cxt %>% mutate(T = yrwk_to_T(YR,WK))
cxt = cxt %>% inner_join(noaaavg,by="T")

# Add distance to Salton Sea
seabound = read_csv("CountOverlays/saltonseaboundary.csv")
cxt = cxt %>% rowwise() %>% mutate(DIST_TO_SEA = dist_to_sea(LAT,LON,seabound))

# Add column of transformed data
cxt = cxt %>% mutate(TransCXT = YTrans(CXT))

# Count the number of time samples at each trap site and add a column for it
t_samps = cxt %>% select(TRAP_NUM) %>% group_by(TRAP_NUM) %>% count()
plot(t_samps$n)
cxt = inner_join(cxt, t_samps, by="TRAP_NUM")
colnames(cxt)[colnames(cxt)=="n"] = "t_samps"

# Table 1 data
tab1 = cxt %>%
	   summarize("Total" = sum(CXT),
	             "Mean" = mean(CXT),
	             "SE" = sd(CXT),
	             "CV" = CV(CXT),
	             "TransCV" = CV(TransCXT))
tab1

# Figure 4
fig4a = cxt %>%
		filter(t_samps >= 31) %>%
	    select(TRAP_NUM,CXT) %>%
	    group_by(TRAP_NUM) %>%
	    summarize_all(funs("Mean" = mean,
				           "Var" = var))
lm4a = lm(fig4a$Var ~ fig4a$Mean)
summary(lm4a)
fig4b = cxt %>%
		#filter(t_samps >= 30) %>%
	    select(TRAP_NUM,TransCXT) %>%
	    group_by(TRAP_NUM) %>%
	    summarize_all(funs("Mean" = mean,
				           "Var" = var))
lm4b = lm(fig4b$Var ~ fig4b$Mean)
summary(lm4b)
dev.new()
par(mfrow=c(2,1))
plot(fig4a$Mean,fig4a$Var,pch=20)
plot(fig4b$Mean,fig4b$Var,pch=20)

# Confidence Intervals for various estimates of central tendency
# untransformed data
ndat = length(cxt$CXT)
df = ndat-1
alpha = 0.05
lcl = tab1$Mean + tab1$SE/sqrt(ndat)*qt(alpha/2,df)
ucl = tab1$Mean - tab1$SE/sqrt(ndat)*qt(alpha/2,df)
c(lcl,tab1$Mean,ucl)
# backtransformed data
cent = mean(cxt$TransCXT)
lcl = cent + sd(cxt$TransCXT)/sqrt(ndat)*qt(alpha/2,df)
ucl = cent - sd(cxt$TransCXT)/sqrt(ndat)*qt(alpha/2,df)
YInvTrans(c(lcl,cent,ucl))
# geometric mean
gm_mean(cxt$CXT)

# Figure 5
brksa = c(0,727,1636,2545,3455,4364,5273,6182,7091,7937)
brksb = -1:10
dev.new()
par(mfrow=c(2,1))
fig5a = hist(cxt$CXT,breaks=brksa,freq=TRUE)
fig5b = hist(cxt$TransCXT,breaks=brksb,freq=TRUE)
skewness(cxt$CXT,method="moment")[1]
kurtosis(cxt$CXT,method="moment")[1]
skewness(cxt$TransCXT,method="moment")[1]
kurtosis(cxt$TransCXT,method="moment")[1]

# Figure 6
fig6 = cxt %>%
	   select(TRAP_NUM,CXT) %>%
	   group_by(TRAP_NUM) %>%
	   summarize_all(funs("GeoMean" = gm_mean))
geomeans = sort(fig6$GeoMean,decreasing=TRUE)
dev.new()
barplot(geomeans,names.arg=1:63)

# ANOVA
library("car")
lmA = lm(TransCXT ~ factor(TRAP_NUM) + factor(MO) + factor(YR),data=cxt)
Anova(lmA,type=2)
aovall = aov(TransCXT~factor(TRAP_NUM),data=cxt)
library("agricolae")
snkall = SNK.test(aovall,"factor(TRAP_NUM)",console=FALSE)
head(snkall$groups)

# Figure 7
library("tidyverse")
month_names = c("F","M","A","M","J","J","A","S","O","N")
month_trans = tibble(MO=formatC(2:11, width = 2, format = "d", flag = "0"),LABEL=month_names)
fig7a = cxt %>%
		filter(YR == "1994") %>%
	    select(MO,CXT) %>%
	    group_by(MO) %>%
	    summarize_all(funs("GeoMean" = gm_mean))
fig7a[nrow(fig7a) + 1,] = list("02",NaN)
fig7a[nrow(fig7a) + 1,] = list("03",NaN)
fig7a[nrow(fig7a) + 1,] = list("11",NaN)
fig7a = fig7a[order(fig7a$MO),]
fig7a = add_column(fig7a,YR="1994")
fig7b = cxt %>%
		filter(YR == "1995") %>%
	    select(MO,CXT) %>%
	    group_by(MO) %>%
	    summarize_all(funs("GeoMean" = gm_mean))
fig7b = add_column(fig7b,YR="1995")
fig7c = cxt %>%
		filter((YR == "1995") & (as.numeric(WK) %in% seq(9,45,by=4))) %>%
	    select(MO,CXT) %>%
	    group_by(MO) %>%
	    summarize_all(funs("GeoMean" = gm_mean))
fig7c = add_column(fig7c,YR="1995/A")
fig7d = cxt %>%
	    select(MO,CXT) %>%
	    group_by(MO) %>%
	    summarize_all(funs("GeoMean" = gm_mean))
fig7d = add_column(fig7d,YR="All")
fig7 = bind_rows(fig7a,fig7b,fig7c,fig7d)
fig7 = inner_join(fig7, month_trans, by="MO")
dev.new()
ggplot(data=fig7, mapping=aes(x=MO, y=GeoMean, group=YR)) + geom_line()

# Table 2
cxt$SEASON = apply(cxt,1,FUN=season)
tab2 = cxt %>%
       select(SEASON,CXT,TransCXT) %>%
	   group_by(SEASON) %>%
	   summarize("No. Trap Nights" = n(),
				 "Mean" = mean(TransCXT),
				 "Geometric Mean" = gm_mean2(CXT),
				 "CL" = CL(TransCXT),
				 "CV" = CV(TransCXT))
tab2 %>% gather(var,val,2:ncol(tab2)) %>% spread(SEASON,val)
spr = cxt %>% filter(SEASON=="spring") %>% select(TRAP_NUM,TransCXT)
aovspr = aov(TransCXT~factor(TRAP_NUM),data=spr)
snkspr = SNK.test(aovspr,"factor(TRAP_NUM)", console=FALSE)
head(snkspr$groups)
sum = cxt %>% filter(SEASON=="summer") %>% select(TRAP_NUM,TransCXT)
aovsum = aov(TransCXT~factor(TRAP_NUM),data=sum)
snksum = SNK.test(aovsum,"factor(TRAP_NUM)", console=FALSE)
head(snksum$groups,7)
fal = cxt %>% filter(SEASON=="fall") %>% select(TRAP_NUM,TransCXT)
aovfal = aov(TransCXT~factor(TRAP_NUM),data=fal)
snkfal = SNK.test(aovfal,"factor(TRAP_NUM)", console=FALSE)
head(snkfal$groups,7)

# Exploratory data analysis
lmB = lm(TransCXT ~ factor(MO) +
		            factor(YR) +
					TOBS +
					I(SLTMRSH/TOTAL) +
					I(DKPND/TOTAL) +
					I(RCRP/TOTAL) +
					I(GRP/TOTAL) +
					I(CIT/TOTAL) +
					I(DAT/TOTAL) +
					I(PST/TOTAL) +
					I(FSH/TOTAL) +
					TOTAL +
					LAT +
					LON,
		 data=cxt)
summary(lmB)
require("car")
print(vif(lmB))	# Calculates the variance inflation factors

lmC=lm(TransCXT~factor(TRAP_NUM)+factor(MO)+factor(YR)+TOBS,data=cxt)
summary(lmC)
print(vif(lmC))

# Notes on multicollinearity:
#	TMAX conflicts with factor(MO) and TOBS
#	factor(TRAP_NUM) conflicts with LAT, LON, and the habitat covariates
#	including all the habitat types as raw covariates causes problems, best to leave out DESERT

lmD=lm(I(CXT+1)~factor(TRAP_NUM)+factor(MO)+factor(YR)+TOBS,data=cxt)
library("MASS")
lambda_curve = boxcox(lmD)
opt_power = lambda_curve$x[lambda_curve$y==max(lambda_curve$y)]
opt_power

# Notes on optimal transform:
#	Box-Cox shows that the transform used in Reisen is near optimal

dev.new()
par(mfrow=c(2,2))
#1. Residuals vs. predicted values
yhat = fitted(lmC)
lmC_resid = resid(lmC)
plot(yhat,lmC_resid,xlab="Predicted Values",ylab="Residuals")
#2. Histogram or box plot of residuals
boxplot(lmC_resid)
title("Residual Boxplot")
#3. Q-Q plot (need car package)
qqPlot(lmC)
title("Model C")
#6. Serial autocorrelations (need forecast package)
library("forecast")
dev.new()
Acf(residuals(lmC))
#7. Durbin-Watson test (need lmtest package)
library("lmtest")
dwtest(lmC)

#PRESS statistic for both models (need plyr package and source the .r files)
#model_fit_stats(lmC)

# Use car package
# Influential Observations
# Cook's D plot
# identify D values > 4/(n-k-1) 
dev.new()
cutoff = 4/(nrow(cxt)-length(lmC$coefficients)-2) 
plot(lmC_resid, which=4, cook.levels=cutoff)
influencePlot(lmC,
			  id.method="noteworthy",
			  id.n=4,
			  main="Influence Plot",
			  sub="Circle size is proportial to Cook's D")
