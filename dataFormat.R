# To Do:
#   Add on site locations
#   Add NOAA Weather data
#   Add Habitat proportions


#Import data [Note: Use data set that includes the "T" column!] ####
mosDat <- read_excel("C:/Users/Jacob/Dropbox/Grad School/2018-2019/Spring/Consulting Class/mosDat.xlsx", 
                     col_types = c("numeric","numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric"))

# Dealing with Temperatures####
library(readxl)
NOAAsubset <- read_excel("C:/Users/Jacob/Dropbox/Grad School/2018-2019/Spring/Consulting Class/Data/NOAAsubset.xlsx")
dim(NOAAsubset)
#View(NOAAsubset)

#first find the two week averages, which is based on the time period "T"
temps<-aggregate(NOAAsubset[, 2:3], list(NOAAsubset$T), mean, na.rm=TRUE)
temps

#rename the time column to match the notation in mosDat
colnames(temps)[colnames(temps)=="Group.1"] <- "T"
temps

temps.mosDat<- data.frame(matrix(ncol = 2, nrow = length(mosDat$T)))
colnames(temps.mosDat) <- c("TMAX", "TOBS")
for (i in c(1:length(mosDat$T))){
  temps.mosDat[i,]<-temps[which(temps$T==mosDat$T[i]),c(2,3)]
}

#These are the temperatures for each observation.
#This is meant to be added onto mosDat
temps.mosDat

#Max number of revisits?####
m=0
for(i in c(1:64)){
  if (length(mosDat[mosDat$TRAP_NUM == i,1]$TRAP_NUM) > m){
    m <-i
  }
}

m #max number of revisits is 33

#Reformat time + count data####

#First build empty data frame to contain the counts and times
num.traps<-64
num.vis<-33
ys<-rep("y",num.vis)
ys<-paste(ys, c(1:num.vis), sep="")
ts<-paste(rep("t",num.vis),c(1:num.vis),sep="")
dat_ready<-data.frame(matrix(ncol = 67, nrow = num.traps))
colnames(dat_ready)<-c("trap_num",ts,ys)
dat_ready[,1]<-c(1:num.traps)
dat_ready


#Loop formats same as mallard.data
for (i in c(1:num.traps)){
  d<-mosDat[mosDat$TRAP_NUM == i,c(1:6)] #Grab all data assoc with trap i
  t<-d$T
  y<-d$CXT
  length(t)<-num.vis #These lines add "NA" to time, count vec's if needed
  length(y)<-num.vis
  dat_ready[i,c(2:67)]<-t(c(t,y))
}
dat_ready

#Add site locations####
library(readr)
TrapSites <- read_csv("TrapSites.csv")


# Write to CSV using R####
write.csv(dat_ready, file = "D:\\Desktop\\dat_ready.csv",row.names=FALSE, na="")

