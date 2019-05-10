# To Do:
#   Add NOAA Weather data


#Import data [Note: Use data set that includes the "T" column!] ####
library("readxl")
file_in_joe = "mosDat.xlsx"
file_in_jac = "C:/Users/Jacob/Dropbox/Grad School/2018-2019/Spring/Consulting Class/mosDat.xlsx"
mosDat <- read_excel(file_in_jac, 
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

mosDat<-cbind(mosDat,temps.mosDat)
#Max number of revisits?####
m=0
for(i in c(1:64)){
  if (length(mosDat[mosDat$TRAP_NUM == i,1]$TRAP_NUM) > m){
    m <-i
  }
}

m #max number of revisits is 33

#Reformat time + count data####

#First build empty data frame to contain the counts, temps, and times
num.traps<-64
num.vis<-33

ys<-paste(rep("y",num.vis), c(1:num.vis), sep="")
ts<-paste(rep("t",num.vis), c(1:num.vis),sep="")
temp.maxs<-paste(rep("temp.max",num.vis),c(1:num.vis),sep="")
temp.obss<-paste(rep("temp.obs",num.vis),c(1:num.vis),sep="")
dat_ready<-data.frame(matrix(ncol = 133, nrow = num.traps))
colnames(dat_ready)<-c("trap_num",ts,ys,temp.maxs,temp.obss)
dat_ready[,1]<-c(1:num.traps)
dat_ready
dim(dat_ready)

# Loop formats same as mallard.data
for (i in c(1:num.traps)){
  d<-mosDat[mosDat$TRAP_NUM == i,] #Grab all data assoc with trap i
  
  t<-d$T
  y<-d$CXT
  temp.m<-d$TMAX
  temp.o<-d$TOBS
  
  length(t)<-num.vis #These lines add "NA" to time, count vec's if needed
  length(y)<-num.vis
  length(temp.m)<-num.vis
  length(temp.o)<-num.vis
  
  dat_ready[i,c(2:133)]<-t(c(t,y,temp.m,temp.o))
}
dat_ready

# Add site locations
library(readr)
TrapSites <- read_csv("TrapSites.csv")
dat_ready <- merge(TrapSites,dat_ready,all=TRUE,by.x="Name",by.y="trap_num")
colnames(dat_ready)[colnames(dat_ready)=="Name"] <- "trap_num"
colnames(dat_ready)[colnames(dat_ready)=="X"] <- "lon"
colnames(dat_ready)[colnames(dat_ready)=="Y"] <- "lat"
dat_ready$description <- NULL

# Add habitat data
habitats <- read_excel("Habitats from Hugh L.xls")
habitats$X__1 <- NULL
habitats$X__2 <- NULL
habitats <- habitats[-c(63:72),]
habitats[63,] <- list("0064",156.46,0,0,0,0,0,0,0,0,156.46)
habitats$TRAP_NUM <- as.numeric(habitats$TRAP_NUM)
dat_ready <- merge(habitats,dat_ready,all=TRUE,by.x="TRAP_NUM",by.y="trap_num")
colnames(dat_ready)[colnames(dat_ready)=="TRAP_NUM"] <- "trap_num"

# Write to CSV using R####
<<<<<<< HEAD
write.csv(dat_ready, file = "D:\\Desktop\\dat_ready.csv",row.names=FALSE, na="")

=======
file_out_joe = "dat_ready.csv"
file_out_jac = "D:\\Desktop\\dat_ready.csv"
write.csv(dat_ready, file = file_out_jac,row.names=FALSE, na="")
>>>>>>> a228b68959fc0d388af40a43c22bff2d776aa92b
