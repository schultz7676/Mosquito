library(readxl)
library(corrplot)

#First inport data. If dealing with the raw data make sure the first row is deleted, as well as the misc rows that come after the main data set (so delete row 1 and anything after row 1973)
mosDat <- read_excel("C:/Users/Jacob/Dropbox/Grad School/2018-2019/Spring/Consulting Class/mosDat.xlsx", 
                     col_types = c("numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric", "numeric", 
                                   "numeric", "numeric"))

summary(mosDat)
#Note: zero detections of CSINC, one NA in PC.
# This command helps us find NA in PC:
which(is.na(mosDat$PC))
#it is row 638. As a temp workaround, assign this as zero so we can look at correlations.
#(note the other columns in this row are also zero. Seems like a reasonable patch for now)
mosDat[638,13]<-0


#define some subset of the data and build a cor matrix
x<-mosDat[,c(5:17)] #entire data set
x<-subset(mosDat,WK==16)[,c(5:17)] #looks at just the 16th week. Adjust subset condition as needed

x<-subset(mosDat,MO==4)[,c(5:17)]               #Data from April 1994 and 1995
x<-subset(mosDat,MO==4 & YR==1994)[,c(5:17)]    #Data from April 1994


cor(x)
corrplot(cor(x))
