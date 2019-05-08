library("readxl")
library("dplyr")

filename = "C:/Users/Joe/Documents/PSU/Mosquito/9495co2t.xlsx"
cxt = read_excel(filename,skip=1)
rows = dim(cxt)[1]
cxt = cxt[1:(rows-5),]
cxt = subset(cxt,select=-c(X__1))
rows = dim(cxt)[1]

CV = function(dat){ return(sd(dat)/mean(dat)) }
YTrans = function(dat) { return(log(dat + 1)) }

cxt = cxt %>% mutate(All_Species = select(.,CXT:CSINC) %>% rowSums(na.rm=TRUE))
cxt = cxt %>% mutate(TransCXT = YTrans(CXT))
t_samps = cxt %>% select(TRAP_NUM) %>% group_by(TRAP_NUM) %>% count()
plot(t_samps$n)
cxt = inner_join(cxt, t_samps, by="TRAP_NUM")

#Table 1 data
tab1 = cxt %>%
	   summarize("Total" = sum(CXT),
	             "Mean" = mean(CXT),
	             "SE" = sd(CXT),
	             "CV" = CV(CXT))

#Figure 4
fig4a = cxt %>%
	    select(TRAP_NUM,CXT) %>%
	    group_by(TRAP_NUM) %>%
	    summarize_all(funs("Mean" = mean,
				           "Var" = var))
dev.new()
par(mfrow=c(2,1))
plot(fig4a$Mean,fig4a$Var,pch=20)
fig4b = cxt %>%
	    select(TRAP_NUM,TransCXT) %>%
	    group_by(TRAP_NUM) %>%
	    summarize_all(funs("Mean" = mean,
	                       "Var" = var))
plot(fig4b$Mean,fig4b$Var,pch=20)
