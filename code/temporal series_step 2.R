# analyses to start thinking in the GLMM temporal series
# the temporal scale of the autorregression. 
data <- read.csv("data/famar-Biom_123_05_15c.csv", header=T, sep=",")
data$S.date <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 
#we are going to limit the database to the natural growing years
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)

write.csv(data, file="data/famar_clean.csv")

#order by season than make biological sense
print(levels(data$SEASON))
data$SEASON <- factor(data$SEASON, levels(data$SEASON)[c(4,2,3,1)])

data_cn1 <- subset(data, data$SITE=="CN1")
data_cn2 <- subset(data, data$SITE=="CN2")
data_cn3 <- subset(data, data$SITE=="CN3")

#Can I use replicates to analyze partial correlation in temporal series?
#If not, average to have a single value per time.
library(plyr)

#L.DW---- 
data_cn1_ave <- ddply(data_cn1, 
                       c("S.date"),
                      summarise,
                      L.DW=quantile(L.DW, c(0.5), na.rm=T, type=8))
acfRes_cn1 <- acf(data_cn1_ave$L.DW)
pacfRes_cn1 <- pacf(data_cn1_ave$L.DW) 

data_cn2_ave <- ddply(data_cn2, 
                      c("S.date"),
                      summarise,
                      L.DW=quantile(L.DW, c(0.5), na.rm=T, type=8))
acfRes_cn2 <- acf(data_cn2_ave$L.DW)
pacfRes_cn2 <- pacf(data_cn2_ave$L.DW) 

data_cn3_ave <- ddply(data_cn3, 
                      c("S.date"),
                      summarise,
                      L.DW=quantile(L.DW, c(0.5), na.rm=T, type=8))
acfRes_cn3 <- acf(data_cn3_ave$L.DW)
pacfRes_cn3 <- pacf(data_cn3_ave$L.DW) 

#most likely we need an autocorrelation structure of COR AR 2 with L.DW

#BG.DW----
data_cn1_ave <- ddply(data_cn1, 
                      c("S.date"),
                      summarise,
                      BG.DW=quantile(BG.DW, c(0.5), na.rm=T, type=8))
acfRes_cn1 <- acf(data_cn1_ave$BG.DW)
pacfRes_cn1 <- pacf(data_cn1_ave$BG.DW) 

data_cn2_ave <- ddply(data_cn2, 
                      c("S.date"),
                      summarise,
                      BG.DW=quantile(BG.DW, c(0.5), na.rm=T, type=8))
acfRes_cn2 <- acf(data_cn2_ave$BG.DW)
pacfRes_cn2 <- pacf(data_cn2_ave$BG.DW) 

data_cn3_ave <- ddply(data_cn3, 
                      c("S.date"),
                      summarise,
                      BG.DW=quantile(BG.DW, c(0.5), na.rm=T, type=8))
acfRes_cn3 <- acf(data_cn3_ave$BG.DW)
pacfRes_cn3 <- pacf(data_cn3_ave$BG.DW)

#Partial autocorrelogram suggests that either no auto temporal cor structure of if so
#put the one of the previous time. 

#Ri.DW----
data_cn1_ave <- ddply(data_cn1, 
                      c("S.date"),
                      summarise,
                      Ri.DW=quantile(Ri.DW, c(0.5), na.rm=T, type=8))
acfRes_cn1 <- acf(data_cn1_ave$Ri.DW)
pacfRes_cn1 <- pacf(data_cn1_ave$Ri.DW) 

data_cn2_ave <- ddply(data_cn2, 
                      c("S.date"),
                      summarise,
                      Ri.DW=quantile(Ri.DW, c(0.5), na.rm=T, type=8))
acfRes_cn2 <- acf(data_cn2_ave$Ri.DW)
pacfRes_cn2 <- pacf(data_cn2_ave$Ri.DW) 

data_cn3_ave <- ddply(data_cn3, 
                      c("S.date"),
                      summarise,
                      Ri.DW=quantile(Ri.DW, c(0.5), na.rm=T, type=8))
acfRes_cn3 <- acf(data_cn3_ave$Ri.DW)
pacfRes_cn3 <- pacf(data_cn3_ave$Ri.DW)

#Partial autocorrelogram suggests that either no auto temporal cor structure of if so
#put the one of the previous time. 

#Ro.DW----
data_cn1_ave <- ddply(data_cn1, 
                      c("S.date"),
                      summarise,
                      Ro.DW=quantile(Ro.DW, c(0.5), na.rm=T, type=8))
acfRes_cn1 <- acf(data_cn1_ave$Ro.DW)
pacfRes_cn1 <- pacf(data_cn1_ave$Ro.DW) 

data_cn2_ave <- ddply(data_cn2, 
                      c("S.date"),
                      summarise,
                      Ro.DW=quantile(Ro.DW, c(0.5), na.rm=T, type=8))
acfRes_cn2 <- acf(data_cn2_ave$Ro.DW)
pacfRes_cn2 <- pacf(data_cn2_ave$Ro.DW) 

data_cn3_ave <- ddply(data_cn3, 
                      c("S.date"),
                      summarise,
                      Ro.DW=quantile(Ro.DW, c(0.5), na.rm=T, type=8))
acfRes_cn3 <- acf(data_cn3_ave$Ro.DW)
pacfRes_cn3 <- pacf(data_cn3_ave$Ro.DW)
#Partial autocorrelogram suggests that either no auto temporal cor structure of if so
#put the one of the previous time. 

#dens.TOT----

data_cn1_ave <- ddply(data_cn1, 
                      c("S.date"),
                      summarise,
                      dens.TOT=quantile(dens.TOT, c(0.5), na.rm=T, type=8))
acfRes_cn1 <- acf(data_cn1_ave$dens.TOT)
pacfRes_cn1 <- pacf(data_cn1_ave$dens.TOT) 

data_cn2_ave <- ddply(data_cn2, 
                      c("S.date"),
                      summarise,
                      dens.TOT=quantile(dens.TOT, c(0.5), na.rm=T, type=8))
acfRes_cn2 <- acf(data_cn2_ave$dens.TOT)
pacfRes_cn2 <- pacf(data_cn2_ave$dens.TOT) 

data_cn3_ave <- ddply(data_cn3, 
                      c("S.date"),
                      summarise,
                      dens.TOT=quantile(dens.TOT, c(0.5), na.rm=T, type=8))
acfRes_cn3 <- acf(data_cn3_ave$dens.TOT)
pacfRes_cn3 <- pacf(data_cn3_ave$dens.TOT)
#Partial autocorrelogram suggests that either no auto temporal cor structure of if so
#put the one of the previous time or two times ago. 
