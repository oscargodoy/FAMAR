# first exploratory analyses FAMAR Oct 8th 
# do a correlogram
#load the database
data <- read.csv("data/famar-Biom_123_05_15c.csv", header=T, sep=",")
#we are going to limit the database to the natural growing years
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)
#put the correct format
data$S.date <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 

#We need to check for outliers

#1. leaf dry weight (L.DW)
dotchart(data$L.DW, col=data$SITE)
plot(data$S.date, data$L.DW, col=data$SITE)
identify(data$S.date, data$L.DW) #this is to kno

#This variable is ok

#2. below ground dry weight (BG.DW)
dotchart(data$BG.DW, col=data$SITE)
plot(data$S.date, data$BG.DW, col=data$SITE)
identify(data$S.date, data$BG.DW) #this is to know

#This variable is ok

#3. Rhizome dry weight  (Ri.DW)
dotchart(data$Ri.DW, col=data$SITE)
plot(data$S.date, data$Ri.DW, col=data$SITE)
identify(data$S.date, data$Ri.DW) #this is to know

#This variable is ok

#3. Root dry weight  (Ro.DW)
dotchart(data$Ro.DW, col=data$SITE)
plot(data$S.date, data$Ro.DW, col=data$SITE)
identify(data$S.date, data$Ro.DW) #this is to know

#This variable is ok

#3. total density   (dens.TOT)
dotchart(data$dens.TOT, col=data$SITE)
plot(data$S.date, data$dens.TOT, col=data$SITE)
identify(data$S.date, data$dens.TOT) #this is to know

#This variable is ok

#select the variables with max nuber of observations
data_sep <- data[,c(4,17,11,41)]

library(corrgram)
library(Hmisc)
corrgram(data_sep, order=TRUE, lower.panel=panel.conf,
         upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density,pch=16, 
         main="Correlograma variables FAMAR")


data_sep <- data[,c(4,11,13,15:24,26:36,39)]

lm1 <- lm(L.DW ~ EzSITE.mol_m2d.q2 + SEASON, data=data)
plot(data$EzSITE.mol_m2d.q2, data$L.DW, col=data$SITE)

data_sub <- subset(data, data$SEASON=="SUMMER")
lm2 <- lm(L.DW ~ EzSITE.mol_m2d.q2 + SITE, data=data_sub)