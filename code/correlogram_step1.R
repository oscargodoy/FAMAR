# first exploratory analyses FAMAR Oct 8th 
# do a correlogram
#load the database
data <- read.csv("data/famar-Biom_123_05_15d.csv", header=T, sep=",")
#put the correct format
data$S.date2 <- 
data$S.date2 <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 

#select the variables
data
data_sep <- data[,c(4,11,13,15:24,26:36,39)]

library(corrgram)
library(Hmisc)
corrgram(data_sep, order=TRUE, lower.panel=panel.conf,
         upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density,pch=16, 
         main="Correlograma variables FAMAR")

#We need to check for outliers

#1. leaf dry weight (L.DW)
dotchart(data$L.DW, col=data$SITE)


plot(data$S.date, data$L.DW, col=data$SITE)
identify(data$S.date, data$L.DW) #this is to know

data_check <- data[c(253,254,255,256,257,511,513),]
