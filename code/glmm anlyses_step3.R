#load the data and start with the analyses
library(reshape2)
library(nlme)

data <- read.csv("data/famar-Biom_123_05_15.csv", header=T, sep=",")
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)
#put the correct format
data$S.date <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 

#create a new variable to assign a specific plot
data$PLOT <- paste(data$SITE, data$R.replicates, sep="_") 

#We are going to establish a linear relationship between biomass
#production and the predictors by doing a log of biomass
#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

#1.leaf dry weight (L.DW)----
model_ldw1 <- lme(log(L.DW) ~ SITE, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
#According to AIC an autocorrelation structure with CorARMA (p=2, q=3) is the best (lowest AIC)
#Now I am adding another potential predictors. 
model_ldw2 <- lme(log(L.DW) ~ SITE + NH4.uM.mean, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

model_ldw3 <- lme(log(L.DW) ~ SITE + NH4.uM.q2 + NO3.uM.q2 + PO4.uM.q2, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
# This last model is the one that work best 
AIC(model_ldw1, model_ldw2, model_ldw3)

summary(model_ldw3)

#Model checking plot 
residuals <- resid(model_ldw3)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_ldw3)
qqnorm(model_ldw3, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)


#2.below ground dry weight (BG.DW)----
#built from the experience of the previous model with L.DW
model_bgdw1 <- lme(log(BG.DW) ~ SITE + NO3.uM.q2, data=data, random=~1|PLOT, control=lCtr, correlation= corAR1(),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
AIC(model_bgdw1)
summary(model_bgdw1)

#Model checking plot 
residuals <- resid(model_bgdw1)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_bgdw1)
qqnorm(model_bgdw1, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#3.rhizome dry weight (Ri.DW)----
model_ridw1 <- lme(log(Ri.DW) ~ SITE + NO3.uM.q2, data=data, random=~1|PLOT, control=lCtr, correlation= corAR1(),
                   weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
AIC(model_ridw1)
summary(model_ridw1)

#Model checking plot 
residuals <- resid(model_ridw1)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_ridw1)
qqnorm(model_ridw1, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#4.rhizome dry weight (Ro.DW)----
data2<-subset(data, data$Ro.DW>0) #doing so to eliminate a row with Ro.Dw value equal to zero.
model_rodw1 <- lme(log(Ro.DW) ~ SITE + NO2.uM.q2 + NO3.uM.q2, data=data2, random=~1|PLOT, control=lCtr, correlation= corAR1(),
                   weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
AIC(model_rodw1)
summary(model_rodw1)

#Model checking plot 
residuals <- resid(model_rodw1)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_rodw1)
qqnorm(model_rodw1, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#5.total density (dens.TOT)----
model_densTOT1 <- lme(log(dens.TOT) ~ SITE + NH4.uM.q2 + NO2.uM.q2 + NO3.uM.q2, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | YEAR), method='ML',na.action=na.omit)
#sometimes is better to modelate the weights of year rather than the season to control for heterocedasticity
AIC(model_densTOT1)
summary(model_densTOT1)

#Model checking plot 
residuals <- resid(model_densTOT1)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_densTOT1)
qqnorm(model_densTOT1, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)


##Final analyses ----

#load the data and start with the analyses
library(reshape2)
library(nlme)

#data <- read.csv("data/famar-Biom_123_05_15_nut.c.csv", header=T, sep=",")
data <- read.csv("data/FAMAR_punching123_05_15_nut.c.csv", header=T, sep=",")
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)
#put the correct format
data$S.date <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 

lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

model_ldw3 <- lme(log(L.DW) ~ SITE + NH4.uM.q2 + NO3.uM.q2 + PO4.uM.q2, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
# This last model is the one that work best 
AIC(model_ldw1, model_ldw2, model_ldw3)