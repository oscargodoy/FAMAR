#load the data and start with the analyses
library(reshape2)
library(nlme)
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(DHARMa)

data <- read.csv("data/famar-Biom_123_05_15_nut.c.csv", header=T, sep=",")
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)
#put the correct format
data$S.date <- as.POSIXct(strptime(data$S.date, format="%d/%m/%Y")) 

#create a new variable to assign a specific plot
data$PLOT <- paste(data$SITE, data$R.replicates, sep="_") 

data$fELEVATION <- as.factor(data$ELEVATION) 

#We are going to establish a linear relationship between biomass
#production and the predictors by doing a log of biomass
#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

#1.leaf dry weight (L.DW)----
model_ldw1 <- lme(log(L.DW) ~ fELEVATION, data=data, random=~1|R.replicates, control=lCtr, correlation= corARMA(p=2, q=3), 
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
#According to AIC an autocorrelation structure with CorARMA (p=2, q=3) is the best (lowest AIC)
#Now I am adding another potential predictors. 
model_ldw2 <- lme(log(L.DW) ~ fELEVATION + NH4.uM.q2, data=data, random=~1|R.replicates, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

model_ldw3 <- lme(log(L.DW) ~ YEAR*fELEVATION + EzSITE.mol_m2d.q2 + NO3.uM.q2 + NH4.uM.q2, data=data, random=~1|R.replicates/ELEVATION, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 


model_ldw4 <- lm(log(L.DW) ~ YEAR*fELEVATION + EzSITE.mol_m2d.q2 + NO3.uM.q2 + NH4.uM.q2, data=data,na.action=na.omit) 

model_ldw5 <- lme(log(BG.DW) ~ YEAR*fELEVATION + EzSITE.mol_m2d.q2 + NO3.uM.q2, data=data, random=~1|R.replicates/ELEVATION, control=lCtr, correlation= corARMA(p=1, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
# This last model is the one that work best 
AIC(model_ldw1, model_ldw2, model_ldw3, model_ldw4)

summary(model_ldw5)
plot(model_ldw3) # check that residuals are ok 


#plot some of the data according to the best fitted model

#we are going to plot ldw as a function of nitrogen of the three elevations. 

p1<-ggplot(data, aes(x = YEAR, y = NH4.uM.q2, color = fELEVATION) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15, aes(fill = fELEVATION))

p2<-ggplot(data, aes(x = YEAR, y = log(L.DW), color = fELEVATION) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15, aes(fill = fELEVATION))

grid.arrange(p1,p2, ncol=2)



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
model_densTOT1 <- lme(log(dens.TOT) ~ YEAR*fELEVATION + NO3.uM.q2, data=data, random=~1|R.replicates/fELEVATION, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | YEAR), method='ML',na.action=na.omit)
#sometimes is better to modelate the weights of year rather than the season to control for heterocedasticity
AIC(model_densTOT1)
summary(model_densTOT1)

ggplot(data, aes(x = YEAR, y = dens.TOT, color = fELEVATION) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15, aes(fill = fELEVATION))

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

bio <- read.csv("data/famar-Biom_123_05_15_minVAR.csv", header=T, sep=",")
bio <- subset(bio, bio$YEAR>2005)
bio <- subset(bio, bio$YEAR<2015)
table(bio$YEAR)
table(bio$SITE)
#put the correct format
bio$S.date <- as.POSIXct(strptime(bio$S.date, format="%d/%m/%Y")) 

lCtr <- lmeControl(maxIter = 5000, msMaxIter = 5000, tolerance = 1e-8, niterEM = 2500, msMaxEval = 2000)

#5.1 LDW----
model_ldw <- lme(log(L.DW) ~ SITE*YEAR*NH4.uM.mean + NO3.uM.mean, data=bio, random=~1|R.replicates/SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
# NOTE: SEASON NEEDS to be in the following order WINTER, SPRING, SUMMER, FALL. 
# This last model is the one that work best 
summary(model_ldw)
anova(model_ldw)
####ES MUY IMPORTANTE VER LA RELACION ENTRE SITIO Y A˜NO Y LUEGO LO SUYO CON NUTRIENTES.

#Model checking plot 
residuals <- resid(model_ldw)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_ldw)
qqnorm(model_ldw, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#5.2 BGDW----
model_bgdw <- lme(log(BG.DW) ~ ELEVATION + YEAR + NO3.uM.mean, data=bio, random=~1|R.replicates/SITE, control=lCtr, correlation= corAR1(),
                   weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 

# This last model is the one that work best 
summary(model_bgdw)


#Model checking plot 
residuals <- resid(model_bgdw)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_bgdw)
qqnorm(model_bgdw, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#5.3 dens.TOT----
model_densTOT <- lme(log(dens.TOT) ~ ELEVATION*NO3.uM.mean + NO2.uM.mean + NH4.uM.mean, data=bio, random=~1|R.replicates/SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                     weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 

# This last model is the one that work best 
summary(model_densTOT)


#Model checking plot 
residuals <- resid(model_densTOT)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_densTOT)
qqnorm(model_densTOT, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#5.4 SDW----
model_sdw <- lme(log(SDW) ~ ELEVATION*NO3.uM.mean + NO2.uM.mean + NH4.uM.mean, data=bio, random=~1|R.replicates/SITE, control=lCtr, correlation= corAR1(),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 


# This last model is the one that work best 
summary(model_sdw)

#Model checking plot 
residuals <- resid(model_sdw)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_sdw)
qqnorm(model_sdw, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#load the data and start with the analyses
library(reshape2)
library(nlme)
punch <- read.csv("data/FAMAR_punching123_05_15_nut.c.csv", header=T, sep=",")
punch <- subset(punch, punch$YEAR>2005)
punch <- subset(punch, punch$YEAR<2015)
table(punch$YEAR)
table(punch$SITE)
#put the correct format
punch$S.date <- as.POSIXct(strptime(punch$S.date, format="%d/%m/%Y")) 

#6.1 LGR----
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

model_lgr <- lme(log(LGR) ~ YEAR + ELEVATION.m * NO2.uM.mean, data=punch, random=~1|R.replicate/SITE, control=lCtr, correlation= corARMA(p=0, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 

punch$fELEVATION <- as.factor(punch$ELEVATION.m)

# This last model is the one that work best 
summary(model_lgr)


ggplot(punch, aes(x = YEAR, y = log(LGR), color = fELEVATION) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15, aes(fill = fELEVATION))

#Model checking plot 
residuals <- resid(model_lgr)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_lgr)
qqnorm(model_lgr, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)







#6.2 LLRW----
#Log of 0 is equal to -INF so we need to put a value over zero for these cases
punch$LLRw[punch$LLRw == 0] <- 0.01
model_llrw <- lme(log(LLRw) ~  YEAR*ELEVATION.m, data=punch, random=~1|R.replicate/SITE, control=lCtr, correlation= corARMA(p=0, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 



summary(model_llrw)


ggplot(punch, aes(x = YEAR, y = log(LLRw), color = fELEVATION) ) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15, aes(fill = fELEVATION))



#Model checking plot 
residuals <- resid(model_llrw)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_llrw)
qqnorm(model_llrw, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)
#This is the model in which normality assumption are less clear. 

#6.3 TurnR----
#For this variable is already normally distributed therefore no need to log-transform it. 
model_turnr <- lme(TurnR ~ ELEVATION.m*NH4.uM.mean + NO2.uM.mean, data=punch, random=~1|R.replicate/SITE, control=lCtr, correlation= corARMA(p=0, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
summary(model_turnr)

#Model checking plot 
residuals <- resid(model_turnr)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_turnr)
qqnorm(model_turnr, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#6.4 SDW----
#For this variable is already normally distributed therefore no need to log-transform it. 
model_sdw <- lme(log(SDW) ~ YEAR + ELEVATION.m*NO2.uM.mean, data=punch, random=~1|R.replicate/SITE, control=lCtr, correlation= corAR1(),
                 weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

summary(model_sdw)

#Model checking plot 
residuals <- resid(model_sdw)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_sdw)
qqnorm(model_sdw, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#6.5 PI----
#For this variable is already normally distributed therefore no need to log-transform it. 
model_pi <- lme(log(PI) ~ ELEVATION.m + NO2.uM.mean + NH4.uM.mean + PO4.uM.mean, data=punch, random=~1|R.replicate/SITE, control=lCtr, correlation= corARMA(p=0, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)
summary(model_pi)

#Model checking plot 
residuals <- resid(model_pi)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_pi)
qqnorm(model_pi, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)
