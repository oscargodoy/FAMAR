#load the data and start with the analyses
library(reshape2)
library(nlme)

data <- read.csv("data/famar_clean.csv", header=T, sep=",")
#We are going to establish a linear relationship between biomass
#production and the predictors by doing a log of biomass
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

#1.leaf dry weight (L.DW)----
model_ldw1 <- lme(log(L.DW) ~ SITE, data=data, random=~1|REPLICATE, control=lCtr, correlation= corAR1(),
                  method='REML',na.action=na.omit)

model_ldw2 <- lme(log(L.DW) ~ SITE + NH4.uM.mean, data=data, random=~1|REPLICATE, control=lCtr, correlation= corAR1(),
                  method='REML',na.action=na.omit)

model_ldw3 <- lme(log(L.DW) ~ SITE + NH4.uM.mean + NO3.uM.mean, data=data, random=~1|REPLICATE, control=lCtr, correlation= corAR1(),
                  method='REML',na.action=na.omit)

AIC(model_ldw1, model_ldw2, model_ldw3)
#Model checking plot 
residuals <- resid(model_ldw3)
plot(fitted(model_ldw3), residuals)
abline(0,0)
plot(fitted(model_ldw3), log(data$L.DW)) # this lines does not work yet
qqnorm(residuals)
qqline(residuals)

#2.below ground dry weight (BG.DW)----

#model1 <- lme(log(L.DW) ~ SITE, data=data, random=~1|REPLICATE, control=lCtr, correlation=corAR1(form=~S.date), method='REML', na.action=na.omit)

