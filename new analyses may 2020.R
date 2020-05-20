#load the data and start with the analyses

library(reshape2)
library(nlme)
library(ggplot2)
library(MuMIn)
library(multcomp)
library(sjPlot)
library(sjmisc)


# 1. Demographic variables ####
data <- read.csv("data/FAMAR_biom_FINAL.csv", header=T, sep=";")
data <- subset(data, data$YEAR>=2005)
data <- subset(data, data$YEAR<=2015)
table(data$YEAR)
table(data$SITE)

#create a new variable to assign a specific plot
data$PLOT <- paste(data$SITE, data$R.replicates, sep="_") 
data2 <- na.omit(data) #remove NA to later perform model selection.

#We are going to establish a linear relationship between biomass
#production and the predictors by doing a log of biomass
#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

# L.DW: Aboveground (hojas)----
#We first model the random part 
model_ldw1 <- lme(log(L.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

model_ldw2 <- lme(log(L.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

model_ldw3 <- lme(log(L.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

model_ldw4 <- lme(log(L.DW) ~ 1, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit)

AIC(model_ldw1, model_ldw2, model_ldw3, model_ldw4)
#model 3 perform best 

options(na.action = "na.fail")

model_ldw <- lme(log(L.DW) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_ldw, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_ldw_best <- lme(log(L.DW) ~ SITE + SST.C + 1, data = data2, 
                      random = ~1 | PLOT, correlation = corARMA(p = 2, q = 2), 
                      weights = varIdent(form = ~1 | SEASON), method = "REML", 
                      control = lCtr)
summary(model_ldw_best)

#Model checking plot 
residuals <- resid(model_ldw_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_ldw_best)
plot(model_ldw_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_ldw_best, Subject ~ resid(.))
qqnorm(model_ldw_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_ldw_best, linfct = mcp(SITE = "Tukey")))
# all the SITES CN1, CN2, and CN3 are differet.

#plot the main results.
plot_model(model_ldw_best, type = "pred", terms = c("SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

# BG.DW: Belowground (raices y rizomas)----
#We first model the random part 
model_bgdw1 <- lme(log(BG.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_bgdw2 <- lme(log(BG.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_bgdw3 <- lme(log(BG.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_bgdw4 <- lme(log(BG.DW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_bgdw1, model_bgdw2, model_bgdw3, model_bgdw4)
#model 1 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_bgdw <- lme(log(BG.DW) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_bgdw, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_bgdw_best <- lme(log(BG.DW) ~ SITE + SST.C + YEAR + SITE*SST.C + SITE*YEAR + SST.C*YEAR + 1, data = data2, 
                      random = ~1 | PLOT, correlation = corARMA(p = 1, q = 1), 
                      weights = varIdent(form = ~1 | SEASON), method = "REML", 
                      control = lCtr)
summary(model_bgdw_best)

#Model checking plot 
residuals <- resid(model_bgdw_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_bgdw_best)
plot(model_bgdw_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_bgdw_best, YEAR ~ resid(.))
qqnorm(model_bgdw_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_bgdw_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_bgdw_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction
lsmeans(model_bgdw_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction
# significant differences between CN3 with CN2 and CN1, there is a difference happening in 2010, 

#plot the main results.
plot_model(model_bgdw_best, type = "pred", terms = c("YEAR", "SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

#dens.TOT: Shoot density ----
#We first model the random part 
model_denstot1 <- lme(log(dens.TOT) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_denstot2 <- lme(log(dens.TOT) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_denstot3 <- lme(log(dens.TOT) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_denstot4 <- lme(log(dens.TOT) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_denstot1, model_denstot2, model_denstot3, model_denstot4)
#model 3 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_denstot <- lme(log(dens.TOT) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_denstot, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_denstot_best <- lme(log(dens.TOT) ~ SITE + SST.C + YEAR + SITE*SST.C + SITE*YEAR + SST.C*YEAR + 1, data = data2, 
                       random = ~1 | PLOT, correlation = corARMA(p = 2, q = 2), 
                       weights = varIdent(form = ~1 | SEASON), method = "REML", 
                       control = lCtr)
summary(model_denstot_best)

#Model checking plot 
residuals <- resid(model_denstot_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_denstot_best)
plot(model_denstot_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_denstot_best, YEAR ~ resid(.))
qqnorm(model_denstot_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_denstot_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_denstot_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction
lsmeans(model_denstot_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction
# marginally differences between CN3 with CN2 and CN1, there is a difference happening in 2010, 

#plot the main results.
plot_model(model_denstot_best, type = "pred", terms = c("YEAR", "SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

#SDW: Shoot size ----
#Esta variable se ha calculado como L.DW/dens.TOT

#We first model the random part 
model_sdw1 <- lme(log(SDW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                      weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_sdw2 <- lme(log(SDW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                      weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_sdw3 <- lme(log(SDW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                      weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_sdw4 <- lme(log(SDW) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                      weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_sdw1, model_sdw2, model_sdw3, model_sdw4)
#model 3 perform best 

#We then model the fixed part
options(na.action = "na.fail")

model_sdw <- lme(log(SDW) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                     weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_sdw, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_sdw_best <- lme(log(SDW) ~ SITE + SST.C + SITE*SST.C + 1, data = data2, 
                          random = ~1 | PLOT, correlation = corARMA(p = 2, q = 2), 
                          weights = varIdent(form = ~1 | SEASON), method = "REML", 
                          control = lCtr)
summary(model_sdw_best)

#Model checking plot 
residuals <- resid(model_sdw_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_sdw_best)
plot(model_sdw_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_sdw_best, YEAR ~ resid(.))
qqnorm(model_sdw_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_sdw_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_sdw_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction
# differences across site at SST.C equal to 57.78 

#plot the main results.
plot_model(model_sdw_best, type = "pred", terms = c("SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

# rm(list = ls())
# 2. Dynamic variables ####
data <- read.csv("data/FAMAR_punching_FINAL.csv", header=T, sep=";")
data <- subset(data, data$YEAR>=2005)
data <- subset(data, data$YEAR<=2015)
table(data$YEAR)
table(data$SITE)

#create a new variable to assign a specific plot
data$PLOT <- paste(data$SITE, data$R.replicates, sep="_") 
data2 <- na.omit(data) #remove NA to later perform model selection.

#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

# LGR: Leaf growth rate (mgDW/shoot/d) ----
#We first model the random part 
model_lgr1 <- lme(log(LGR) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_lgr2 <- lme(log(LGR) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_lgr3 <- lme(log(LGR) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_lgr4 <- lme(log(LGR) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_lgr1, model_lgr2, model_lgr3, model_lgr4)
#model 1 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_lgr <- lme(log(LGR) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_lgr, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_lgr_best <- lme(log(LGR) ~ SITE + SST.C + YEAR + SITE*SST.C + SITE*YEAR + SST.C*YEAR + SITE*SST.C*YEAR + 1, 
                      data = data2, random = ~1 | PLOT, correlation = corARMA(p = 1, q = 1), 
                      weights = varIdent(form = ~1 | SEASON), method = "REML", 
                      control = lCtr)
summary(model_lgr_best)

#Model checking plot 
residuals <- resid(model_lgr_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_lgr_best)
plot(model_lgr_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_lgr_best, YEAR ~ resid(.))
qqnorm(model_lgr_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_lgr_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_lgr_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction
lsmeans(model_lgr_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction
# differences between Cn1 and CN2

#plot the main results.
plot_model(model_lgr_best, type = "pred", terms = c("YEAR", "SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

# LLRw: Leaf loss rate (mgDW/shoot/d). ----
#We first model the random part 
model_llrw1 <- lme(log(LLRw +0.1) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)
## add the 0.1 to avoid problems with the log itself. 

model_llrw2 <- lme(log(LLRw +0.1) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_llrw3 <- lme(log(LLRw +0.1) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_llrw4 <- lme(log(LLRw +0.1) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_llrw1, model_llrw2, model_llrw3, model_llrw4)
#model 1 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_llrw <- lme(log(LLRw +0.1) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_llrw, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_llrw_best <- lme(log(LLRw +0.1) ~ SITE + SST.C + YEAR + SITE*SST.C + SITE*YEAR + SST.C*YEAR + 1, 
                      data = data2, random = ~1 | PLOT, correlation = corARMA(p = 1, q = 1), 
                      weights = varIdent(form = ~1 | SEASON), method = "REML", 
                      control = lCtr)
summary(model_llrw_best)

#Model checking plot 
residuals <- resid(model_llrw_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_llrw_best)
plot(model_llrw_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_llrw_best, YEAR ~ resid(.))
qqnorm(model_llrw_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_llrw_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_llrw_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction
lsmeans(model_llrw_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction


#plot the main results.
plot_model(model_llrw_best, type = "pred", terms = c("YEAR", "SST.C[50,57,69]", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 

# PI: plastochrone index (d) ----
#We first model the random part 
model_pi1 <- lme(log(PI) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=1),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_pi2 <- lme(log(PI) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=1, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_pi3 <- lme(log(PI) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_pi4 <- lme(log(PI) ~ 1, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_pi1, model_pi2, model_pi3, model_pi4)
#model 4 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_pi <- lme(log(PI) ~ SITE*YEAR*SST.C, data=data2, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_pi, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_pi_best <- lme(log(PI) ~ SITE + YEAR + SITE*YEAR + 1, 
                       data = data2, random = ~1 | PLOT, correlation = corARMA(p = 2, q = 3), 
                       weights = varIdent(form = ~1 | SEASON), method = "REML", 
                       control = lCtr)
summary(model_pi_best)

#Model checking plot 
residuals <- resid(model_pi_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_pi_best)
plot(model_pi_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_pi_best, YEAR ~ resid(.))
qqnorm(model_pi_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_pi_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_pi_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction

#plot the main results.
plot_model(model_pi_best, type = "pred", terms = c("YEAR", "SITE"))
#Values of SST.C, corresponds to 1quartile, median, and third quartile. 
#PI fit does not work

# 3. Dynamic variables at the prairie level ####

# rm(list = ls())

data <- read.csv("data/FAMAR_meadow_FINAL.csv", header=T, sep=";")
data <- subset(data, data$YEAR>=2005)
data <- subset(data, data$YEAR<=2015)
table(data$YEAR)
table(data$SITE)

data2 <- na.omit(data) #remove NA to later perform model selection.

#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

# PROD.S (gDW.m2.season). Producción 'bruta' ----
#We first model the random part 
model_prods1 <- lme(log(Prod.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=1),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_prods2 <- lme(log(Prod.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_prods3 <- lme(log(Prod.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_prods4 <- lme(log(Prod.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=1),
                 weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_prods1, model_prods2, model_prods3, model_prods4)
#model 3 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_prods <- lme(log(Prod.S) ~ SITE*YEAR*SST.C, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_prods, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_prods_best <- lme(log(Prod.S) ~ SITE + SST.C + YEAR + SITE*SST.C + SITE*YEAR + SST.C*YEAR  + 1, 
                     data = data2, random = ~1 | SITE, correlation = corARMA(p = 2, q = 2), 
                     weights = varIdent(form = ~1 | SEASON), method = "REML", 
                     control = lCtr)
summary(model_prods_best)

#Model checking plot 
residuals <- resid(model_prods_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_prods_best)
plot(model_prods_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_prods_best, YEAR ~ resid(.))
qqnorm(model_prods_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
summary(glht(model_prods_best, linfct = mcp(SITE = "Tukey"))) #tukey for a single factor
lsmeans(model_prods_best, list(pairwise ~ SITE*YEAR), adjust = "tukey") #tukey for an interaction
lsmeans(model_prods_best, list(pairwise ~ SITE*SST.C), adjust = "tukey") #tukey for an interaction

#plot the main results.
plot_model(model_prods_best, type = "pred", terms = c("YEAR", "SITE"))
plot_model(model_prods_best, type = "pred", terms = c("SST.C[50,57,69]", "SITE"))


# LOSS.S (gDW.m2.season). Perdidas ----
#We first model the random part 
model_losss1 <- lme(log(Loss.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=1),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_losss2 <- lme(log(Loss.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_losss3 <- lme(log(Loss.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_losss4 <- lme(log(Loss.S) ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=3),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_losss1, model_losss2, model_losss3)
#model 3 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_losss <- lme(log(Loss.S) ~ SITE*YEAR*SST.C, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_losss, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_losss_best <- lme(log(Loss.S) ~  SST.C + YEAR + SST.C*YEAR  + 1, 
                        data = data2, random = ~1 | SITE, correlation = corARMA(p = 2, q = 2), 
                        weights = varIdent(form = ~1 | SEASON), method = "REML", 
                        control = lCtr)
summary(model_losss_best)

#Model checking plot 
residuals <- resid(model_losss_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_losss_best)
plot(model_losss_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_losss_best, YEAR ~ resid(.))
qqnorm(model_losss_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
lsmeans(model_losss_best, list(pairwise ~ SST.C*YEAR), adjust = "tukey") #tukey for an interaction

#plot the main results.
plot_model(model_losss_best, type = "pred", terms = c("YEAR", "SITE"))
plot_model(model_losss_best, type = "pred", terms = c("YEAR","SST.C[50,57,69]"))

# NET.S (gDW.m2.season). Producción 'neta' ----

#We first model the random part 
model_nets1 <- lme(Net.S ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=1),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_nets2 <- lme(Net.S ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=1, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_nets3 <- lme(Net.S ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

model_nets4 <- lme(Net.S ~ 1, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=3),
                    weights = varIdent(form= ~ 1 | SEASON), method='REML',na.action=na.omit)

AIC(model_nets1, model_nets2, model_nets3)
#model 3 perform best 

#We then model the fixed part 
options(na.action = "na.fail")

model_nets <- lme(Net.S ~ SITE*YEAR*SST.C, data=data2, random=~1|SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                   weights = varIdent(form= ~ 1 | SEASON), method='REML')
ms2 <- dredge(model_nets, trace = TRUE, rank = "AICc", REML = FALSE)
(attr(ms2, "rank.call"))
fmList <- get.models(ms2, 1:4)
summary(model.avg(fmList))

#Best model 
model_nets_best <- lme(Net.S ~  SST.C + 1, 
                        data = data2, random = ~1 | SITE, correlation = corARMA(p = 2, q = 2), 
                        weights = varIdent(form = ~1 | SEASON), method = "REML", 
                        control = lCtr)
summary(model_nets_best)

#Model checking plot 
residuals <- resid(model_nets_best)
hist(residuals, prob=TRUE, col="darkgray")
Range = seq(min(residuals), max(residuals), length = length(residuals))
Norm = dnorm(Range, mean = mean(residuals), sd = sd(residuals))
lines(Range, Norm, col = "blue", lwd = 2)

plot(model_nets_best)
plot(model_nets_best, resid(., type = "p") ~ fitted(.) | SITE, abline = 0)
plot(model_nets_best, YEAR ~ resid(.))
qqnorm(model_nets_best, ~ranef(., level=1))
qqnorm(residuals)
qqline(residuals)

#Perform tukey test 
#no fixed factor selected 
#plot the main results.
plot_model(model_nets_best, type = "pred", terms = "SST.C")




