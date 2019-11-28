#load the data and start with the analyses
library(reshape2)
library(nlme)
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(DHARMa)

data <- read.csv("data/meadow.year.summary.csv", header=T, sep=",")
data <- subset(data, data$YEAR>2005)
data <- subset(data, data$YEAR<2015)
table(data$YEAR)
table(data$SITE)

#create a new variable to assign a specific plot
data$PLOT <- paste(data$SITE, data$R.replicates, sep="_") 

data$fELEVATION <- as.factor(data$fELEVATION) 

#We are going to establish a linear relationship between biomass
#production and the predictors by doing a log of biomass
#here are some parameters to help the models converge
lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

#1.Gross.Prod.q2----
model_gp1 <- lme(Gross.Prod.q2 ~ YEAR, data=data, control=lCtr, random = ~ 1|R.replicates, correlation= corARMA(p=0, q=1), 
                   method='ML')

Model does not work because I need replication, does not work with the means. !!