#do some graphic representation of the main models. 
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

lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

model_ldw <- lme(log(L.DW) ~ SITE + NH4.uM.q2 + NO3.uM.q2 + PO4.uM.q2, data=data, random=~1|PLOT, control=lCtr, correlation= corARMA(p=2, q=3),
                  weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 

n_dat <- expand.grid(SITE=unique(data$SITE),
                      NH4=c(min(data$NH4.uM.q2),
                            max(data$NH4.uM.q2)),
                      NO3=c(min(data$NO3.uM.q2),
                            max(data$NO3.uM.q2)),
                     PO4=c(min(data$PO4.uM.q2),
                           max(data$PO4.uM.q2)))

library(ggplot2)
p <- ggplot(data, aes(x=NH4.uM.q2, y=log(L.DW), colour=SITE))

+
  geom_point(size=3) +
  geom_line(aes(y=predict(model_ldw), group=PLOT, size="Plots"))+
  geom_line(data=n_dat, aes(y=predict(model_ldw, level=0, newdata=n_dat), size="SITES")) +
  scale_size_manual(name="Predictions", values=c("plots"=0.5, "SITES"=3)) +
  theme_bw(base_size=22) 
print(p)

library(ggplot2)
ggplot(tempEf,aes(TRTYEAR, r, group=interaction(site, Myc), col=site, shape=Myc )) + 
  facet_grid(~N) +
  geom_line(aes(y=fit, lty=Myc), size=0.8) +
  geom_point(alpha = 0.3) + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()
