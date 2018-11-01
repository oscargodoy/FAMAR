##Final analyses ----

#load the data and start with the analyses
library(reshape2)
library(nlme)
library(visreg)

bio <- read.csv("data/famar-Biom_123_05_15_minVAR.csv", header=T, sep=",")
bio <- subset(bio, bio$YEAR>2005)
bio <- subset(bio, bio$YEAR<2015)
table(bio$YEAR)
table(bio$SITE)
#put the correct format
bio$S.date <- as.POSIXct(strptime(bio$S.date, format="%d/%m/%Y")) 

lCtr <- lmeControl(maxIter = 5000, msMaxIter = 5000, tolerance = 1e-8, niterEM = 2500, msMaxEval = 2000)

#Example with LDW
model_ldw <- lme(log(L.DW) ~ SITE*NH4.uM.mean + NO3.uM.mean, data=bio, random=~1|R.replicates/SITE, control=lCtr, correlation= corARMA(p=2, q=2),
                 weights = varIdent(form= ~ 1 | SEASON), method='ML',na.action=na.omit) 
summary(model_ldw)
visreg(model_ldw, "NH4.uM.mean", by="SITE", gg=TRUE, xlab="NH4 concentration (uM)", 
       ylab="Leaf dry weight (g) log transformed", layout=c(3,1), fill.par=list(col="#008DFF33"))


#####ESTO ES SOLO PARA RECUPERAR SI ACASO. 

#Obtained fitted data 
fit_model_ldw <- effect(term="SITE*NH4.uM.q2", 
                        xlevels= list(NH4.uM.q2=c(as.numeric(summary(data$NH4.uM.q2)[2]), 
                                                  as.numeric(summary(data$NH4.uM.q2)[5]))), 
                        mod=model_ldw)
fit_model_ldw <- as.data.frame(fit_model_ldw)
fit_model_ldw
fit_model_ldw$NH4.uM.q2 <- as.factor(fit_model_ldw$NH4.uM.q2)
ggplot(fit_model_ldw, aes(x=SITE, y=fit, colour=NH4.uM.q2, group=NH4.uM.q2)) +
        geom_point() + 
         geom_line(size=1.2) +
         geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=NH4.uM.q2),alpha=0.3) + 
         labs(title = "Leaf dry weight", x= "SITES", y="Leaf Dry Biomass (g)", color="Median NH4", fill="Median NH4") + theme_classic() + theme(text=element_text(size=20))

n_dat <- expand.grid(SITE=unique(data$SITE),
                      NH4=c(min(data$NH4.uM.q2),
                            max(data$NH4.uM.q2)),
                      NO3=c(min(data$NO3.uM.q2),
                            max(data$NO3.uM.q2)),
                     PO4=c(min(data$PO4.uM.q2),
                           max(data$PO4.uM.q2)))




+
  geom_point(size=3) +
  geom_line(aes(y=predict(model_ldw), group=PLOT, size="Plots"))+
  geom_line(data=n_dat, aes(y=predict(model_ldw, level=0, newdata=n_dat), size="SITES")) +
  scale_size_manual(name="Predictions", values=c("plots"=0.5, "SITES"=3)) +
  theme_bw(base_size=22) 
print(p)

##Eje y representar ldw, eje x nutrientes, y cada sitio con sus puntos, linea media y barras de error.

