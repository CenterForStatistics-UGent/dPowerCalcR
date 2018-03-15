# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Contribution (relative) of sources of variation to
# the total variance
#
# Figures are given as Figures 6 and 7 in the 
# Supplementary Material of the manuscript


#plot for relative quantification

source('CNVparams.R')
between.var <- c()
for(i in 1:length(pois.var)){
	between.var[i] <- max(0,CN.var[i]-pois.var[i])
}
df.var <- data.frame(Contribution=c(between.var/CN.var,1-between.var/CN.var),Source=c(rep("Between-replicate",length(between.var)),rep("Stochastic",length(between.var))))
df.var <- df.var[c(order(between.var/CN.var),order(between.var/CN.var)+length(between.var)),]
df.var$ID <- c(1:length(between.var),1:length(between.var))
df.var$Contribution <- round(df.var$Contribution*100)

library(ggplot2)
ggplot(data=df.var, aes(x=ID,y=Contribution,fill=Source))+geom_bar(stat="identity",width=1)+xlab("Sample")+ylab("Contribution (%)")



#plot for absolute quantification
source('ABSparams.R')

between.var <- c()
for(i in 1:length(est.var)){
	between.var[i] <- max(0,est.varJones[i]-exp.varJones[i])
}

between.var2 <- c()
for(i in 1:length(est.varLievens)){
	between.var2[i] <- max(0,est.varLievens[i]-exp.varLievens[i])
}

between.var <- c(between.var,between.var2)
total.var <- c(est.var,est.varLievens)
between.var <- (between.var/total.var)[order(between.var/total.var)]
between.var <- between.var[-391] #remove 0 val
df.var.abs <- data.frame(Contribution=c(between.var,1-between.var),Source=c(rep("Between-replicate",length(between.var)),rep("Stochastic",length(between.var))))
df.var.abs$ID <- c(1:length(between.var),1:length(between.var))
df.var.abs$Contribution <- round(df.var.abs$Contribution*100)
ggplot(data=df.var.abs, aes(x=ID,y=Contribution,fill=Source))+geom_bar(stat="identity",width=1)+xlab("Sample")+ylab("Contribution (%)")
