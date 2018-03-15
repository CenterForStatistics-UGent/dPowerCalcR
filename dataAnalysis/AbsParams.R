# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Calculating parameters affecting the power in empirical data
# of digital PCR absolute quantification experiments
#
# Summary statistics of this data appear in Table 1
# of the main manuscript

# read in the Lievens data and format the data
# data preprocessed with AbsProcessingLievens.R

DataLievens <- read.csv("LievensData.csv")
Unique <- paste(DataLievens$Plate,DataLievens$Target,DataLievens$Conc,DataLievens$Primers,DataLievens$Probe,DataLievens$Anneal,DataLievens$Touchdown,DataLievens$Enhancer,DataLievens$Cycles,DataLievens$Sonication)
DataLievens <- data.frame(Unique,Total=DataLievens$Total,Negative=DataLievens$Negative)

# remove reactions with very few partitions
DataLievens <- DataLievens[DataLievens$Total>1000,]

# check distribution of total number of partitions
hist(DataLievens$Total)

#check unique experiments
Unique <- unique(DataLievens$Unique)

# calculate parameters affecting the power for each
# unique experiment (replicates are considered the same experiment)
est <- c()
est.var <- c()
exp.var <- c()
n.repl <- c()
pA.CV <- c()
for(i in 1:length(Unique)){
	data.subset <- DataLievens[DataLievens$Unique==Unique[i],]
	pA.CV[i] <- sd(data.subset[,3]/data.subset[,2])/mean(data.subset[,3]/data.subset[,2])
	est <- c(est,-log(data.subset[,3]/data.subset[,2]))
	est.var <- c(est.var, rep(var(-log(data.subset[,3]/data.subset[,2])),nrow(data.subset)))
	exp.var <- c(exp.var, (1-exp(log(data.subset[,3]/data.subset[,2])))/(data.subset[,2]*exp(log(data.subset[,3]/data.subset[,2]))))
	n.repl <- c(n.repl, nrow(data.subset))
}

n.replLievens <- n.repl

est <- est[complete.cases(est.var)]
exp.var <- exp.var[complete.cases(est.var)]
est.var <- est.var[complete.cases(est.var)]
interrepLievens <- est.var-exp.var
est.varLievens <- est.var
exp.varLievens <- exp.var



# read in the Jones data
# data preprocessed with AbsProcessingJones.R

DataJones <- read.csv("JonesData.csv")
Unique <- unique(DataJones$Sample)
est <- c()
est.var <- c()
exp.var <- c()
n.repl <- c()
for(i in 1:length(Unique)){
	data.subset <- DataJones[DataJones$Sample==Unique[i],]
	est <- c(est,-log(data.subset[,3]/data.subset[,2]))
	est.var <- c(est.var, rep(var(-log(data.subset[,3]/data.subset[,2])),nrow(data.subset)))
	exp.var <- c(exp.var, (1-exp(log(data.subset[,3]/data.subset[,2])))/(data.subset[,2]*exp(log(data.subset[,3]/data.subset[,2]))))
	n.repl <- c(n.repl, nrow(data.subset))
}
est.varJones <- est.var
exp.varJones <- exp.var


CV <- c()
for(i in 1:length(Unique)){
	data.subset <- DataJones[DataJones$Sample==Unique[i],]
	est <- c(-log(data.subset[,3]/data.subset[,2]))
	est.var <- c(var(-log(data.subset[,3]/data.subset[,2])))
	CV <- c(CV, sqrt(est.var)/mean(est))	
}

#summary statistics
#number of partitions
min(c(DataJones$Total,DataLievens$Total))
median(c(DataJones$Total,DataLievens$Total))
max(c(DataJones$Total,DataLievens$Total))
quantile(c(DataJones$Total,DataLievens$Total),0.05)
#fraction of negatives
frac.neg <-c(DataLievens$Negative/DataLievens$Total,DataJones$Negative/DataJones$Total)
frac.neg <- frac.neg[frac.neg!=1]
min(frac.neg)
median(frac.neg)
max(frac.neg)
quantile(frac.neg, 0.95)
#number of replicates
min(c(n.replLievens,n.repl))
median(c(n.replLievens,n.repl))
max(c(n.replLievens,n.repl))
quantile(c(n.replLievens,n.repl), 0.05)
#interrep
min(c(interrepLievens,est.var-exp.var))
median(c(interrepLievens,est.var-exp.var))
max(c(interrepLievens,est.var-exp.var))
quantile(c(interrepLievens,est.var-exp.var), 0.95)