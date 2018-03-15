# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Calculating parameters affecting the power in empirical data
# of digital PCR relative quantification experiments
#
# Summary statistics of this data appear in Table 2
# of the main manuscript

# read in the Vynck data and format the data
# data from supplementary of Vynck et al. (2016). Flexible
# analysis of digital PCR experiments. Biomol Detect Quantif,
# 9, 1-13.

data <- read.csv("NIPDmerged.csv",sep=";")
head(data)

# function to calculate the stochastic variance
calcPoisVar <- function(data.subset,ref.subset){
	pA <- mean(unlist(unname(c(data.subset[2]/rowSums(data.subset)))), na.rm = TRUE)
	pB <- mean(unlist(unname(c(ref.subset[2]/rowSums(ref.subset)))), na.rm = TRUE)
	VarPA<-pA*(1-pA)/mean(unname(rowSums(data.subset)), na.rm = TRUE)
	VarPB<-pB*(1-pB)/mean(unname(rowSums(ref.subset)), na.rm = TRUE)
	delta <- log(pA)/log(pB)
	NB <- 2
	poisVar<-(NB*delta/(pA*log(pA)))^2*VarPA+(NB*delta^2/(log(pA)*pB))^2*VarPB
	return(poisVar)
}

# calculate parameters affecting the power for each
# of the replicated experiments
n <- c()
frac.neg <- c()
RQ <- c()
CN <- c()
CN.mean <- c()
CN.var <- c()
pois.var <- c()
n.repl <- c()

for(i in 1:14){
	for(j in 1:14){
		cols <- c(i*2-1, i*2) + 1
		rows <- c((j*3-2):(j*3))
		data.subset <- data[rows, cols]
		ref.subset <- data[rows, 28:29]
		n <- c(n, unname(rowSums(data.subset)))
		frac.neg <- c(frac.neg, unlist(unname(c(data.subset[2]/rowSums(data.subset)))))
		CN.subset <- log(unlist(unname(c(data.subset[2]/rowSums(data.subset)))))/log(unlist(unname(c(ref.subset[2]/rowSums(ref.subset)))))*2
		CN <- c(CN, CN.subset)
		CN.mean <- c(CN.mean, mean(CN.subset, na.rm = TRUE))
		CN.var <- c(CN.var, var(CN.subset, na.rm = TRUE))
		pois.var <- c(pois.var, calcPoisVar(data.subset, ref.subset))
		n.repl <- c(n.repl, sum(complete.cases(data.subset[,1])))
	}
}

# combine data in data frame
n <- n[complete.cases(n)]
frac.neg <- frac.neg[complete.cases(frac.neg)]
frac.neg <- frac.neg[frac.neg != 1]
CN.df <- data.frame(mean=CN.mean, var=CN.var)
CN.df <- CN.df[complete.cases(CN.var), ]
pois.var <- pois.var[complete.cases(CN.var)]

CN.var <- CN.var[complete.cases(CN.var)]
pois.var <- pois.var[CN.var != 0]

CN.df <- CN.df[CN.var != 0, ]
CN.var <- CN.var[CN.var != 0]


# summary stats
min(n)
median(n)
max(n)
quantile(n, 0.05)
min(frac.neg)
median(frac.neg)
max(frac.neg)
quantile(frac.neg,0.95)
min(n.repl)
median(n.repl)
max(n.repl)
quantile(n.repl,0.05)
betweenvar <- CN.var - pois.var
min(betweenvar)
median(betweenvar)
max(betweenvar)
quantile(betweenvar, 0.95)