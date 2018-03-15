# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Thresholding of the absolute quantification data
# Lievens data


n <- c()
n.neg <- c()

# Dataset.csv as obtained from https://github.com/Gromgorgel/ddPCR
Data <- read.csv("Dataset.csv")

for(i in 2:ncol(Data)){
	cat(i,"\n")
	fluo <- as.numeric(as.character(Data[11:nrow(Data), i]))
	res <- cloudy(fluo)
	n <- c(n, sum(res$droplets[1:2]))
	n.neg <- c(n.neg, res$droplets[2])
}

Plate <- as.character(unlist(unname(Data[1,-1])))
Target <- as.character(unlist(unname(Data[2,-1])))
Conc <- as.character(unlist(unname(Data[3,-1])))
Primers <- as.character(unlist(unname(Data[4,-1])))
Probe <- as.character(unlist(unname(Data[5,-1])))
Anneal <- as.character(unlist(unname(Data[6,-1])))
Touchdown <- as.character(unlist(unname(Data[7,-1])))
Enhancer <- as.character(unlist(unname(Data[8,-1])))
Cycles <- as.character(unlist(unname(Data[9,-1])))
Sonication <- as.character(unlist(unname(Data[10,-1])))

LievensData <- data.frame(Plate=Plate,Target=Target,Conc=Conc,Primers=Primers,Probe=Probe,Anneal=Anneal,Touchdown=Touchdown,Enhancer=Enhancer,Cycles=Cycles,Sonication=Sonication,Total=n,Negative=n.neg)
write.csv(LievensData,"LievensData.csv",row.names=F)