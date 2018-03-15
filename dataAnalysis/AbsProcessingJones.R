# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Thresholding of the absolute quantification data
# Jones data

n <- c()
n.neg <- c()

# data as obtained from
# https://github.com/jacobhurst/definetherain

setwd("definetherain/10e0")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="10")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

setwd("definetherain/10e1")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="10")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

setwd("definetherain/10e2")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="09")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

setwd("definetherain/10e3")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="09")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

setwd("definetherain/10e4")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="08")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

setwd("definetherain/10e5")
Data<-readDpcr(dataType=1, testPattern="11",newPattern="08")
conc <- thresholdTrypsteen(getwd(),negativeData=Data[[1]]$testData,newData=Data[[1]]$newData,reps = 100)

n.neg <- c(n.neg, conc$nNegative[-1])
n <- c(n, conc$nValid[-1])

write.csv(data.frame(Sample=rep(c("10e0","10e1","10e2","10e3","10e4","10e5"),each=4),Total=n,Negative=n.neg),"JonesData.csv",row.names=F)