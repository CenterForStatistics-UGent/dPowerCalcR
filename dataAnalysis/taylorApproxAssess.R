# Author: M Vynck
# Last update: March 15, 2018
#
# Contents:
# Assessment of how well the Taylor approximations work
# (by simulation)
#
# Results given as Figures 1 and 2 in the Supplementary
# Material of the manuscript

set.seed(11)
N <- 15000
nsims <- 100
pA <- seq(0.0001,0.9999,0.001)
lambda <- -log(pA)
est <- c()
var.est <- c()
for(i in 1:length(pA)){
	for(j in 1:nsims){
		data <- rpois(N,lambda[i])
		est[j] <- -log(mean(data<1))
	}
	var.est[i] <- var(est)
	
}

var.taylor <- (1-pA)/(pA*N)
plot(pA,log(var.est),type="l",xlab="Negative fraction of partitions",ylab="Log variance of avg. number of copies per partition",lty=1,col="grey",lwd=1)
lines(pA,log(var.taylor),lty=2,lwd=0.4)
legend("bottomleft",c("Empirical","Taylor approximation"),lty=c(1,2),col=c("grey","black"),lwd=c(1,0.4))

plot(1:length(var.est),var.est/var.taylor)


set.seed(11)
nsims <- 500
reps <- 3
mean.est <- c()
est.mean <- c()
est.mean.temp <- c()
mean.est.temp <- c()
for(i in 1:length(pA)){
	for(j in 1:nsims){
		data <- rnorm(reps,pA[i],pA[i]/100)
		est.mean.temp[j] <- mean(1/data)
		mean.est.temp[j] <- 1/mean(data)
	}	
	mean.est[i] <- mean(mean.est.temp)
	est.mean[i] <- mean(est.mean.temp)
}
plot(pA,lowess(pA,100*(est.mean)/(mean.est)-100)$y,type="l",xlab="Negative fraction of partitions",ylab="Difference (%)",ylim=c(-1,6))

set.seed(12)
nsims <- 500
reps <- 3
mean.est1 <- c()
est.mean1 <- c()
est.mean.temp <- c()
mean.est.temp <- c()
for(i in 1:length(pA)){
	for(j in 1:nsims){
		data <- rnorm(reps,pA[i],pA[i]/20)
		est.mean.temp[j] <- mean(1/data)
		mean.est.temp[j] <- 1/mean(data)
	}	
	mean.est1[i] <- mean(mean.est.temp)
	est.mean1[i] <- mean(est.mean.temp)
}
lines(pA,lowess(pA,100*(est.mean1)/(mean.est1)-100)$y,lty=2)

set.seed(13)
nsims <- 500
reps <- 3
mean.est2 <- c()
est.mean2 <- c()
est.mean.temp <- c()
mean.est.temp <- c()
for(i in 1:length(pA)){
	for(j in 1:nsims){
		data <- rnorm(reps,pA[i],pA[i]/10)
		est.mean.temp[j] <- mean(1/data)
		mean.est.temp[j] <- 1/mean(data)
	}	
	mean.est2[i] <- mean(mean.est.temp)
	est.mean2[i] <- mean(est.mean.temp)
}
lines(pA,lowess(pA,100*(est.mean2)/(mean.est2)-100)$y,lty=3)

set.seed(14)
nsims <- 500
reps <- 3
mean.est3 <- c()
est.mean3 <- c()
est.mean.temp <- c()
mean.est.temp <- c()
for(i in 1:length(pA)){
	for(j in 1:nsims){
		data <- rnorm(reps,pA[i],pA[i]/5)
		est.mean.temp[j] <- mean(1/data)
		mean.est.temp[j] <- 1/mean(data)
	}	
	mean.est3[i] <- mean(mean.est.temp) 
	est.mean3[i] <- mean(est.mean.temp)
}
lines(pA,lowess(pA,100*(est.mean3)/(mean.est3)-100)$y,lty=4)

legend("topleft",c("CV 1%","CV 5%","CV 10%","CV 20%"),lty=c(1,2,3,4))