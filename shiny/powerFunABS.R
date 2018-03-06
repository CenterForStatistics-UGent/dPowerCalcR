################################################################
####           Functions                                    ####
################################################################


varABS<-function(N,NNegA,delta,r,VarBetween,Vp) {
	# This function calculates the variance of the CNV estimate based on r replicates
	#
	# based on estimate of between replicate variance
	# based on first order Taylor expansion
	# delta: effect under H1
	# N: total number of droplets
	# NNegA: number of negatives for the target
	# r: number of replicates
	# VarBetween: variance between replicates
	# Vp: partition volume
	#
	pA<-NNegA/N
	VarPA<-pA*(1-pA)/N
	vPoisson<-VarPA/(pA*Vp)^2
	v<-(vPoisson+VarBetween)/r
	return(v)
}

varBetween <- function(N,Npos,Vp){
	VarTotal<-var(-log(1-Npos/N)/Vp)
	NposA <- mean(Npos)
	N <- mean(N)
	pA<-1-NposA/N
	VarPA<-pA*(1-pA)/N
	vPoisson<-VarPA/(pA*Vp)^2
	v<-VarTotal-vPoisson
	if(is.na(v)){ v <- -1}
	return(v)	
}


power.ddPCR.abs.fixed<-function(delta=1,r=3,NNeg.fr=0.2,Vp=0.85,alpha=0.05,ylim=c(0.5,1),N=20000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {

	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	delta=delta+1
	c = -log(NNeg.fr)/Vp
	c.0 = c/delta
	sd.ABS<-sapply(N,function(n,NNegA.fr,delta,Vp) {
		sqrt(varABS(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,Vp=Vp))}, NNegA.fr=NNeg.fr,delta=delta,Vp=Vp)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(c-c.0)/sd.ABS))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(c-c.0)/sd.ABS))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(c-c.0)/sd.ABS)+pnorm(qnorm(alpha/2)-(c-c.0)/sd.ABS))
	}
	return(pwr)
}

power.ddPCR.abs<-function(delta=1,r=3,NNeg.fr=0.2,Vp=0.85,alpha=0.05,ylim=c(0.5,1),Nmin=5000,Nmax=20000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	delta=delta+1
	c = -log(NNeg.fr)/Vp
	c.0 = c/delta
	N<-seq(Nmin,Nmax,by=(Nmax-Nmin)/1000) # sequence of total number of droplets for which the power will be calculated
	sd.ABS<-sapply(N,function(n,NNegA.fr,delta,Vp) {
		sqrt(varABS(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,Vp=Vp))}, NNegA.fr=NNeg.fr,delta=delta,Vp=Vp)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(c-c.0)/sd.ABS))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(c-c.0)/sd.ABS))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(c-c.0)/sd.ABS)+pnorm(qnorm(alpha/2)-(c-c.0)/sd.ABS))
	}
	if(add) {
		lines(N,pwr,lty=lty)
	}
	else {
		plot(N,pwr,xlab="Total number of partitions",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.abs.N<-function(delta=1,r=3,NNeg.fr=0.2,Vp=0.85,alpha=0.05,ylim=c(0.5,1),N=15000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	delta=delta+1
	NNeg.fr<-seq(0.001,0.999,by=0.001) # sequence of fractions of negatives
	c = -log(NNeg.fr)/Vp
	c.0 = c/delta
	sd.ABS<-sapply(NNeg.fr,function(NNegA.fr,n,delta,Vp) {
		sqrt(varABS(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,Vp=Vp))}, n=N,delta=delta,Vp=Vp)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(c-c.0)/sd.ABS))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(c-c.0)/sd.ABS))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(c-c.0)/sd.ABS)+pnorm(qnorm(alpha/2)-(c-c.0)/sd.ABS))
	}
	#print(paste("maximum power = ",max(pwr),", negative fraction = ",which(pwr==max(pwr))/length(pwr)))
	if(add) {
		lines(NNeg.fr,pwr,lty=lty)
	}
	else {
		plot(NNeg.fr,pwr,xlab="Negative fraction of partitions",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.abs.N.max<-function(delta.low=1.01,delta.high=1.50,step=0.01,r=3,N=14000,Vp=0.85,alpha=0.05,ylim=c(0,1),add=F,lty=1,alternative="two.sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	NNeg.fr<-seq(0.01,0.99,by=0.001)
	delta<-seq(delta.low,delta.high,step)
	optimal.frac<-array(0,length(delta))
	pwr<-array(0,length(delta))
	for(i in 1:length(delta)){
		sd.ABS<-sapply(NNeg.fr,function(NNegA.fr,n,delta,Vp) {
			sqrt(varABS(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=0.0031,Vp=Vp))}, n=N,delta=delta[i],Vp=Vp)
		pwr.temp<-(1-pnorm(qnorm(1-alpha)-(delta[i]-1)/sd.ABS))
		pwr[i]<-max(pwr.temp) # 
		optimal.frac[i]<-which(pwr.temp==max(pwr.temp))/length(pwr.temp)
	}
	par(mfrow=c(1,2))
	plot(delta,pwr,xlab="CNV",ylab="Maximum power achievable",type="l",main="",ylim=ylim,lty=lty)
	plot(delta,optimal.frac,xlab="CNV",ylab="Optimal negative fraction",type="l",main="",ylim=c(0.3,0.35),lty=lty)
	par(mfrow=c(1,1))
}

power.ddPCR.abs.r<-function(delta=1,NNeg.fr=0.2,Vp=0.85,alpha=0.05,ylim=c(0.5,1),N=15000,rmin=2,rmax=10,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	delta <- delta+1
	c = -log(NNeg.fr)/Vp
	c.0 = c/delta
	r<-rmin:rmax # sequence of number of replicates for which the power is computed
	sd.ABS<-sapply(r,function(n,N,NNegA.fr,delta,Vp) {
		sqrt(varABS(N,NNegA=N*NNegA.fr,delta=delta,r=n,VarBetween=interrep,Vp=Vp))}, N=N,NNegA.fr=NNeg.fr,delta=delta,Vp=Vp)
	if(alternative=="greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(c-c.0)/sd.ABS))
	} else if(alternative=="less"){
		pwr<-(pnorm(qnorm(alpha)-(c-c.0)/sd.ABS))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(c-c.0)/sd.ABS)+pnorm(qnorm(alpha/2)-(c-c.0)/sd.ABS))
	}
	if(add) {
		lines(r,pwr,lty=lty)
	}
	else {
		plot(r,pwr,xlab="Number of replicates",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.abs.VB<-function(delta=1,NNeg.fr=0.2,Vp=0.85,alpha=0.05,ylim=c(0.5,1),N=15000,r=3,interrepmin=0,interrepmax=0.01,add=F,lty=1,alternative="Two-sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	delta=delta+1
	c = -log(NNeg.fr)/Vp
	c.0 = c/delta
	VB<-seq(interrepmin,interrepmax,(interrepmax-interrepmin)/1000) # sequence of number of replicates for which the power is computed
	sd.ABS<-sapply(r,function(n,N,NNegA.fr,delta,Vp) {
		sqrt(varABS(N,NNegA=N*NNegA.fr,delta=delta,r=n,VarBetween=VB,Vp=Vp))}, N=N,NNegA.fr=NNeg.fr,delta=delta,Vp=Vp)
	if(alternative=="greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(c-c.0)/sd.ABS))
	} else if(alternative=="less"){
		pwr<-(pnorm(qnorm(alpha)-(c-c.0)/sd.ABS))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(c-c.0)/sd.ABS)+pnorm(qnorm(alpha/2)-(c-c.0)/sd.ABS))
	}
	if(add) {
		lines(VB,pwr,lty=lty)
	}
	else {
		plot(VB,pwr,xlab="Between-replicate variability",ylab="Power",type="l",main="",ylim=ylim,lty=lty,xlim=c(0,max(VB)))
	}
}


