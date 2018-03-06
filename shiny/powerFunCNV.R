################################################################
####           Functions                                    ####
################################################################


# variance of CNV based on first order taylor expansion. This is the variance that I derived using a first order Taylor expansion for Var(X/Y) with X and Y independently distributed. This is equivalent to the CI-based formulation of the approximation used by Bio-Rad

varCNV2<-function(N,NNegA,delta,r,VarBetween,NB) {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	pA<-NNegA/N
	pB<-pA^(1/delta)
	VarPA<-pA*(1-pA)/N
	VarPB<-pA^(1/delta)*(1-pA^(1/delta))/N
	vPoisson<-(NB*delta/(pA*log(pA)))^2*VarPA+(NB*delta^2/(log(pA)*pB))^2*VarPB
	v<-(vPoisson+VarBetween)/r
	return(v)
}

varBetweenCNV <- function(data,NB=1){
	data.target <- data[data[,3]==TRUE,1:2]
	data.ref <- data[data[,3]==FALSE,1:2]
	vTotal <- var(log(1-data.target[,1]/data.target[,2])/log(1-data.ref[,1]/data.ref[,2]))*NB^2
	vPois <- calcPoisVar(data.target,data.ref)*NB^2
	v <- vTotal - vPois
	return(v)	
}

calcPoisVar <- function(data.subset,ref.subset){
	pA <- mean(1-(data.subset[,1]/(data.subset[,2])), na.rm = TRUE)
	pB <- mean(1-(ref.subset[,1]/(ref.subset[,2])), na.rm = TRUE)
	VarPA<-pA*(1-pA)/mean(data.subset[,2], na.rm = TRUE)
	VarPB<-pB*(1-pB)/mean(ref.subset[,2], na.rm = TRUE)
	delta <- log(pA)/log(pB)
	NB <- 1
	poisVar<-(NB*delta/(pA*log(pA)))^2*VarPA+(NB*delta^2/(log(pA)*pB))^2*VarPB
	return(poisVar)
}


power.ddPCR.fixed<-function(delta=1,r=3,NNeg.fr=0.2,NB=2,alpha=0.05,ylim=c(0.5,1),N=20000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
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
	CN.delta <- NB*delta 
	sd.CNV<-sapply(N,function(n,NNegA.fr,delta,NB) {
		sqrt(varCNV2(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,NB=NB))}, NNegA.fr=NNeg.fr,delta=delta,NB=NB)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(CN.delta-NB)/sd.CNV))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(CN.delta-NB)/sd.CNV))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(CN.delta-NB)/sd.CNV)+pnorm(qnorm(alpha/2)-(CN.delta-NB)/sd.CNV))
	}
	return(pwr)
}

power.ddPCR<-function(delta=1,r=3,NNeg.fr=0.2,NB=2,alpha=0.05,ylim=c(0.5,1),Nmin=5000,Nmax=20000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
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
	CN.delta <- NB*delta 
	N<-seq(Nmin,Nmax,by=(Nmax-Nmin)/1000) # sequence of total number of droplets for which the power will be calculated
	sd.CNV<-sapply(N,function(n,NNegA.fr,delta,NB) {
		sqrt(varCNV2(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,NB=NB))}, NNegA.fr=NNeg.fr,delta=delta,NB=NB)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(CN.delta-NB)/sd.CNV))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(CN.delta-NB)/sd.CNV))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(CN.delta-NB)/sd.CNV)+pnorm(qnorm(alpha/2)-(CN.delta-NB)/sd.CNV))
	}
	if(add) {
		lines(N,pwr,lty=lty)
	}
	else {
		plot(N,pwr,xlab="Total number of partitions",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.N<-function(delta=1,r=3,NNeg.fr=0.2,NB=2,alpha=0.05,ylim=c(0.5,1),N=15000,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
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
	CN.delta <- NB*delta 
	NNeg.fr<-seq(0.01,0.99,by=0.001) # sequence of fractions of negatives
	sd.CNV<-sapply(NNeg.fr,function(NNegA.fr,n,delta,NB) {
		sqrt(varCNV2(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=interrep,NB=NB))}, n=N,delta=delta,NB=NB)
	if(alternative=="Greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(CN.delta-NB)/sd.CNV))
	} else if(alternative=="Less"){
		pwr<-(pnorm(qnorm(alpha)-(CN.delta-NB)/sd.CNV))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(CN.delta-NB)/sd.CNV)+pnorm(qnorm(alpha/2)-(CN.delta-NB)/sd.CNV))
	}
	#print(paste("maximum power = ",max(pwr),", negative fraction = ",which(pwr==max(pwr))/length(pwr)))
	if(add) {
		lines(NNeg.fr,pwr,lty=lty)
	}
	else {
		plot(NNeg.fr,pwr,xlab="Negative fraction of partitions",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.N.max<-function(delta.low=1.01,delta.high=1.50,step=0.01,r=3,N=14000,NB=2,alpha=0.05,ylim=c(0,1),add=F,lty=1,alternative="two.sided") {
	#
	# delta: effect under H1
	# NNegA.fr: fraction of negatives
	# r: number of replicates
	# Vp: partition volume	 
	# alpha: significance level (one-sided)
	# add=F: add=T --> line added to existing graph
	# lty=1: line type to be used for plotting
	#
	CN.delta <- NB*delta 
	NNeg.fr<-seq(0.01,0.99,by=0.001)
	delta<-seq(delta.low,delta.high,step)
	optimal.frac<-array(0,length(delta))
	pwr<-array(0,length(delta))
	for(i in 1:length(delta)){
		sd.CNV<-sapply(NNeg.fr,function(NNegA.fr,n,delta,NB) {
			sqrt(varCNV2(n,NNegA=n*NNegA.fr,delta=delta,r=r,VarBetween=0.0031,NB=NB))}, n=N,delta=delta[i],NB=NB)
		pwr.temp<-(1-pnorm(qnorm(1-alpha)-(delta[i]-1)/sd.CNV))
		pwr[i]<-max(pwr.temp) # 
		optimal.frac[i]<-which(pwr.temp==max(pwr.temp))/length(pwr.temp)
	}
	par(mfrow=c(1,2))
	plot(delta,pwr,xlab="CNV",ylab="Maximum power achievable",type="l",main="",ylim=ylim,lty=lty)
	plot(delta,optimal.frac,xlab="CNV",ylab="Optimal negative fraction",type="l",main="",ylim=c(0.3,0.35),lty=lty)
	par(mfrow=c(1,1))
}

power.ddPCR.r<-function(delta=1,NNeg.fr=0.2,NB=2,alpha=0.05,ylim=c(0.5,1),N=15000,rmin=2,rmax=10,interrep=0.0031,add=F,lty=1,alternative="Two-sided") {
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
	CN.delta <- NB*delta 
	r<-rmin:rmax # sequence of number of replicates for which the power is computed
	sd.CNV<-sapply(r,function(n,N,NNegA.fr,delta,NB) {
		sqrt(varCNV2(N,NNegA=N*NNegA.fr,delta=delta,r=n,VarBetween=interrep,NB=NB))}, N=N,NNegA.fr=NNeg.fr,delta=delta,NB=NB)
	if(alternative=="greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(CN.delta-NB)/sd.CNV))
	} else if(alternative=="less"){
		pwr<-(pnorm(qnorm(alpha)-(CN.delta-NB)/sd.CNV))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(CN.delta-NB)/sd.CNV)+pnorm(qnorm(alpha/2)-(CN.delta-NB)/sd.CNV))
	}
	if(add) {
		lines(r,pwr,lty=lty)
	}
	else {
		plot(r,pwr,xlab="Number of replicates",ylab="Power",type="l",main="",ylim=ylim,lty=lty)
	}
}


power.ddPCR.VB<-function(delta=1,NNeg.fr=0.2,NB=2,alpha=0.05,ylim=c(0.5,1),N=15000,r=3,interrepmin=0,interrepmax=0.01,add=F,lty=1,alternative="Two-sided") {
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
	CN.delta <- NB*delta 
	VB<-seq(interrepmin,interrepmax,(interrepmax-interrepmin)/1000) # sequence of number of replicates for which the power is computed
	sd.CNV<-sapply(r,function(n,N,NNegA.fr,delta,NB) {
		sqrt(varCNV2(N,NNegA=N*NNegA.fr,delta=delta,r=n,VarBetween=VB,NB=NB))}, N=N,NNegA.fr=NNeg.fr,delta=delta,NB=NB)
	if(alternative=="greater"){
		pwr<-(1-pnorm(qnorm(1-alpha)-(CN.delta-NB)/sd.CNV))
	} else if(alternative=="less"){
		pwr<-(pnorm(qnorm(alpha)-(CN.delta-NB)/sd.CNV))
	} else {
		pwr<-(1-pnorm(qnorm(1-alpha/2)-(CN.delta-NB)/sd.CNV)+pnorm(qnorm(alpha/2)-(CN.delta-NB)/sd.CNV))
	}
	if(add) {
		lines(VB,pwr,lty=lty)
	}
	else {
		plot(VB,pwr,xlab="Between-replicate variability",ylab="Power",type="l",main="",ylim=ylim,lty=lty,xlim=c(0,max(VB)))
	}
}


