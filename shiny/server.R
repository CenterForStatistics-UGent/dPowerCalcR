library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)
source("powerFunCNV.R")
source("powerFunABS.R")
numSamples<-5

server<-function(input, output)({

paramsABS <- eventReactive(input$calcResultsABS,{
		list(replicates=as.numeric(input$replicates),
			partitions=as.numeric(input$partitions),
			interrep=as.numeric(input$interrep),
			effect=as.numeric(input$effect),
			fraction=as.numeric(input$fraction),
			significance=as.numeric(input$significance),
			Vp=as.numeric(input$Vp),
			alternative=input$Alt,
			minpart=as.numeric(input$minpart),
			maxpart=as.numeric(input$maxpart),
			minrep=as.numeric(input$minrep),
			maxrep=as.numeric(input$maxrep),
			mineffect=as.numeric(input$mineffect),
			maxeffect=as.numeric(input$maxeffect),
			minrepl=as.numeric(input$minrepl),
			maxrepl=as.numeric(input$maxrepl))
	}
)

paramsREL <- eventReactive(input$calcResultsREL,{
		list(replicates=as.numeric(input$replicatesREL),
			partitions=as.numeric(input$partitionsREL),
			interrep=as.numeric(input$interrepREL),
			effect=as.numeric(input$effectREL),
			fraction=as.numeric(input$fractionREL),
			significance=as.numeric(input$significanceREL),
			NB=as.numeric(input$NBREL),
			alternative=input$AltREL,
			minpart=as.numeric(input$minpartREL),
			maxpart=as.numeric(input$maxpartREL),
			minrep=as.numeric(input$minrepREL),
			maxrep=as.numeric(input$maxrepREL),
			mineffect=as.numeric(input$mineffectREL),
			maxeffect=as.numeric(input$maxeffectREL),
			minrepl=as.numeric(input$minreplREL),
			maxrepl=as.numeric(input$maxreplREL))
	}
)

output$textABS <- renderText({
	pwr <- power.ddPCR.abs.fixed(delta=paramsABS()[["effect"]],
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1)
	)
	paste0("The power for the specified parameters is: ", round(pwr*100,2), "%.")
})

output$textREL <- renderText({
	pwr <- power.ddPCR.fixed(delta=paramsREL()[["effect"]],
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1)
	)
	paste0("The power for the specified parameters is: ", round(pwr*100,2), "%.")
})

output$plotABS <- renderPlot({
	par(mfrow=c(2,2))
	#effect size / number of partitions
	power.ddPCR.abs(delta=paramsABS()[["mineffect"]],
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			Nmin=paramsABS()[["minpart"]],
			Nmax=paramsABS()[["maxpart"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.abs(delta=1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]]),
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			Nmin=paramsABS()[["minpart"]],
			Nmax=paramsABS()[["maxpart"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.abs(delta=1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]]),
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			Nmin=paramsABS()[["minpart"]],
			Nmax=paramsABS()[["maxpart"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.abs(delta=paramsABS()[["maxeffect"]],
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			Nmin=paramsABS()[["minpart"]],
			Nmax=paramsABS()[["maxpart"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsABS()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]])*100,"%"),
					paste0(paramsABS()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")

	power.ddPCR.abs.N(delta=paramsABS()[["mineffect"]],
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.abs.N(delta=1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]]),
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.abs.N(delta=1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]]),
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.abs.N(delta=paramsABS()[["maxeffect"]],
			r=paramsABS()[["replicates"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsABS()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]])*100,"%"),
					paste0(paramsABS()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	
	
	power.ddPCR.abs.r(delta=paramsABS()[["mineffect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			rmin=paramsABS()[["minrep"]],
			rmax=paramsABS()[["maxrep"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.abs.r(delta=1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]]),
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			rmin=paramsABS()[["minrep"]],
			rmax=paramsABS()[["maxrep"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.abs.r(delta=1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]]),
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			rmin=paramsABS()[["minrep"]],
			rmax=paramsABS()[["maxrep"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.abs.r(delta=paramsABS()[["maxeffect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			rmin=paramsABS()[["minrep"]],
			rmax=paramsABS()[["maxrep"]],
			interrep=paramsABS()[["interrep"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsABS()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsABS()[["mineffect"]]+paramsABS()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsABS()[["mineffect"]]+2*paramsABS()[["maxeffect"]])*100,"%"),
					paste0(paramsABS()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	

	power.ddPCR.abs.VB(delta=paramsABS()[["effect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			r=paramsABS()[["minrep"]],
			interrepmin=paramsABS()[["minrepl"]],
			interrepmax=paramsABS()[["maxrepl"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.abs.VB(delta=paramsABS()[["effect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			r=floor(1/3*(2*paramsABS()[["minrep"]]+paramsABS()[["maxrep"]])),
			interrepmin=paramsABS()[["minrepl"]],
			interrepmax=paramsABS()[["maxrepl"]],
			alternative=paramsABS()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.abs.VB(delta=paramsABS()[["effect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			r=ceiling(1/3*(paramsABS()[["minrep"]]+2*paramsABS()[["maxrep"]])),
			interrepmin=paramsABS()[["minrepl"]],
			interrepmax=paramsABS()[["maxrepl"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.abs.VB(delta=paramsABS()[["effect"]],
			NNeg.fr=paramsABS()[["fraction"]],
			Vp=paramsABS()[["Vp"]],
			alpha=paramsABS()[["significance"]],
			N=paramsABS()[["partitions"]],
			r=paramsABS()[["maxrep"]],
			interrepmin=paramsABS()[["minrepl"]],
			interrepmax=paramsABS()[["maxrepl"]],
			alternative=paramsABS()[["alternative"]],
			add=T,lty=4)
	legend.char <- c((paramsABS()[["minrep"]]),
					floor(1/3*(2*paramsABS()[["minrep"]]+paramsABS()[["maxrep"]])),
					ceiling(1/3*(paramsABS()[["minrep"]]+2*paramsABS()[["maxrep"]])),
					(paramsABS()[["maxrep"]]))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Replicates")	
	
})

output$plotCNV <- renderPlot({
	par(mfrow=c(2,2))
	#effect size / number of partitions
	power.ddPCR(delta=paramsCNV()[["mineffect"]],
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			Nmin=paramsCNV()[["minpart"]],
			Nmax=paramsCNV()[["maxpart"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR(delta=1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]]),
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			Nmin=paramsCNV()[["minpart"]],
			Nmax=paramsCNV()[["maxpart"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR(delta=1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]]),
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			Nmin=paramsCNV()[["minpart"]],
			Nmax=paramsCNV()[["maxpart"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=3)
	power.ddPCR(delta=paramsCNV()[["maxeffect"]],
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			Nmin=paramsCNV()[["minpart"]],
			Nmax=paramsCNV()[["maxpart"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsCNV()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(paramsCNV()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")

	power.ddPCR.N(delta=paramsCNV()[["mineffect"]],
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.N(delta=1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]]),
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.N(delta=1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]]),
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.N(delta=paramsCNV()[["maxeffect"]],
			r=paramsCNV()[["replicates"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsCNV()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(paramsCNV()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	
	
	power.ddPCR.r(delta=paramsCNV()[["mineffect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			rmin=paramsCNV()[["minrep"]],
			rmax=paramsCNV()[["maxrep"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.r(delta=1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]]),
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			rmin=paramsCNV()[["minrep"]],
			rmax=paramsCNV()[["maxrep"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.r(delta=1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]]),
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			rmin=paramsCNV()[["minrep"]],
			rmax=paramsCNV()[["maxrep"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.r(delta=paramsCNV()[["maxeffect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			rmin=paramsCNV()[["minrep"]],
			rmax=paramsCNV()[["maxrep"]],
			interrep=paramsCNV()[["interrep"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsCNV()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsCNV()[["mineffect"]]+paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsCNV()[["mineffect"]]+2*paramsCNV()[["maxeffect"]])*100,"%"),
					paste0(paramsCNV()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	

	power.ddPCR.VB(delta=paramsCNV()[["effect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			r=paramsCNV()[["minrep"]],
			interrepmin=paramsCNV()[["minrepl"]],
			interrepmax=paramsCNV()[["maxrepl"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.VB(delta=paramsCNV()[["effect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			r=floor(1/3*(2*paramsCNV()[["minrep"]]+paramsCNV()[["maxrep"]])),
			interrepmin=paramsCNV()[["minrepl"]],
			interrepmax=paramsCNV()[["maxrepl"]],
			alternative=paramsCNV()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.VB(delta=paramsCNV()[["effect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			r=ceiling(1/3*(paramsCNV()[["minrep"]]+2*paramsCNV()[["maxrep"]])),
			interrepmin=paramsCNV()[["minrepl"]],
			interrepmax=paramsCNV()[["maxrepl"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.VB(delta=paramsCNV()[["effect"]],
			NNeg.fr=paramsCNV()[["fraction"]],
			NB=paramsCNV()[["NB"]],
			alpha=paramsCNV()[["significance"]],
			N=paramsCNV()[["partitions"]],
			r=paramsCNV()[["maxrep"]],
			interrepmin=paramsCNV()[["minrepl"]],
			interrepmax=paramsCNV()[["maxrepl"]],
			alternative=paramsCNV()[["alternative"]],
			add=T,lty=4)
	legend.char <- c((paramsCNV()[["minrep"]]),
					floor(1/3*(2*paramsCNV()[["minrep"]]+paramsCNV()[["maxrep"]])),
					ceiling(1/3*(paramsCNV()[["minrep"]]+2*paramsCNV()[["maxrep"]])),
					(paramsCNV()[["maxrep"]]))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Replicates")	
	
})

output$plotREL <- renderPlot({
	par(mfrow=c(2,2))
	#effect size / number of partitions
	power.ddPCR(delta=paramsREL()[["mineffect"]],
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			Nmin=paramsREL()[["minpart"]],
			Nmax=paramsREL()[["maxpart"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR(delta=1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]]),
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			Nmin=paramsREL()[["minpart"]],
			Nmax=paramsREL()[["maxpart"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR(delta=1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]]),
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			Nmin=paramsREL()[["minpart"]],
			Nmax=paramsREL()[["maxpart"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=3)
	power.ddPCR(delta=paramsREL()[["maxeffect"]],
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			Nmin=paramsREL()[["minpart"]],
			Nmax=paramsREL()[["maxpart"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsREL()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]])*100,"%"),
					paste0(paramsREL()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")

	power.ddPCR.N(delta=paramsREL()[["mineffect"]],
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.N(delta=1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]]),
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.N(delta=1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]]),
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.N(delta=paramsREL()[["maxeffect"]],
			r=paramsREL()[["replicates"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsREL()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]])*100,"%"),
					paste0(paramsREL()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	
	
	power.ddPCR.r(delta=paramsREL()[["mineffect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			rmin=paramsREL()[["minrep"]],
			rmax=paramsREL()[["maxrep"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.r(delta=1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]]),
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			rmin=paramsREL()[["minrep"]],
			rmax=paramsREL()[["maxrep"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.r(delta=1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]]),
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			rmin=paramsREL()[["minrep"]],
			rmax=paramsREL()[["maxrep"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.r(delta=paramsREL()[["maxeffect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			rmin=paramsREL()[["minrep"]],
			rmax=paramsREL()[["maxrep"]],
			interrep=paramsREL()[["interrep"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=4)
	legend.char <- c(paste0(paramsREL()[["mineffect"]]*100,"%"),
					paste0(1/3*(2*paramsREL()[["mineffect"]]+paramsREL()[["maxeffect"]])*100,"%"),
					paste0(1/3*(paramsREL()[["mineffect"]]+2*paramsREL()[["maxeffect"]])*100,"%"),
					paste0(paramsREL()[["maxeffect"]]*100,"%"))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Effect size")	

	power.ddPCR.VB(delta=paramsREL()[["effect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			r=paramsREL()[["minrep"]],
			interrepmin=paramsREL()[["minrepl"]],
			interrepmax=paramsREL()[["maxrepl"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1))
	power.ddPCR.VB(delta=paramsREL()[["effect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			r=floor(1/3*(2*paramsREL()[["minrep"]]+paramsREL()[["maxrep"]])),
			interrepmin=paramsREL()[["minrepl"]],
			interrepmax=paramsREL()[["maxrepl"]],
			alternative=paramsREL()[["alternative"]],
			ylim=c(0,1),
			add=T,lty=2)
	power.ddPCR.VB(delta=paramsREL()[["effect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			r=ceiling(1/3*(paramsREL()[["minrep"]]+2*paramsREL()[["maxrep"]])),
			interrepmin=paramsREL()[["minrepl"]],
			interrepmax=paramsREL()[["maxrepl"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=3)
	power.ddPCR.VB(delta=paramsREL()[["effect"]],
			NNeg.fr=paramsREL()[["fraction"]],
			NB=paramsREL()[["NB"]],
			alpha=paramsREL()[["significance"]],
			N=paramsREL()[["partitions"]],
			r=paramsREL()[["maxrep"]],
			interrepmin=paramsREL()[["minrepl"]],
			interrepmax=paramsREL()[["maxrepl"]],
			alternative=paramsREL()[["alternative"]],
			add=T,lty=4)
	legend.char <- c((paramsREL()[["minrep"]]),
					floor(1/3*(2*paramsREL()[["minrep"]]+paramsREL()[["maxrep"]])),
					ceiling(1/3*(paramsREL()[["minrep"]]+2*paramsREL()[["maxrep"]])),
					(paramsREL()[["maxrep"]]))	
	legend("bottomright",legend=legend.char,lty=c(1:4), title="Replicates")	
	
})



#########################
# BETWEEN-REP VARIATION #
#########################

values <- reactiveValues(data=data.frame(Positives = c(rep(0,(numSamples))), Total = c(rep(0,(numSamples)))))

adjust.data <- function(df_){
  	n <- nrow(df_)
  	if(n<input$nsample){
		df_[(n+1):input$nsample,1]<-0
  		df_[(n+1):input$nsample,2]<-0
  		row.names(df_)<-1:nrow(df_)
  	}else{
	  	df_<- df_[1:input$nsample,]
  	}
    df_$Positives<-as.numeric(df_$Positives)
    df_$Total<-as.numeric(df_$Total)	
    return(df_)
}

observe({
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
	}
})

observeEvent(input$nsample,({
	
	if(!is.null(input$hot)){
		values$data <- hot_to_r(input$hot)
	}
  	values$data<-adjust.data(values$data)
})
)


output$hot <- renderRHandsontable({
        rhandsontable(values$data,readOnly=F, useTypes= TRUE) %>%
      		hot_table(highlightCol = TRUE, highlightRow = FALSE) %>%
  			hot_col("Positives", format = "0") %>%
  			hot_col("Total", format = "0")
})   


interrepvarABS <- eventReactive(input$calcInterrepABS,{
		list(data = values$data,Vp=as.numeric(input$interrepVp))
	}
)

output$interrepABS <- renderText({
	data.interrep <- as.data.frame(interrepvarABS()[["data"]])
	betweenrep <- varBetween(data.interrep[,2], data.interrep[,1], interrepvarABS()[["Vp"]])
	if(betweenrep<0){
		paste0("There is no evidence for between-replicate variation.")
	} else {
		paste0("The between-replicate variation is: ", betweenrep, ".")
	}
})




valuesCNV <- reactiveValues(data=data.frame(Positives = c(rep(0,(numSamples*2))), Total = c(rep(0,(numSamples*2))), Target = c(rep(TRUE,numSamples),rep(FALSE,numSamples)), stringsAsFactors=FALSE))

adjust.dataCNV <- function(df_){
  	n <- nrow(df_)
  	if(n<(input$nsampleCNV*2)){
		df_[(n+1):(input$nsampleCNV*2),1]<-0
  		df_[(n+1):(input$nsampleCNV*2),2]<-0
  		df_[1:input$nsampleCNV,3]<-TRUE
  		df_[(input$nsampleCNV+1):(input$nsampleCNV*2),3]<-FALSE
  		row.names(df_)<-1:nrow(df_)
  	}else{
	  	df_<- df_[c(1:(input$nsampleCNV),(n/2+1):(n/2+input$nsampleCNV)),]
  		df_[1:input$nsampleCNV,3]<-TRUE
  		df_[(input$nsampleCNV+1):(input$nsampleCNV*2),3]<-FALSE	  	
   		row.names(df_)<-1:nrow(df_) 		
  	}
    df_$Positives<-as.numeric(df_$Positives)
    df_$Total<-as.numeric(df_$Total)	
    df_$Target<-as.logical(df_$Target)
    return(df_)
}

observe({
	if(!is.null(input$hotCNV)){
		valuesCNV$data <- hot_to_r(input$hotCNV)
	}
})

observeEvent(input$nsampleCNV,({
	
	if(!is.null(input$hotCNV)){
		valuesCNV$data <- hot_to_r(input$hotCNV)
	}
  	valuesCNV$data<-adjust.dataCNV(valuesCNV$data)
})
)


output$hotCNV <- renderRHandsontable({
        rhandsontable(valuesCNV$data,readOnly=F, useTypes= TRUE) %>%
      		hot_table(highlightCol = TRUE, highlightRow = FALSE) %>%
  			hot_col("Positives", format = "0") %>%
  			hot_col("Total", format = "0") %>%
  			hot_col("Target")
})   


interrepvarCNV <- eventReactive(input$calcInterrepCNV,{
		list(data = valuesCNV$data,NB=as.numeric(input$interrepNB))
	}
)

output$interrepCNV <- renderText({
	data.interrep <- as.data.frame(interrepvarCNV()[["data"]])
	betweenrep <- varBetweenCNV(data.interrep, interrepvarCNV()[["NB"]])
	if(betweenrep<0){
		paste0("There is no evidence for between-replicate variation.")
	} else {
		paste0("The between-replicate variation is: ", betweenrep, ".")
	}
})








valuesREL <- reactiveValues(data=data.frame(Positives = c(rep(0,(numSamples*2))), Total = c(rep(0,(numSamples*2))), Target = c(rep(TRUE,numSamples),rep(FALSE,numSamples)), stringsAsFactors=FALSE))

adjust.dataREL <- function(df_){
  	n <- nrow(df_)
  	if(n<(input$nsampleREL*2)){
		df_[(n+1):(input$nsampleREL*2),1]<-0
  		df_[(n+1):(input$nsampleREL*2),2]<-0
  		df_[1:input$nsampleREL,3]<-TRUE
  		df_[(input$nsampleREL+1):(input$nsampleREL*2),3]<-FALSE
  		row.names(df_)<-1:nrow(df_)
  	}else{
	  	df_<- df_[c(1:(input$nsampleREL),(n/2+1):(n/2+input$nsampleREL)),]
  		df_[1:input$nsampleREL,3]<-TRUE
  		df_[(input$nsampleREL+1):(input$nsampleREL*2),3]<-FALSE	  	
   		row.names(df_)<-1:nrow(df_) 		
  	}
    df_$Positives<-as.numeric(df_$Positives)
    df_$Total<-as.numeric(df_$Total)	
    df_$Target<-as.logical(df_$Target)
    return(df_)
}

observe({
	if(!is.null(input$hotREL)){
		valuesREL$data <- hot_to_r(input$hotREL)
	}
})

observeEvent(input$nsampleREL,({
	
	if(!is.null(input$hotREL)){
		valuesREL$data <- hot_to_r(input$hotREL)
	}
  	valuesREL$data<-adjust.dataREL(valuesREL$data)
})
)


output$hotREL <- renderRHandsontable({
        rhandsontable(valuesREL$data,readOnly=F, useTypes= TRUE) %>%
      		hot_table(highlightCol = TRUE, highlightRow = FALSE) %>%
  			hot_col("Positives", format = "0") %>%
  			hot_col("Total", format = "0") %>%
  			hot_col("Target")
})   


interrepvarREL <- eventReactive(input$calcInterrepREL,{
		list(data = valuesREL$data,NB=as.numeric(input$interrepNB))
	}
)

output$interrepREL <- renderText({
	data.interrep <- as.data.frame(interrepvarREL()[["data"]])
	betweenrep <- varBetweenCNV(data.interrep, interrepvarREL()[["NB"]])
	if(betweenrep<0){
		paste0("There is no evidence for between-replicate variation.")
	} else {
		paste0("The between-replicate variation is: ", betweenrep, ".")
	}
})



output$reportGenerateABS <- downloadHandler(
	filename = function() {
		paste(paste("report",strftime(strptime(date(),format="%a %b %d %H:%M:%S %Y"),format="%Y%m%d%H%M%S"), sep=''), sep = '.', 'pdf')
	},
	content = function(file) {
		src <- normalizePath('reportABS.Rmd')
		dat <- paramsABS()
		# temporarily switch to the temp dir, in case you do not have write
		# permission to the current working directory
		owd <- setwd(tempdir())
		on.exit(setwd(owd))
		file.copy(src, 'reportABS.Rmd',overwrite=TRUE)
		library(rmarkdown) 
		outReport <- rmarkdown::render('reportABS.Rmd', pdf_document(),params=list(ls()))
		file.rename(outReport, file)
	}	
)



output$reportGenerateCNV <- downloadHandler(
	filename = function() {
		paste(paste("report",strftime(strptime(date(),format="%a %b %d %H:%M:%S %Y"),format="%Y%m%d%H%M%S"), sep=''), sep = '.', 'pdf')
	},
	content = function(file) {
		src <- normalizePath('reportCNV.Rmd')
		dat <- paramsCNV()
		# temporarily switch to the temp dir, in case you do not have write
		# permission to the current working directory
		owd <- setwd(tempdir())
		on.exit(setwd(owd))
		file.copy(src, 'reportCNV.Rmd',overwrite=TRUE)
		library(rmarkdown) 
		outReport <- rmarkdown::render('reportCNV.Rmd', pdf_document(),params=list(ls()))
		file.rename(outReport, file)
	}	
)

output$reportGenerateREL <- downloadHandler(
	filename = function() {
		paste(paste("report",strftime(strptime(date(),format="%a %b %d %H:%M:%S %Y"),format="%Y%m%d%H%M%S"), sep=''), sep = '.', 'pdf')
	},
	content = function(file) {
		src <- normalizePath('reportREL.Rmd')
		dat <- paramsREL()
		# temporarily switch to the temp dir, in case you do not have write
		# permission to the current working directory
		owd <- setwd(tempdir())
		on.exit(setwd(owd))
		file.copy(src, 'reportREL.Rmd',overwrite=TRUE)
		library(rmarkdown) 
		outReport <- rmarkdown::render('reportREL.Rmd', pdf_document(),params=list(ls()))
		file.rename(outReport, file)
	}	
)

})
