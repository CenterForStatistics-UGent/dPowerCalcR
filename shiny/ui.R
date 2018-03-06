#install.packages("shiny")
#install.packages("shinydashboard")
#install.packages("rhandsontable")
#install.packages("shinyBS")
#install.packages("rmarkdown")
#install.packages("ggplot2")
library(rhandsontable)
library(shinydashboard)
library(shinyBS)
library(shiny)
library(ggplot2)
library(rmarkdown)

ui<-shinyUI(dashboardPage(
    skin="red",
	dashboardHeader(title = "dPowerCalcR"),
	dashboardSidebar(
	sidebarMenu(
    	menuItem("Power calculation", tabName = "power", icon = icon("calculator"),
    		menuSubItem("Absolute quantification", tabName = "abs"),
   	 		menuSubItem("Relative quantification", tabName = "rel"),
   	 		menuSubItem("Copy number", tabName = "cnv")),
   	 	menuItem("Between-replicate variation", tabName = "interrep", icon = icon("bullseye"),	
    		menuSubItem("Absolute quantification", tabName = "interrepABS"),
   	 		menuSubItem("Relative quantification", tabName = "interrepREL"),
   	 		menuSubItem("Copy number", tabName = "interrepCNV")),   	 	
       	menuItem("Help", tabName = "help", icon = icon("question-circle"),
    		menuSubItem("Citation", tabName = "cite"),
   	 		menuSubItem("FAQ", tabName = "faq"))
    )

),
dashboardBody(
	tabItems(

	################################
	# ABS
	################################

    tabItem(tabName = "abs",
    	h2("Power calculation for absolute quantification"),
		fluidRow(
			box(
				#HTML("<hr>"),
				title="Power parameters",
				width=6,
				height=770,
				textInput("partitions","Number of partitions",value="15000"),
				textInput("replicates","Number of replicates",value="3"),
				textInput("interrep","Between-replicate variation",value="0.001"),
				textInput("effect","Effect size",value="0.1"),
				textInput("Vp","Partition volume",value="1"),
				sliderInput("fraction", "Fraction of negatives", 0, 1, 0.2),
				sliderInput("significance", "Significance level", 0, 0.1, 0.05),
				selectInput("Alt","Hypothesis",choices=c("Two-sided","Less","Greater")),
				HTML("<br>"),
				actionButton("calcResultsABS", "Calculate")				
			),
			box(
				title="Plot parameters",
				width=6,
				height=770,
				textInput("minpart","Minimum number of partitions",value="1000"),
				textInput("maxpart","Maximum number of partitions",value="20000"),
				textInput("minrep","Minimum number of replicates",value="2"),
				textInput("maxrep","Maximum number of replicates",value="10"),
				textInput("mineffect", "Minimal effect size", value="0.05"),
				textInput("maxeffect", "Maximal effect size", value="0.2"),
				textInput("minrepl", "Minimal between-replicate variation", value="0"),
				textInput("maxrepl", "Maximal between-replicate variation", value="0.01")
			)
	   ),
	   fluidRow(
			box(
				title="Results",width=12,
				textOutput("textABS"),
				plotOutput("plotABS",height="1000px")

			)
		),
	   fluidRow(
			box(title="Generate report",width=12,
	        	textInput("reportTitle","Title of the report",value="Digital PCR experiment design report"),
	        	downloadButton("reportGenerateABS","Generate report")
    	)      
      )
	),
	
	
	################################
	# REL QUANT 
	################################
	tabItem(tabName = "rel",
    	h2("Power calculation for relative quantification"),
		fluidRow(
			box(
				#HTML("<hr>"),
				title="Power parameters",
				width=6,
				height=770,
				textInput("partitionsREL","Number of partitions",value="15000"),
				textInput("replicatesREL","Number of replicates",value="3"),
				textInput("interrepREL","Between-replicate variation",value="0.001"),
				textInput("effectREL","Effect size",value="0.1"),
				sliderInput("fractionREL", "Fraction of negatives for target", 0, 1, 0.2),
				sliderInput("significanceREL", "Significance level", 0, 0.1, 0.05),
				#textInput("NBREL","Expected relative quantity",value="1"),
				selectInput("AltREL","Hypothesis",choices=c("Two-sided","Less","Greater")),
				HTML("<br>"),
				actionButton("calcResultsREL", "Calculate")				
			),
			box(
				title="Plot parameters",
				width=6,
				height=770,
				textInput("minpartREL","Minimum number of partitions",value="1000"),
				textInput("maxpartREL","Maximum number of partitions",value="20000"),
				textInput("minrepREL","Minimum number of replicates",value="2"),
				textInput("maxrepREL","Maximum number of replicates",value="10"),
				textInput("mineffectREL", "Minimal effect size", value="0.05"),
				textInput("maxeffectREL", "Maximal effect size", value="0.2"),
				textInput("minreplREL", "Minimal between-replicate variation", value="0"),
				textInput("maxreplREL", "Maximal between-replicate variation", value="0.01")
			)
	   ),
	   fluidRow(
			box(
				title="Results",width=12,
				textOutput("textREL"),
				plotOutput("plotREL",height="1000px")
			)
		),
	   fluidRow(
			box(title="Generate report",width=12,
	        	textInput("reportTitle","Title of the report",value="Digital PCR experiment design report"),
	        	downloadButton("reportGenerateREL","Generate report")
    	)      
      )
	),

	
	################################
	# COPY NUMBER 
	################################
	tabItem(tabName = "cnv",
    	h2("Power calculation for copy number variation"),
		fluidRow(
			box(
				#HTML("<hr>"),
				title="Power parameters",
				width=6,
				height=770,
				textInput("partitionsCNV","Number of partitions",value="15000"),
				textInput("replicatesCNV","Number of replicates",value="3"),
				textInput("interrepCNV","Between-replicate variation",value="0.001"),
				textInput("effectCNV","Effect size",value="0.1"),
				sliderInput("fractionCNV", "Fraction of negatives for target", 0, 1, 0.2),
				sliderInput("significanceCNV", "Significance level", 0, 0.1, 0.05),
				textInput("NB","Expected number of copies",value="2"),
				selectInput("AltCNV","Hypothesis",choices=c("Two-sided","Less","Greater")),
				HTML("<br>"),
				actionButton("calcResultsCNV", "Calculate")				
			),
			box(
				title="Plot parameters",
				width=6,
				height=770,
				textInput("minpartCNV","Minimum number of partitions",value="1000"),
				textInput("maxpartCNV","Maximum number of partitions",value="20000"),
				textInput("minrepCNV","Minimum number of replicates",value="2"),
				textInput("maxrepCNV","Maximum number of replicates",value="10"),
				textInput("mineffectCNV", "Minimal effect size", value="0.05"),
				textInput("maxeffectCNV", "Maximal effect size", value="0.2"),
				textInput("minreplCNV", "Minimal between-replicate variation", value="0"),
				textInput("maxreplCNV", "Maximal between-replicate variation", value="0.01")
			)
	   ),
	   fluidRow(
			box(
				title="Results",width=12,
				textOutput("textCNV"),
				plotOutput("plotCNV",height="1000px")
			)
		),
	   fluidRow(
			box(title="Generate report",width=12,
	        	textInput("reportTitle","Title of the report",value="Digital PCR experiment design report"),
	        	downloadButton("reportGenerateCNV","Generate report")
    	)      
      )
	),	
	
	
	tabItem(tabName = "interrepABS",
		h2("Calculate between-replicate variation for absolute quantities"),
    	fluidRow(
    		box(
		 		sliderInput("nsample", "Number of samples", 2, 100, 5),
		 		textInput("interrepVp","Partition volume", "1"),
				rHandsontableOutput("hot"),
				HTML("<br>"),
				actionButton("calcInterrepABS", "Calculate")		 	      
      	)),
	   fluidRow(
			box(
				textOutput("interrepABS")
			))      	
    ),	
	
	tabItem(tabName = "interrepCNV",
		h2("Calculate between-replicate variation for copy numbers"),
    	fluidRow(
    		box(
		 		sliderInput("nsampleCNV", "Number of samples", 2, 100, 5),
		 		textInput("interrepNB","Expected number of copies", "2"),
				rHandsontableOutput("hotCNV"),
				HTML("<br>"),
				actionButton("calcInterrepCNV", "Calculate")		 	      
      	)),
	   fluidRow(
			box(
				textOutput("interrepCNV")
			))      	
    ),		

	tabItem(tabName = "interrepREL",
		h2("Calculate between-replicate variation for relative quantities"),
    	fluidRow(
    		box(
		 		sliderInput("nsampleREL", "Number of samples", 2, 100, 5),
				rHandsontableOutput("hotREL"),
				HTML("<br>"),
				actionButton("calcInterrepREL", "Calculate")		 	      
      	)),
	   fluidRow(
			box(
				textOutput("interrepREL")
			))      	
    ),	
    	
	tabItem(tabName = "cite",
    	box(

		 	h4("Citation"),
		 	p("When using this application to design your dPCR experiments, please cite:"),
		 	p("Vynck, M. et al. (2018). On the design of digital PCR experiments."),
		 	p("Full text available at ..."),
		 	p("More information on digital PCR data analysis on http://dpcr.ugent.be.")		 	      
      	)
    ),
    tabItem(tabName = "faq",
    	box(
			h4("Where can I find the underlying methodology?"),
			p("The full text of our paper is available at ..."),
		 	h4("Questions? Feel free to contact me at Matthijs.Vynck@UGent.be.")
      	)
	)
  )
)
)
)
