#################################
# Load required libraries
#################################
library(shiny)
library(dplyr)
library(ggplot2)
library(bslib)
library(shinythemes)
library(DT)
library(devtools)
library(pairwiseAdonis)
library(pegas)
library(poppr)
library(vegan)

#################################
# User interface
#################################
ui <- fluidPage(theme = shinytheme('yeti'),
  titlePanel("MDM Calculator"),
    br(),
  "Assuming a stepwise mutation model, mutational distances from the calculated
  modal haplotype (MDM) is defined as the sum of the absolute number of repeat differences
  at each locus that an individual haplotype deviates from the population-specific modal
  haplotype within a particular haplogroup. MDM from the global modal haplotype takes on the same
  assumptions but calculates the modal haplotype from all input haplotypes. Input data must
  have two columns: HG indicating haplogroups and Pop indicating population. See example data for 
  exact data format.",
      br(),
      br(),
  sidebarLayout(
    sidebarPanel(
      fileInput("myfileinput", "Input CSV file", 
      	accept = "text/csv",
      	placeholder = "No file selected"),
      checkboxInput("header", "Header", TRUE),
      selectInput("populationselect", 
      	label = "Population", 
      	choice=character(0)),
      selectInput("haplogroupselect", 
      	label = "Haplogroup", 
      	choice=character(0)),
      downloadButton("downloadPlot", "Download pdf"),
      br(),
      br(),
      a(href="example_data.csv", "Download example data", download=NA, target="_blank"),
    ),

    mainPanel(

    	tabsetPanel(type = "tabs",
    		tabPanel("MDM Plot", 
    			h3(textOutput("caption")),
    			plotOutput("plotview"),
    			dataTableOutput("neighbors")
    			),
    		tabPanel("Modal haplotype", 
    			h3("Modal haplotype"), 
    			dataTableOutput("modeout")
    			),
    		tabPanel("Global MDM",
					h3("MDM from global modal haplotype"),
					plotOutput("popplot")
    			),
    		tabPanel("PCoA", 
    			h3("PCoA"), 
    			plotOutput("pcoa"),
    			dataTableOutput("perma")
    			),
    		)
    )
  )
)
#################################
# Server instance
#################################
server <- function(input, output, session) {
	#Reactive to store loaded data
	reactives <- reactiveValues(       
        mydata = NULL      
    	) 

	#Observe file being selected
	observeEvent(input$myfileinput, {    
		#Store loaded data in reactive
		reactives$mydata <- read.csv(file = input$myfileinput$datapath)
		#Update select input
		updateSelectInput(session, 
			inputId = 'populationselect', 
			label = 'Population', 
			choices  = reactives$mydata$Pop)
		updateSelectInput(session, 
			inputId = 'haplogroupselect', 
			label = 'Haplogroup', 
			choices  = reactives$mydata$HG) 
		updateSelectInput(session, 
			inputId = 'test', 
			label = 'test') 
		})

#################################
# Save Input 
#################################
   mydata <- reactive({
		req(input$myfileinput, input$header, file.exists(input$myfileinput$datapath))
		read.csv(input$myfileinput$datapath, header = input$header)
		})

#################################
# Select HG and Pop From Input 
#################################
	hgdat <- reactive({
	# pull data selected by user from dataframe
		hgsplit <- mydata() %>% 
			filter(Pop == input$populationselect & HG == input$haplogroupselect)
		# remove any columns with missing values
		dat <- hgsplit[,colSums(is.na(hgsplit)) == 0]
		# only include numeric columns for downstream processing
		dat <- select_if(dat, is.numeric)
		})

#################################
# Sample Number in Plot 
#################################
	subxlab <- reactive({
		subxlab <- paste("(n = ", dim(hgdat())[1], ")", sep="")
		})

#################################
# Modal haplotype calculation
#################################
	modehap <- reactive({		
		# function to calculate modal value 
		getmode <- function(v) {
			uniqv <- unique(v)
			uniqv[which.max(tabulate(match(v, uniqv)))]
		}
		# find modal value across all columns
		mode <- sapply(as.data.frame(hgdat()), function(x) getmode(x))
		mode		
		})

#################################
# MDM Calculation
#################################
	mdm <- reactive({
		# empty vector to populate intragroup mutational distance scores  
		x <- c()
		# calculate number of mutations from modal haplotype 
		for(row in 1:nrow(hgdat())){
			temp <- rbind(hgdat()[row,], modehap()) 
			temp <- temp - as.list(temp[1,])
			x <- c(x, sum(abs(temp[2,])))
		}
		# save as dataframe for plotting
		x <- as.data.frame(x)
		})

#################################
# MDM Histogram Function
#################################
	hgplot <- reactive({
		# plot
		hgplot <- ggplot(mdm(), aes(x=x)) + 
			geom_histogram(aes(y=after_stat(density)), 
				binwidth=1, 
				colour="white", 
				fill="grey") + 
			theme_minimal() + 
			xlab(paste("Mutational steps from mode", subxlab())) + 
			ylab("Frequency") + 
			geom_density(bw=1, lty=2)
		hgplot
		})

#################################
# Neighbor haplotypes Function
#################################
		neighbors <- reactive({
			prop.table(table(mdm() <= 1))
			})

#################################
# All Populations Density Plot 
#################################
populationPlot <- reactive({
	# pull data selected by user from dataframe
		hgsplit <- mydata() %>% 
			filter(HG == input$haplogroupselect)
		# remove last column with HG info
		hgclean <- hgsplit[,-ncol(hgsplit)]
		# remove columns with missing values
		dat <- hgclean[,colSums(is.na(hgclean)) == 0]
		# save population vector for plotting
		Population <- dat$Pop
		# only include numeric columns for downstream processing
		dat <- select_if(dat, is.numeric)
		# get number of individuals
		subxlab <- paste("(n = ", dim(dat[1]), ")", sep="")
		# calculate modal haplotype
		getmode <- function(v) {
			uniqv <- unique(v)
			uniqv[which.max(tabulate(match(v, uniqv)))]
		}
		mode <- sapply(as.data.frame(dat), function(x) getmode(x))
		mode
		# empty vector to populate intragroup mutational distance scores  
		x <- c()
		# calculate number of mutations from modal haplotype 
		for(row in 1:nrow(dat)){
			temp <- rbind(dat[row,], mode) 
			temp <- temp - as.list(temp[1,])
			x <- c(x, sum(abs(temp[2,])))
		}
		# save as dataframe for plotting
		x <- as.data.frame(cbind(Population, x))
		ggplot(x, aes(x=as.numeric(x), group=Population, colour=Population)) + 
		geom_density(alpha=0.5) +
		theme_minimal() +
		xlab(paste("Mutational steps from mode", subxlab)) +
		ylab("frequency")
		})

#################################
# PCoA Plot Function
#################################
		pcoa <- reactive({
			# pull data selected by user from dataframe
			hgsplit <- mydata() %>% 
				filter(HG == input$haplogroupselect)
			groups <- df2genind(hgsplit, sep="\t", pop=hgsplit$Pop, ploidy=1)
			# make sure there are no missing data
			groups <- missingno(groups, type="loci")
			# generate genpop object
			obj <- genind2genpop(groups)
			ca1 <- dudi.coa(as.data.frame(obj$tab), scannf=FALSE, nf=3)
			barplot(ca1$eig, main="Eigenvalues")
			s.label(ca1$li, csub = 2, sub="CA 1-2")
			add.scatter.eig(ca1$eig, nf = 3, xax = 1, yax = 2, posi = "bottomright")
			})


#################################
# PERMANOVA
#################################
	perma <- reactive({
		# pull data selected by user from dataframe
		hgsplit <- mydata() %>% 
			filter(HG == input$haplogroupselect)
		haps <- df2genind(hgsplit, sep="\t", pop=hgsplit$Pop, ploidy=1)
		haps <- missingno(haps, type="loci")
		# for now assume all loci are quadnucleotide repeats for testing
		ssr <- rep(4, nLoc(haps))
		# permanova analysis
		bdist <- bruvo.dist(haps, replen=ssr)
		permares <- as.data.frame(pairwise.adonis(bdist, hgsplit$Pop))
		permares
		})

#################################
# Outputs
#################################
	formulaText <- reactive({
		paste(input$populationselect, input$haplogroupselect)
		})
	# Return the formula text for printing as a caption ----
	output$caption <- renderText({
		formulaText()
		})
	# output object to print the modal haplotype
	output$modeout <- renderDataTable({
		temp <- as.data.frame(cbind(rownames(as.data.frame(modehap())), modehap()))
		colnames(temp) <- c("Locus", "Repeat Count")
		temp
		})
	# output object to print neighbor haplotypes
	output$neighbors <- renderDataTable({
		datatable(as.data.frame(neighbors()),
			caption = "Proportion Neighbor Haplotypes (MDM <= 1)",
			rownames = FALSE, 
			colnames = c("Is.neighbor", "Frequency"))
		})

	# Return results from permanova analysis
	output$perma <- renderDataTable({
		datatable(as.data.frame(perma()))
		})
	# output object to view the plot
	output$plotview <- renderPlot({hgplot()})
	output$popplot <- renderPlot({populationPlot()})
	output$pcoa <- renderPlot({pcoa()})
	# download the plot    
 	output$downloadPlot <- downloadHandler(
    filename = function(){
    	paste(input$populationselect, input$haplogroupselect, ".pdf", sep="")
    	},
    content = function(file) {
      pdf(file)
      print(hgplot())
      print(populationPlot())
      dev.off()
    }
  )

}

#################################
# Launch app
#################################
shinyApp(ui, server)