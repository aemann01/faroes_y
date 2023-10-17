#################################
# Load required libraries
#################################
library(shiny)
library(dplyr)
library(ggplot2)
library(bslib)
#################################
# User interface
#################################
ui <- fluidPage(
  titlePanel("MDM Calculator"),
  "This calculator requires your input to be comma separated with a .csv extension. 
  Input data file needs to have two columns: 
  one with the haplogroup information called HG and another with population information called Pop",
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
      downloadButton("downloadPlot", "Download pdf")
    ),

    mainPanel(
		h3(textOutput("caption")),
		plotOutput("plotview"),
		h3("Population specific modal haplotype"),
		dataTableOutput("modeout"),
		h3("Population distance from global modal haplotype"),
		plotOutput("popplot")

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
		ggplot(x, aes(x=x, group=Population, colour=Population)) + 
		geom_density() +
		theme_minimal() +
		xlab(paste("Mutational steps from mode", subxlab)) +
		ylab("frequency")
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
		colnames(temp) <- c("Locus", "Repeat Number")
		temp
		})
	# output object to view the plot
	output$plotview <- renderPlot({hgplot()})
	output$popplot <- renderPlot({populationPlot()})
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