# faroes_y

## Faroes project mutational distance from modal haplotype

```R
library(dplyr)
library(ggplot2)

hgs <- read.table("faroes_hg_data.txt", header=T, row.names=1)

# split by population and haplogroup
hgsplit <- hgs %>% 
	group_split(Pop,HG)

for(group in 1:length(hgsplit)){
	# pull group, remove any columns with missing values, get label for mutational distance histogram
	dat <- hgsplit[[group]]
	dat <- dat[,colSums(is.na(dat)) == 0]
	title_lab <- paste(dat$Pop[1], dat$HG[1])
	# only include numeric columns for downstream processing
	dat <- select_if(dat, is.numeric)
	print(paste("Calculating mutational distance histograms for: ", title_lab))
	print(paste("Number of individuals: ", dim(dat)[1]))
	subxlab <- paste("(n = ", dim(dat)[1], ")", sep="")
	subxlab

	# function to calculate modal value 
	getmode <- function(v) {
	   uniqv <- unique(v)
	   uniqv[which.max(tabulate(match(v, uniqv)))]
	}
	# find modal value across all columns
	mode <- sapply(as.data.frame(dat), function(x) getmode(x))
  	# empty vector to populate intragroup mutational distance scores  
  	x <- c()

  	# calculate number of mutations from modal haplotype 
  	for(row in 1:nrow(dat)){
		temp <- rbind(dat[row,], mode) 
		temp <- temp - as.list(temp[1,])
		x <- c(x, sum(abs(temp[2,])))
	}
	
	# save as dataframe for plotting
	x <- as.data.frame(x)

	# figure name 
	figout <- paste(title_lab, "_hist.pdf", sep="")
	figout <- gsub(" ", "_", figout)

	# plot mutational distance histogram
	pdf(figout)
	ggplot(x, aes(x=x)) + 
		geom_histogram(aes(y=after_stat(density)), binwidth=1, colour="black", fill="white") + 
		theme_minimal() + 
		xlab(paste("Mutational steps from mode", subxlab)) + 
		ylab("Frequency") + 
		ggtitle(title_lab) +
		geom_density(bw=1, lty=2) 
	dev.off()

}

```

Test shiny app

```R
library(shiny)
# run the application
runApp("/Users/mann/github/faroes_y/mdmHistogram")

# save to server
library(rsconnect)
rsconnect::deployApp("/Users/mann/github/faroes_y/mdmHistogram/")
```

