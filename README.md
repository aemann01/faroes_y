# faroes_y

Shiny app: https://aemann01.shinyapps.io/mdmhistogram/

Under construction! :🚧:

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

Make radial graphs for haplogroup frequencies in each population

```R
install.packages("fmsb")
install.packages("tidyverse")
library(tidyverse)
library(fmsb)
library(ggplot2)

hgs <- c("R1a", "R1b", "I1")

dat <- read.table("faroes_hg_data.txt", sep="\t", header=T, row.names=1)

# radar charts of haplogroup frequencies

for (i in 1:length(unique(dat$Pop))){
	temp <- dat[dat$Pop == unique(dat$Pop)[i],]
	grouped <- temp %>% group_by(HG) %>% summarise(count = n())
	maximum <- max(grouped$count)
	# now collapse other haplogroups into single category
	other <- grouped[!(grouped$HG %in% hgs),]
	x <- cbind("Other", sum(other$count))
	x <- as.data.frame(x)
	colnames(x) <- c("HG", "count")
	x <- rbind(grouped[grouped$HG %in% hgs,], x)
	df <- as.data.frame(rbind(x$HG, x$count))
	colnames(df) <- df[1,]
	df <- df[-1,]
	df <- rbind(rep(maximum, length(x$HG)), rep(0, length(x$HG)), df)
	df2 <- lapply(df, as.numeric)
	pdf(paste(unique(dat$Pop)[i], "_radarHG.pdf", sep=""))
	radarchart(as.data.frame(df2), pfcol=rgb(0.1,0.1,0.1,0.4) , plwd=4)
	dev.off()
}

# classic pie charts of haplogroup frequencies

for (i in 1:length(unique(dat$Pop))){
	temp <- dat[dat$Pop == unique(dat$Pop)[i],]
	grouped <- temp %>% group_by(HG) %>% summarise(count = n())
	maximum <- max(grouped$count)
	# now collapse other haplogroups into single category
	other <- grouped[!(grouped$HG %in% hgs),]
	x <- cbind("Other", sum(other$count))
	x <- as.data.frame(x)
	colnames(x) <- c("HG", "count")
	x <- rbind(grouped[grouped$HG %in% hgs,], x)
	pdf(paste(unique(dat$Pop)[i], "_pieHG.pdf", sep=""))
	pie(as.numeric(x$count), labels=x$HG)
	dev.off()
}
```

Upset plot of shared haplotypes within major haplogroups

```R
install.packages("UpSetR")
library(UpSetR)

# only consider loci that are found in all populations
sub <- dat[, apply(dat, 2, function(x) !any(is.na(x)))]

# this loop is broken somehow, keep for now but run sep
for (i in 1:length(hgs)){
	df <- sub[sub$HG == hgs[[i]],]
	# drop haplogroup column
	df <- df[1:length(df)-1]
	# get unique haplotypes per population
	df <- unique(df)
	# haplotypes as a single column
	haps <- as.data.frame(cbind(df$Pop, unite(col="haptype", df[,-1], sep="_")))
	haps <- tibble(haps)
	# create presence absence pivot table
	df2 <- haps %>% 
	  mutate(present = 1) %>% # create a dummy column
	  pivot_wider(names_from = `df$Pop`, values_from = present) %>% 
	  mutate_at(vars(Faroes, Iceland, Denmark, Ireland, Norway, Sweden), ~ifelse(is.na(.), 0, 1))
	# format for upset plot
	forup <- data.frame(df2[,-1])
	pdf(paste(hgs[[i]], "_upset.pdf", sep=""))
	upset(forup, order.by="freq")
	dev.off()
}

##### R1a
df <- sub[sub$HG == hgs[[1]],]
# drop haplogroup column
df <- df[1:length(df)-1]
# get unique haplotypes per population
df <- unique(df)
# haplotypes as a single column
haps <- as.data.frame(cbind(df$Pop, unite(col="haptype", df[,-1], sep="_")))
haps <- tibble(haps)
# create presence absence pivot table
df2 <- haps %>% 
  mutate(present = 1) %>% # create a dummy column
  pivot_wider(names_from = `df$Pop`, values_from = present) %>% 
  mutate_at(vars(Faroes, Iceland, Denmark, Ireland, Norway, Sweden), ~ifelse(is.na(.), 0, 1))
# reorder
df2 <- df2[,c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes")]
# format for upset plot
forup <- data.frame(df2)
pdf(paste(hgs[[1]], "_upset.pdf", sep=""))
upset(forup, keep.order=T, sets=c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes"))
dev.off()

##### R1b
df <- sub[sub$HG == hgs[[2]],]
# drop haplogroup column
df <- df[1:length(df)-1]
# get unique haplotypes per population
df <- unique(df)
# haplotypes as a single column
haps <- as.data.frame(cbind(df$Pop, unite(col="haptype", df[,-1], sep="_")))
haps <- tibble(haps)
# create presence absence pivot table
df2 <- haps %>% 
  mutate(present = 1) %>% # create a dummy column
  pivot_wider(names_from = `df$Pop`, values_from = present) %>% 
  mutate_at(vars(Faroes, Iceland, Denmark, Ireland, Norway, Sweden), ~ifelse(is.na(.), 0, 1))
# reorder
df2 <- df2[,c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes")]
# format for upset plot
forup <- data.frame(df2)
pdf(paste(hgs[[2]], "_upset.pdf", sep=""))
upset(forup, keep.order=T, sets=c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes"))
dev.off()

##### I1
df <- sub[sub$HG == hgs[[3]],]
# drop haplogroup column
df <- df[1:length(df)-1]
# get unique haplotypes per population
df <- unique(df)
# haplotypes as a single column
haps <- as.data.frame(cbind(df$Pop, unite(col="haptype", df[,-1], sep="_")))
haps <- tibble(haps)
# create presence absence pivot table
df2 <- haps %>% 
  mutate(present = 1) %>% # create a dummy column
  pivot_wider(names_from = `df$Pop`, values_from = present) %>% 
  mutate_at(vars(Faroes, Iceland, Denmark, Ireland, Norway, Sweden), ~ifelse(is.na(.), 0, 1))
# reorder
df2 <- df2[,c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes")]
# format for upset plot
forup <- data.frame(df2)
pdf(paste(hgs[[3]], "_upset.pdf", sep=""))
upset(forup, keep.order=T, sets=c("Denmark", "Sweden", "Norway", "Ireland", "Iceland", "Faroes"))
dev.off()
```

Global population correspondance analysis 

```R
library(vegan)
library(adegenet)

fordist <- sub[,1:dim(sub)[[2]]-1]
groups <- df2genind(fordist, sep="\t", pop=fordist$Pop, ploidy=1)
# make sure there are no missing data
groups <- missingno(groups, type="loci")
# generate genpop object
obj <- genind2genpop(groups)
ca1 <- dudi.coa(as.data.frame(obj$tab), scannf=FALSE, nf=3)
barplot(ca1$eig, main="Eigenvalues")
# display typologies
pdf("correspondance_analysis.pdf")
s.label(ca1$li, csub = 2, sub="CA 1-2")
add.scatter.eig(ca1$eig, nf = 3, xax = 1, yax = 2, posi = "bottomright")
dev.off()
```

Poppr networks

#### TO DO: Why are there no shared haplotypes in these networks?

```R
### R1b haplotype network
x <- sub[sub$HG == "R1b",]
fornet <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fornet, sep="\t", pop=fornet$Pop, ploidy=1)
# using dissimilarity distance metri
genind_sub <- popsub(genind, exclude = character(0))
genind_dist <- diss.dist(genind_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(genind_sub, genind_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
pdf("R1b_network.pdf")
plot_poppr_msn(genind,
               min_span_net,
               inds = "NA",
               mlg = FALSE,
               gadj = 3,
               nodescale = 10,
               palette = funky,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)
dev.off()

### R1a haplotype network
x <- sub[sub$HG == "R1a",]
fornet <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fornet, sep="\t", pop=fornet$Pop, ploidy=1)
# using dissimilarity distance metri
genind_sub <- popsub(genind, exclude = character(0))
genind_dist <- diss.dist(genind_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(genind_sub, genind_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
pdf("R1a_network.pdf")
plot_poppr_msn(genind,
               min_span_net,
               inds = "NA",
               mlg = FALSE,
               gadj = 3,
               nodescale = 10,
               palette = funky,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)
dev.off()

### I1 haplotype network
x <- sub[sub$HG == "I1",]
fornet <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fornet, sep="\t", pop=fornet$Pop, ploidy=1)
# using dissimilarity distance metri
genind_sub <- popsub(genind, exclude = character(0))
genind_dist <- diss.dist(genind_sub, percent = FALSE, mat = FALSE)
min_span_net <- poppr.msn(genind_sub, genind_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
pdf("I1_network.pdf")
plot_poppr_msn(genind,
               min_span_net,
               inds = "NA",
               mlg = FALSE,
               gadj = 3,
               nodescale = 10,
               palette = funky,
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = TRUE,
               layfun = igraph::layout_nicely)
dev.off()
```
