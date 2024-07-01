# faroes_y

Shiny app: https://aemann01.shinyapps.io/mdmhistogram/

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

Diversity estimates

```R
library(poppr)
# reload data to make running this block of code by itself easier 
dat <- read.table("faroes_hg_data.FULL.txt", sep="\t", header=T, row.names=1)
sub <- dat[, apply(dat, 2, function(x) !any(is.na(x)))]
# first across all haplogroups
fordiv <- sub[,1:dim(sub)[[2]]-1]
genind <- df2genind(fordiv, sep="\t", pop=fordiv$Pop, ploidy=1)
# get full genotypic diversity across all haplogroups
fullgenodiv <- poppr(genind, minsamp = 10)
# also get genotypic diversity (simpson's index) after controlling for sample size
N <- fullgenodiv$N
lambda <- fullgenodiv$lambda
round((N/(N-1))*lambda, 2)
# [1] 0.97 1.00 0.98 0.98 0.99 0.99 1.00

# by haplogroup
# R1b
x <- sub[sub$HG == "R1b",]
fordiv <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fordiv, sep="\t", pop=fordiv$Pop, ploidy=1)
fullgenodiv <- poppr(genind)
N <- fullgenodiv$N
lambda <- fullgenodiv$lambda
round((N/(N-1))*lambda, 2)
# [1] 0.94 0.99 0.98 0.97 0.98 1.00 0.99

# R1a
x <- sub[sub$HG == "R1a",]
fordiv <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fordiv, sep="\t", pop=fordiv$Pop, ploidy=1)
fullgenodiv <- poppr(genind)
N <- fullgenodiv$N
lambda <- fullgenodiv$lambda
round((N/(N-1))*lambda, 2)
# [1] 0.89 0.99 1.00 0.95 0.98 1.00 0.98

# I1
x <- sub[sub$HG == "I1",]
fordiv <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fordiv, sep="\t", pop=fordiv$Pop, ploidy=1)
fullgenodiv <- poppr(genind)
N <- fullgenodiv$N
lambda <- fullgenodiv$lambda
round((N/(N-1))*lambda, 2)
# [1] 0.86 0.97 0.93 0.83 0.92 0.96 0.98
```

Poppr networks

```R
library(poppr)
# reload data to make running this block of code by itself easier 
dat <- read.table("faroes_hg_data.FULL.txt", sep="\t", header=T, row.names=1)
sub <- dat[, apply(dat, 2, function(x) !any(is.na(x)))]

### R1b haplotype network
x <- sub[sub$HG == "R1b",]

fornet <- x[,1:dim(x)[[2]]-1]
genind <- df2genind(fornet, sep="\t", pop=fornet$Pop, ploidy=1)
# using dissimilarity distance metric
genind_sub <- popsub(genind, exclude = character(0))
genind_dist <- diss.dist(genind_sub, percent = FALSE, mat = FALSE)

# amova results
Pops <- genind_sub$pop
amova.res <- poppr.amova(genind_sub, ~Pops, dist=genind_dist, method="pegas", within=FALSE)
# $call
# ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

# $results
#                  Df   Sum Sq   Mean Sq
# Between samples   5 193.3712 38.674236
# Within samples  315 594.1491  1.886188
# Total           320 787.5202  2.461001

# $componentsofcovariance
#                                 Sigma         %
# Variations  Between samples 0.7758846  29.14589
# Variations  Within samples  1.8861875  70.85411
# Total variations            2.6620721 100.00000

# $statphi
#                         Phi
# Phi-samples-total 0.2914589
amova.test <- randtest(amova.res)
# Monte-Carlo test
# Call: as.randtest(sim = res, obs = sigma[1])

# Observation: 0.7758846

# Based on 99 replicates
# Simulated p-value: 0.01
# Alternative hypothesis: greater

#       Std.Obs   Expectation      Variance
#  8.776340e+01 -4.489916e-04  7.824739e-05

# get mean distance for whole dataset
mean(diss.dist(genind_sub))
# [1] 4.922002

# now for each population separately
# order
popNames(genind_sub)
# [1] "Faroes"  "Denmark" "Ireland" "Iceland" "Norway"  "Sweden"
mean(diss.dist(popsub(genind_sub,1)))
# [1] 3.541176
mean(diss.dist(popsub(genind_sub,2)))
# [1] 4.00341
mean(diss.dist(popsub(genind_sub,3)))
# [1] 3.758982
mean(diss.dist(popsub(genind_sub,4)))
# [1] 3.79414
mean(diss.dist(popsub(genind_sub,5)))
# [1] 3.038095
mean(diss.dist(popsub(genind_sub,6)))
# [1] 4.4


library(pegas)




min_span_net <- poppr.msn(genind_sub, genind_dist, showplot = FALSE, include.ties = TRUE)

set.seed(69)
pdf("R1b_network.pdf")
plot_poppr_msn(genind,
               min_span_net,
               inds = "NA",
               mlg = FALSE,
               gadj = 3,
               nodescale = 10,
               wscale = TRUE,
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

# amova results
amova.res <- poppr.amova(genind_sub, ~Pop, dist=genind_dist)
# $call
# ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

# $results
#                  Df    Sum Sq   Mean Sq
# Between samples   5  94.50873 18.901747
# Within samples  137 219.87588  1.604933
# Total           142 314.38462  2.213976

# $componentsofcovariance
#                                 Sigma         %
# Variations  Between samples 0.8391384  34.33362
# Variations  Within samples  1.6049334  65.66638
# Total variations            2.4440718 100.00000

# $statphi
#                         Phi
# Phi-samples-total 0.3433362
amova.test <- randtest(amova.res)
# Monte-Carlo test
# Call: as.randtest(sim = res, obs = sigma[1])

# Observation: 0.8391384

# Based on 99 replicates
# Simulated p-value: 0.01
# Alternative hypothesis: greater

#      Std.Obs  Expectation     Variance
# 38.320767792  0.003193868  0.000475868
mean(diss.dist(genind_sub))
# [1] 4.427952
# now for each population separately
# order
popNames(genind_sub)
# [1] "Faroes"  "Denmark" "Ireland" "Iceland" "Norway"  "Sweden"
mean(diss.dist(popsub(genind_sub,1)))
# [1] 2.547489
mean(diss.dist(popsub(genind_sub,2)))
# [1] 4.043011
# mean(diss.dist(popsub(genind_sub,3))) # too few observations for ireland
mean(diss.dist(popsub(genind_sub,4)))
# [1] 3.082353
mean(diss.dist(popsub(genind_sub,5)))
# [1] 4.345455
mean(diss.dist(popsub(genind_sub,6)))
# [1] 4.6

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

# amova results
amova.res <- poppr.amova(genind_sub, ~Pop, dist=genind_dist)
# ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

# $results
#                  Df   Sum Sq   Mean Sq
# Between samples   5 155.3242 31.064845
# Within samples  175 238.7531  1.364304
# Total           180 394.0773  2.189319

# $componentsofcovariance
#                                Sigma         %
# Variations  Between samples 1.073185  44.02831
# Variations  Within samples  1.364304  55.97169
# Total variations            2.437489 100.00000

# $statphi
#                         Phi
# Phi-samples-total 0.4402831
amova.test <- randtest(amova.res)
# Monte-Carlo test
# Call: as.randtest(sim = res, obs = sigma[1])

# Observation: 1.073185

# Based on 99 replicates
# Simulated p-value: 0.01
# Alternative hypothesis: greater

#       Std.Obs   Expectation      Variance
# 63.0740307127 -0.0021611805  0.0002906668
mean(diss.dist(genind_sub))
# [1] 4.378637
# now for each population separately
# order
popNames(genind_sub)
# [1] "Faroes"  "Denmark" "Ireland" "Iceland" "Norway"  "Sweden"
mean(diss.dist(popsub(genind_sub,1)))
# [1] 1.830049
mean(diss.dist(popsub(genind_sub,2)))
# [1] 3.415546
mean(diss.dist(popsub(genind_sub,3))) 
# [1] 3.133333
mean(diss.dist(popsub(genind_sub,4)))
# [1] 2.180995
mean(diss.dist(popsub(genind_sub,5)))
# [1] 2.794872
mean(diss.dist(popsub(genind_sub,6)))
# [1] 3.225108

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

Faroes city haplogroup frequency figures

```R
dat <- read.table("faroes_city_hg_freq.txt", sep="\t", header=T)

pdf("hg_klaksvik.pdf")
pie(dat$Klaksvik, labels=dat$HG)
dev.off()

pdf("hg_torshavn.pdf")
pie(dat$Torshavn, labels=dat$HG)
dev.off()

pdf("hg_tvoroyri.pdf")
pie(dat$Tvoroyri, labels=dat$HG)
dev.off()
```










