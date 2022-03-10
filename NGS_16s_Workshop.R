#set working directory
setwd("C:\\Users\\Jody\\Desktop\\Mothur.win\\mothur")

#####
#Downloading phyloseq package. Note: You will only need to do the next few commands
#the first time you use the R project.Note : it will take a little while for phyloseq and
#associated packages to download and install. Wait until the "stop sign" is gone and the
#">" sign in the console has appeared before continuing
#For all subsequent analyses, start at the "Load Packages" section.

#if (!requireNamespace ("BiocManager", quietly =TRUE))
#+  install.packages("BiocManager")
#+BiocManager::install(version = "3.13")
#+BiocManager::install(c("phyloseq"))

#Load packages
library(phyloseq)

library(lubridate)
library(ggplot2)
library(doBy)
library(reshape2)
library(dplyr)
library(vegan)
library(SPECIES)
library(phangorn)
library(GUniFrac)

theme_set(theme_bw())

#######import file#########

sharedfile = "3rivers.unique.pick.good.filter.unique.pick.an.shared"
taxfile = "3rivers.unique.pick.good.filter.unique.pick.an.0.03.cons.taxonomy"
mapfile = "RiverWater.csv"

mothur_data<-import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)

map <- read.csv(mapfile)
head(map)

map <- sample_data(map)
rownames(map) <- map$SampleID

moth_merge <- merge_phyloseq(mothur_data, map)
moth_merge

colnames(tax_table(moth_merge))
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

moth_sub <- moth_merge
colnames(tax_table(moth_sub))

Mydata <- moth_sub
Mydata_scale <-Mydata

### You can draw bar graphs but to be honest, the ones you can draw in Excel look much nicer

plot_bar(Mydata, fill = "Phylum")

plot_bar(Mydata, fill = "Class")

#lets look at a subset of the data 
#first lets look at just the Proteobacteria (and rename this subset of data "SUBS")

SUBS <- subset_taxa(Mydata, Phylum=="Proteobacteria")
plot_bar(SUBS, fill = "Class")

#Lets go one step futher and look at just the Alphaproteobacteria
#(name this subset as "SUBS2")

SUBS2 <- subset_taxa(Mydata, Class=="Alphaproteobacteria")
plot_bar(SUBS2, fill = "Order")
plot_bar(SUBS2, fill = "Genus", facet_grid=~Order)

#.........HEATMAPs............#

plot_heatmap(Mydata)

#You can look at a subset of your data. In this example, we are looking at the 
#Rhodobacteraceae and have renamed this sub-dataset as "SUBS3"

SUBS3 <- subset_taxa(Mydata, Class=="Chlamydiae")
plot_heatmap(SUBS3)

plot_heatmap(SUBS3, taxa.label = "Genus")

plot_heatmap(SUBS3, taxa.label = "Family", low="orange", high = "firebrick2")

plot_heatmap(SUBS3, taxa.label = "Family", low="orange", high="firebrick2", na.value="white")

#background is grey 

plot_heatmap(SUBS3, taxa.label = "Family", low = "orange", high = "firebrick2", na.value = "grey", title = "MY HEATMAP IS ON FIRE")

#...Alpha Diversity.....#

#Looking at all the indices: 
#"observed", "Chaol", "ACE", "Shannon", "InvSimpson", "Simpson", "Fisher"):

plot_richness(Mydata, measures = NULL)

#If you just want selected indices, then specify the ones you want:
plot_richness(Mydata, measures = c("observed","Chao1", "Shannon", "Simpson"))

#Subsampling your dataset
sample_sums(Mydata)
datasubsampled <- rarefy_even_depth(Mydata, rngseed = TRUE)
sample_sums(datasubsampled)
plot_richness(datasubsampled, measures = c("observed","Chao1", "Shannon", "Simpson"))

#Changing plot colours and symbols
plot1 = plot_richness(datasubsampled, x="River", col="River", measures=c("observed","Chao1", "Simpson"))
p1 <- plot1 + geom_point(aes(shape="location"), size=3) + scale_shape_manual(values = c(8,15,17), name = "Site")
p1

#....Adonis.....#
set.seed(1)

#Calculate the bray-curtis distance matrix and call it "DM" (or anything you like)
DM <- phyloseq::distance(Mydata_scale, method = "bray")
#If you want to see other types of distance matrices are supported i.e. alternatives to Bray-Curtis, use:
#?phyloseq::distance

#make a data frame from the sample data 
sampleDF <- data.frame(sample_data(Mydata))

#Adonis test
adonis(DM ~ River, data = sampleDF) #Refer to workshop present 2020 Part 2.pdf page 55 for interpretation

#....Principle Co-ordinate Analysis (PCoA)....#
MyPCoA <- ordinate(physeq = Mydata_scale, method = "PCoA", distance = "bray")

plot_ordination(physeq = Mydata_scale, ordination = MyPCoA, color = "River", shape = "Location", title = "Anything goes") +
  scale_color_manual(values = c("blue4", "red", "green4", "magenta", "#1919ff", "yellow")) + #Different colours
  geom_point(aes(color = River), size = 2) +
  scale_shape_manual(values = c(8,1,15)) #Different shapes


#......NMDS......#
set.seed(1)

#ordinate
MyNMDS <- ordinate(physeq = Mydata_scale, method = "NMDS", distance = "bray")

#Draw stress plot. You want the R^2 to be as close to 1 as possible
stressplot(MyNMDS)

plot_ordination(physeq = Mydata_scale, ordination = MyNMDS, color = "River", shape = "Location", title = "Anything goes") +
  scale_color_manual(values = c("blue4", "red", "green4", "magenta", "#1919ff", "yellow")) + #Different colours
  geom_point(aes(color = River), size = 4) +
  scale_shape_manual(values = c(8,1,15)) #Different shapes

