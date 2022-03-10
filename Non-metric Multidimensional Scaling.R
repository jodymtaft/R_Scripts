#From https://www.youtube.com/watch?v=QljEeBei-JA
setwd("C:/Users/TaftJ/Documents/Projects/Psammophis leightoni/Psammophis_submissions")

#Attach data
dataNMDS <- read.delim("Morphology_nmds.txt", header=TRUE)

attach(dataNMDS)
str(dataNMDS)
ncol (dataNMDS) #This calls the number of columns in the dataframe

#Install packages
#install.packages("RColorBrewer")
library("vegan")
library("ggplot2")
library("RColorBrewer")



#########################
#subset the dataframe on which to base the ordination (dataframe 1)
data_1 <- dataNMDS[,3:17]

#Identify the columns that contains the descriptive/environmental data (dataframe 2)
data_2 <- dataNMDS[,1:2]


#ordination by NMDS
NMDS <- metaMDS(data_1, distance = "binomial", k = 3, trymax = 200, autotransform = FALSE)

#########################
#Data visualisation
colors <- c("green", "blue", "yellow")
names(colors) = c("leightoni", "namibensis", "trinasalis")

plot(NMDS$points, 
     type="p",
     xlab = "axis 1", ylab = "axis 2",
     main= "Psammophis morphology - nMDS plot",
     pch=19,
     col=c("green", "blue", "yellow"),
     cex=1.2
     )


#########################
#Data visualisation

#Plot ordination so that points are coloured and shaped according to the groups of interest
species <- as.factor(data_2$Species)

plot(NMDS$points, main="Psammophis morphology", col=species, pch=19,  xlab = "axis 1", ylab = "axis 2")

#Connect the points that belong to the same treatment with ordispider
ordihull(NMDS, groups = data_2$Species,  label = FALSE, draw = c("polygon"), col = colors, alpha = 0.2, lty = "dashed")  #lty = line type, alpha = transparency of shape fill

#Add legend
txt <- c("P. leightoni", "P. namibensis", "P. trinsalis")
# txt <- c("Stress = 0.074") #If you need to add the stress component
legend('bottomright', txt, pch=19, col=c("green", "blue", "yellow"),cex=0.8, bty = "n", text.font = 3)



#####################
#Bootstrapping and testing for differences between the groups
fit <- adonis(data_1 ~ Species, data=data_2, permutations=999, method="binomial")
fit


#####################
#Check assumption of homogeneity of multivariate dispersion
distances_data <- vegdist(data_1)
anova(betadisper(distances_data, data_2$Species))

#Export image to TIFF
tiff(file="nMDS_Plot_MCH.svg",
     width=10, height=6, units="in", res=300)
#After this then replot the graph with the legend codes before running dev.off()
dev.off()
