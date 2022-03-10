install.packages("devtools")
library(devtools)

#Next install ENMtools, if you downloaded the file from the website, use install.local, otherwise use install.packages

#install_local("C:\\Users\\Jody\\Downloads\\ENMTools-master\\ENMTools-master") #You need devtools loaded before running this line. Download ENMTools-master from github repository: https://github.com/danlwarren/ENMTools
#install.packages("ENMtools") #You may need an older version of R to do this. 
library(ENMTools)

#............***Run ENMTools***...........

env.files <- list.files(path = "C:/Users/TaftJ/Desktop/WLT/EVs/asc", pattern = "asc$", full.names = TRUE) #File for .asc BioClim variables
env <- stack(env.files)
names(env) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
#The names have to be identical to your variables in the directory, any mistake and your script won't run. 

env <- setMinMax(env) #I'm not sure what this does 

raster.cor.matrix(env) #Runs corrrelation matrix on input environmental files - might take a while

raster.cor.plot(env) #Plots two graphs to depict correlation

#Cluster plot that does multidimensional scaling of the predictor variables and then plots them in 
#a two dimensional space so that more correlated predictors are closer to each other.
#Heatmap that colors pairs of predictors by their Pearson correlation coefficient.

env <- env[[c("bio.01", "bio.03", "bio.05", "bio.06", "bio.08", "bio.10", "bio.12", "bio.17", "bio.19", "TRI.")]] #Call variables you select that aren't correlated. These you select after reading the output heatmap and graph
plot(env)

raster.cor.matrix(env) #Repeat process on selected variables 

raster.cor.plot(env) #Redo script and plot graphs until heatmap is completely blue - no correlation between variables



#....................Raster overlap.......................

aridus = raster("C:/Users/TaftJ/Desktop/minor asc/aridus_LIG.asc")    #Import asc file - change file each time for time slice
cloetei = raster("C:/Users/TaftJ/Desktop/minor asc/cloetei_LIG.asc")  #Import asc file - change file each time for time slice
minor = raster("C:/Users/TaftJ/Desktop/minor asc/minor_LIG.asc")  #Import asc file - change file each time for time slice

raster.overlap(aridus, minor, verbose = FALSE) #raster.overlap calculates overlap between two layers 

raster.overlap(aridus, cloetei, verbose = FALSE)

raster.overlap(cloetei, minor, verbose = FALSE)
