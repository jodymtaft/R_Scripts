#You need the packages dismo, ENMeval, sp and rJava. The maxent.jar file from the GUI must be copied into the 'jarva' folder in the dismo package
#install.packages("dismo")
#install.packages("ENMeval")
#install.packages("sp")
#install.packages("rJava")
#install.packages("spThin") #You'll need this to thin occurrence points by distance     

setwd("C:/Users/TaftJ/Documents/Projects/Cordylus minor/Cordylus_2020/Model Outputs - Taxon/aridus/") #Place to save output files

#..............Reduce spatial auto-correlation in occurrence points...............................

library(spThin)

loc.data <- read.csv( "C:/Users/TaftJ/Desktop/WLT/Sclerophrys_2020.csv") #Make sure the correct RAW occurrence file is selected
                                                                               #A data.frame of occurence locations. It can include several columnns, 
                                                                               #but must include at minimum a column of latitude values, a column of longitude values, and a column of species names.
thin (loc.data , 
      lat.col = "ddlat",          #Name of column of latitude values. Caps sensitive.
      long.col = "ddlong",        #Name of column of longitude values. Caps sensitive.
      spec.col = "Species",        #Name of column of species name. Caps sensitive.
      thin.par = 5,             #Thinning parameter - the distance (in kilometers) that you want records to be separated by
      reps = 1,                 #The number of times to repete the thinning process. Given the random process of removing nearest-neighbors there should be 'rep' number of different sets of coordinates.
      locs.thinned.list.return = FALSE, #TRUE/FALSE - If true, the 'list' of the data.frame of thinned locs resulting from each replication is returned 
      write.files = TRUE,        #TRUE/FALSE - If true, new *.csv files will be written with the thinned locs data
      max.files = 1,             #The maximum number of *csv files to be written based on the thinned data
      out.dir = "C:/Users/TaftJ/Desktop/WLT", #Directory to write new *csv files to
      out.base = "Sclerophrys_2020_5km", #A file basename to give to the thinned datasets created - Change based on km distance
      write.log.file = FALSE,     #TRUE/FALSE create/append log file of thinning run
      log.file = "spatial_thin_log.txt", #Text log file
      verbose = TRUE             #TRUE/FALSE - If true, running details of the function are print at the console
)

library(dismo)
library(ENMeval) #Install ENMeval
library(sp)
library(rJava)
library(rgdal)
library(ecospat)
library(sf)
library(blockCV)
library(dplyr)
library(raster)
library(rangeModelMetadata)


#.......................Set Environmental variables............................................

#files=c()#Here you put your .asc environmental files

path = "C:/Users/TaftJ/Documents/Projects/Cordylus minor/Cordylus_2020/Paleo/WorldClim 1.4 Present 30_ZAF/" #Path to read environmental data
files <- list.files(path, pattern='asc$', full.names=TRUE ) #This line creates a loop to read in all files within 'path' folder
names(files) <- c("bio 1", "bio 2", "bio 3", "bio 4", "bio 5", "bio 6", "bio 7", "bio 8", "bio 9", "bio 10", "bio 11", "bio 12", "bio13", "bio 14", "bio 15", "bio 16", "bio 17", "bio 18", "bio 19", "TRI")
env <- files[c("bio 1", "bio 2", "bio 3", "bio 4", "bio 5", "bio 6", "bio 7", "bio 8", "bio 9", "bio 10", "bio 11", "bio 12", "bio13", "bio 14", "bio 15", "bio 16", "bio 17", "bio 18", "bio 19", "TRI")] #Call variables you select that aren't correlated

raster_stack=stack(env)#stack the environmental files
plot(raster_stack) #Check to see if all layers are imported

#.....................Create Background points................................................

#Now create background points layer either by random with code below (or in vignette) or create a Bias file - see Petford et al 2019 - DOI = 10.1080/21564574.2019.1681524

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: If the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background.
bg <- randomPoints(raster_stack[["bio.1"]], n=10000) #Make sure the raster_stack value is present in environmental 'path' folder e.g. raster_stack[[1]]=bio01.asc
bg <- as.data.frame(bg)

# Check the coverage of background points.
plot(raster_stack[["bio.1"]], legend=FALSE)  #Plot environmental layer used
points(bg, col='blue')                       #Overlay background points onto layer used to check coverage

#Bias<-read.csv()#From Melissa - "I created a bias file and used that as my background points". 

#.....................Read in Occurrence data file............................................

Species<-read.csv("C:/Users/TaftJ/Documents/Projects/Cordylus minor/Cordylus_2020/Maxent/minor_ENMeval_Occ.csv")#Your species occurence data. In the csv file there must be no species name, only coordinates, unlike the Maxent GUI

#.....................Setup Model and Model parameters........................................

eval1 <- ENMevaluate(Species, raster_stack, bg.coords = bg, #species occurence data, raster_stack of environmental variables and your background coords
occ.grp = NULL, #you dont need to worry about this unless you want to choose your own partitioning method
bg.grp = NULL, #As above
RMvalues = seq(0, 4, 0.5), #Valuse to use for regularization paramater, what I have put basically says use between 0 and 4.5 in 0.5 increments
fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), #the feature class combinations to use
#categoricals = c("TRI."), #if you have any categorical environemental variables in the raster_stack include them here 
n.bg = 10000, #Number of background points
method = "jackknife", #Method of data partitioning. Jackknife is appropriate to use when you have 25 or less occurence records. Block method is good when the goal is to project data to a different time period. See vignette and paper for more details.
overlap = FALSE, #If you want to calculate niche overlap between the different results write TRUE
#kfolds = NA, bin.output = TRUE, clamp = TRUE,
parallel = FALSE, progbar = FALSE, updateProgress = FALSE, rasterPreds = NULL, algorithm = "maxent.jar")  #Run in parrellel speeds up computation - see vignette

e.mx.l <- ENMevaluate(occs = Species, envs = raster_stack, bg = bg, 
                      algorithm = 'maxnet', partitions = 'jackknife', 
                      tune.args = list(fc = c("L","LQ","LQH","H"), rm = seq(1, 4 ,0.5)))

e.mx <- ENMevaluate(occs = Species, envs = raster_stack, bg = bg, 
                      algorithm = 'maxent.jar', partitions = 'jackknife', 
                      tune.args = list(fc = c("L","LQ","LQH","H"), rm = seq(1, 4 ,0.5)))


eval1@results #Lists all results of the ENMeval
#write.csv(eval1@results, "All results ENMeval.csv" ) #Writes All results to ENMeval - You dont need this though
e.mx.l@results #Lists all results of the ENMeval
#write.csv(eval1@results, "All results ENMeval.csv" ) #Writes All results to ENMeval - You dont need this though

eval1@results[which(eval1@results$delta.AICc==0),] #returns the model with the smallest change in AIC and thus what is considered the optimum model. Also check the other metrics here, particularly AUC,or10, ormtp. Details in papers.
write.csv(eval1@results[which(eval1@results$delta.AICc==0),], "Optimal Model_aridus_Present_30.csv" ) #Writes output to CSV

e.mx.l@results[which(e.mx.l@results$delta.AICc==0),] #returns the model with the smallest change in AIC and thus what is considered the optimum model. Also check the other metrics here, particularly AUC,or10, ormtp. Details in papers.
write.csv(e.mx.l@results[which(e.mx.l@results$delta.AICc==0),], "Optimal Model_minor_Present_30.csv" ) #Writes output to CSV

res <- eval.results(e.mx.l)

opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

opt.seq <- res %>% 
   filter(or.10p.avg == min(or.10p.avg)) %>% 
   filter(auc.val.avg == max(auc.val.avg))
opt.seq

mod.seq <- eval.models(e.mx.l)[[opt.seq$tune.args]]
mod.seq$betas

plot(mod.seq, type = "cloglog")
dev.off()

pred.seq <- eval.predictions(e.mx.l)[[opt.seq$tune.args]]
plot(pred.seq)

#........

eval1@predictions #details of model predictions

#Now plot the model with the lowest AICc
plot(eval1@predictions[[which(eval1@results$delta.AICc==0)]], main="Relative occurrence rate")

aic.opt <- eval1@models[[which(eval1@results$delta.AICc==0)]] #defining the optimum model

aic.opt #Should open html file that you get when you run classic Maxent for optimum model

aic.opt@results #Model stats of optimum model
write.csv(aic.opt@results, "Optimum Model Stats_Clade_1_Present_30.csv") #Writes optimum model stats for Threshold values

var.importance(aic.opt) #returns the environmental variable contributions of optimum model
write.csv(var.importance(aic.opt), "Variable importance_Clade_1_Present_30.csv")
response(aic.opt)#Plots response curves
response(aic.opt, "bio.01") #If you want to plot a specific response curve

#The following codes plot various metrics related to the models run
eval.plot(eval1@results)

eval.plot(eval1@results, 'avg.test.AUC', var='var.test.AUC', legend = F)

df <- var.importance(aic.opt) #Shows environmental variable with most importance by barplot
barplot(df$permutation.importance, names.arg=df$variable, las=2, ylab="Permutation Importance")

par(mfrow=c(2,3))
eval.plot(eval1@results) 
eval.plot(eval1@results, 'avg.test.AUC', legend = F)
eval.plot(eval1@results, 'avg.diff.AUC', legend = F)
eval.plot(eval1@results, 'avg.test.orMTP', legend = F)
eval.plot(eval1@results, 'avg.test.or10pct', legend = F)
plot(eval1@results$avg.test.AUC, eval1@results$delta.AICc, bg=eval1 @results$features, pch=21, cex= eval1@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
legend("topright", legend=unique(eval1@results$features), pt.bg=eval1@results$features, pch=21)
mtext("Circle size proportional to regularization multiplier value")

#remove plot format
par(mfrow=c(1,1))

#Note that the models are all in raw format, to change to logistic or cumlative use the following script
Logmap <- predict(aic.opt, raster_stack, outputformat=logistic) 
plot(Logmap)#Plot map to see if its correct

#Dont forget to save your map
#Currentraster<-writeRaster(eval1@predictions[[which(eval1@results$delta.AICc==0)]], filename="Minor_Present.tif", format="GTiff", overwrite=TRUE) #This saves to "Documents" Folder
#or
Currentraster<-writeRaster(Logmap, filename="Clade_1_Present_30_ZAF.tif", format="GTiff", overwrite=TRUE)
#Also dont forget to change the filenames

Currentraster<-writeRaster(pred.seq, filename="cloetei_Present_30_ZAF.tif", format="GTiff", overwrite=TRUE)
#Also dont forget to change the filenames


#......***Project to different time and space***.......

palaeofiles=c("C:/Users/TaftJ/Documents/Projects/Cordylus minor/Cordylus_2020/Paleo/WorldClim 1.4 LIG_ZAF")#Load future conditions .asc files
palaeofiles <- list.files(palaeofiles, pattern='asc$', full.names=TRUE ) #This line creates a loop to read in all files within 'path' folder
names(palaeofiles) <- c("bio.01", "bio.02", "bio.03", "bio.04", "bio.05", "bio.06", "bio.07", "bio.08", "bio.09", "bio.10", "bio.11", "bio.12", "bio.13", "bio.14", "bio.15", "bio.16", "bio.17", "bio.18", "bio.19", "TRI.")
palaeo <- palaeofiles[c("bio.01", "bio.02", "bio.03", "bio.04", "bio.05", "bio.06", "bio.07", "bio.08", "bio.09", "bio.10", "bio.11", "bio.12", "bio.13", "bio.14", "bio.15", "bio.16", "bio.17", "bio.18", "bio.19", "TRI.")] #Call variables you select that aren't correlated

raster_stackpalaeo=stack(palaeo)#Stack past/future conditions
palaeo<-predict(pred.seq, raster_stackpalaeo) #convert to logistic
plot(palaeo)
palaeoraster<-writeRaster(palaeo, filename="Clade_1_LIG_30.tif", format="GTiff", overwrite=TRUE) #Change file name to save output map

palaeofiles=c("C:\\Users\\Jody\\Desktop\\Cordylus\\Paleo\\Present_2-5\\mpi\\asc")#Load future conditions .asc files
palaeofiles <- list.files(palaeofiles, pattern='asc$', full.names=TRUE ) #This line creates a loop to read in all files within 'path' folder
names(palaeofiles) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "TRI")
palaeo <- palaeofiles[c("bio1", "bio2", "bio3", "bio6", "bio12", "bio19", "TRI")] #Call variables you select that aren't correlated

raster_stackpalaeo=stack(palaeo)#Stack future conditions
palaeo<-predict(aic.opt, raster_stackpalaeo, outputformat=logistic) #convert to logistic
plot(palaeo)
palaeoraster<-writeRaster(palaeo, filename="Minor_LGM3_MPI_worldclim.tif", format="GTiff", overwrite=TRUE) #Change file name to save output map
