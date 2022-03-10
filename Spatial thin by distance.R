install.packages("spThin") #You'll need this to thin occurrence points by distance     

setwd("C:\\Users\\Jody\\Desktop\\Cordylus\\Model Outputs") #Place to save output files

#..............Reduce spatial auto-correlation in occurrence points...............................

library(spThin)

#The first step is reading in your occurrence file. Make sure you have the columns SPECIES, LAT, LONG 

loc.data <- read.csv( "C:\\Users\\Jody\\Desktop\\Cordylus\\Maxent\\Minor.csv") #Make sure the correct RAW occurrence file is selected

#A data.frame of occurence locations. It can include several columnns, 
#but must include at minimum a column of latitude values, a column of longitude values, and a column of species names.

#I change the values for thin.par to the distance I want to thin by. 
#I also change out.dir for the directory I want to save the file. 
#And don't forget to change the file name using out.base

thin (loc.data , 
      lat.col = "LAT",          #Name of column of latitude values. Caps sensitive.
      long.col = "LONG",        #Name of column of longitude values. Caps sensitive.
      spec.col = "SPECIES",     #Name of column of species name. Caps sensitive.
      thin.par = 10,            #Thinning parameter - the distance (in kilometers) that you want records to be separated by
      reps = 1,                 #The number of times to repete the thinning process. Given the random process of removing nearest-neighbors there should be 'rep' number of different sets of coordinates.
      locs.thinned.list.return = FALSE, #TRUE/FALSE - If true, the 'list' of the data.frame of thinned locs resulting from each replication is returned 
      write.files = TRUE,        #TRUE/FALSE - If true, new *.csv files will be written with the thinned locs data
      max.files = 1,             #The maximum number of *csv files to be written based on the thinned data
      out.dir = "C:\\Users\\Jody\\Desktop\\Cordylus", #Directory to write new *csv files to
      out.base = "Thinned_5km", #A file basename to give to the thinned datasets created - Change based on km distance
      write.log.file = FALSE,     #TRUE/FALSE create/append log file of thinning run
      log.file = "spatial_thin_log.txt", #Text log file
      verbose = TRUE             #TRUE/FALSE - If true, running details of the function are print at the console
)

#Your output file will be in .csv format and can be used straight in Maxent