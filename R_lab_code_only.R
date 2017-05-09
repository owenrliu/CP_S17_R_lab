## ----libraries,message=F-------------------------------------------------
# if need be, you can install the packages you don't have with the command install.packages(), 
# with the package names in quotes:
# install.packages("raster","rgdal","rasterVis","maps","rgeos","dplyr","RColorBrewer")

# Load the libraries into this R session
library(raster)       #Main raster library with nearly all functions used in this analysis
library(rgdal)        #Spatial library - most functions used from rgdal are for vectors (shapefiles)
library(rasterVis)    #Useful for raster visualizations
library(maps)         #Has a database of maps. I use this to add a map to my raster to visualize land boundaries
library(rgeos)        #Need this library for topology operations on geometries
library(dplyr)        #NOT spatial - this is a data wrangling library
library(RColorBrewer) #Also not spatial - used to set the spectral color scheme 


## ----settings------------------------------------------------------------
# rainbow color scheme
cols = rev(colorRampPalette(brewer.pal(11, 'Spectral'))(255)) 

#setting smaller margins for plotting
par(mar=c(2,2,1,1))

## ----threats,results='hide'----------------------------------------------
threats_dir <- 'E:/TA files/CP2017_Owen/R lab development/Threats_data' # Directory where all my files are. THIS WILL BE DIFFERENT FOR YOU
threat_files <- list.files(threats_dir,full.names = T) # List the files in this folder
threat_files # print the file names into the console

## ----all threats import--------------------------------------------------
all_threats <- raster(threat_files[2])

## ----all threats plot----------------------------------------------------
plot(all_threats,col=cols)

## ----basemap-------------------------------------------------------------
# add a landmap to your shapefile. the add=T argument tells R to add it to the existing plot.
# make sure you understand what the other arguments do
plot(all_threats,ext=extent(-130,-110,24,50),col=cols)
map('world',fill=T,add=T,col='gray')

## ----extent--------------------------------------------------------------
plot(all_threats,col=cols,ext=extent(-121,-117,32,35),main="Cumulative Threats") # A good extent for the Santa Barbara Channel

## ----zoom----------------------------------------------------------------
plot(all_threats,col=cols)
# zoom(all_threats,col=cols) #Interactive code. Uncomment to run.

## ----data atts-----------------------------------------------------------
all_threats

## ----hist,warning=F------------------------------------------------------
hist(all_threats,main="Cumulative Threats Frequency")

## ----cellStats-----------------------------------------------------------
cellStats(all_threats,mean)

## ----Species,results='hide'----------------------------------------------
all_spp <- raster("E:/TA files/CP2017_Owen/R lab development/Species_data/ca_curr_sp_rich.tif")
all_spp
plot(all_spp,col=cols)

## ----crop and resample---------------------------------------------------
#?crop see what the crop function does

threats_crop <- crop(all_threats,all_spp) #Crop the threats layer to the same extent at species

## ----resample,message=F--------------------------------------------------

#?resample see what the resample function does
# NOTE: the progress='text' argument is a great tool: it prints out the progress
# of a longer-running function into the console, so you can see how the operation is going

# the method='ngb' argument specifies that we want to use a nearest neighbor algorithm to resample, instead of interpolation
spp_res <- resample(all_spp,threats_crop,method='ngb',progress='text')

## ----stack---------------------------------------------------------------
spp_threat_stack <- stack(threats_crop,spp_res)
plot(spp_threat_stack,col=cols)

## ----spp hist------------------------------------------------------------
hist(spp_res,main="Species Raster Values")

## ----reclass zeroes------------------------------------------------------
# notice that in the following, we are OVERWRITING the original spp_res object.
# This is okay in this instance since we won't be using the old version, but
# often it is better to assign any output of a function to a new variable or object
spp_res <- reclassify(spp_res,rcl=c(-Inf,0,NA))
hist(spp_res,main="Species Raster Values, Zeroes Removed") # did the function do what we were hoping?

## ----top 20 percent------------------------------------------------------
#?quantile what does the quantile function do?
spp_cutoff <- quantile(spp_res,0.8) # Find the value of the 80th percentile
spp_maxVal <- cellStats(spp_res,max) #find the maximum

# Our reclassification matrix. Make sure you know what this is saying
rcl_mat <- c(-Inf,spp_cutoff,0,
            spp_cutoff,spp_maxVal,1)

# Reclassify the species layer
spp_binary <- reclassify(spp_res,rcl=rcl_mat)

## ----reclass plot--------------------------------------------------------
# Because we have binary data now, I want to change the color scheme again
binary_cols <- c("white","firebrick")
plot(spp_binary,col=binary_cols,legend=F,main="Top 20% of Species Richness")
map('world',fill=T,add=T,col='gray')

## ----threat reclass,echo=F-----------------------------------------------
#?quantile what does the quantile function do?
threat_cutoff <- quantile(threats_crop,0.8) # Find the value of the 80th percentile
threat_maxVal <- cellStats(threats_crop,max) #find the maximum

# Our reclassification matrix. Make sure you know what this is saying
rcl_mat <- c(-Inf,threat_cutoff,0,
            threat_cutoff,threat_maxVal,1)

# Reclassify the species layer
threat_binary <- reclassify(threats_crop,rcl=rcl_mat)
plot(threat_binary,col=binary_cols,legend=F,main="Top 20% of Cumulative Threats")
map('world',fill=T,add=T,col='gray')

## ----hotspots------------------------------------------------------------
# the hotspots
hotspots <- overlay(spp_binary,threat_binary,fun=function(x,y){x+y})

# color breakpoints. We need three colors now! (cell values of 0,1,or 2)
brks_hotspots <- seq(0,3,length.out=4) 
hotspot_cols <- c("white","lightblue","firebrick") #

# plot the hotspots!
plot(hotspots,col=hotspot_cols,legend=F,main="Hotspots");map('world',fill=T,add=T,col='gray80')
plot(hotspots,col=hotspot_cols,ext=extent(-121,-117,32,35),main="Hotspots, SB Channel",legend=F)
map('world',fill=T,add=T,col='gray80')

## ----all code, eval=F----------------------------------------------------
### import data ###
all_spp <- raster("E:/TA files/CP2017_Owen/R lab development/Species_data/ca_curr_sp_rich.tif")
all_threats <- raster("E:/TA files/CP2017_Owen/R lab development/Threats_data/full_modelnv.tif")

#### Crop, resample, and reclassify ###
all_spp <- reclassify(all_spp,rcl=c(-Inf,0,NA)) # reclass 0 to NA
threats_crop <- crop(all_threats,all_spp) # crop threats to species
spp_res <- resample(all_spp,threats_crop,method='ngb') # resample species to threat's resolution

#### Function to output a binary raster based on a user-given quantile (default is top 20%) ###
reclassify_topx <- function(rast,quant=0.8) {
  topx <- quantile(rast,quant) #find the 80% quantile of the raster values
  maxVal <- cellStats(rast,max) #find the maximum
  rcl <- c(-Inf,topx,0,
            topx,maxVal,1) # reclassify matrix (see help file for ?reclassify)
  out <- reclassify(rast,rcl=rcl)
  return(out) # returns the new binary raster
}

### Find top 20%, using the code from above. We could easily choose a different quantile here. ###
all_threats_top20 <- reclassify_topx(threats_crop,quant=0.8)
all_spp_top20 <- reclassify_topx(spp_res,quant=0.8)

### overlay and plot the hotspots ###
hotspots <- overlay(all_threats_top20,all_spp_top20,fun=function(x,y){x+y})
brks_hotspots <- seq(0,3,length.out=4) # color breakpoints
hotspot_cols <- c("white","lightblue","firebrick") # colors
plot(hotspots,col=hotspot_cols);map('world',fill=T,add=T,col='gray80')


