#Wallace EcoMod v1.9.3-6
#last-updated:2022-07-25
#author: Isabel de Silva
#purpose: Species distribution modeling prototype via Wallace for North Central
##RISCC project

### Set-up----

#set working directory
setwd('/Users/isabel/Documents/work/cu_boulder/GRA/earthlab/data/sdms_niches/wallacev2/20220725/')

#load dependencies
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(rgdal)
library(tidyverse)
library(wallace)
library(foreign)
library(rgdal)
library(terra)

# Load Species List ----
spe<-read.csv('lc_spp_plants_stateinvasives.csv')

#prep for creating species subfolders in for loop
spec<-spe$SP_GNAME_ECO_UNFORMATTED

#loop through running an SDM for each species
for (i in seq_along(spec)){
  folder<-dir.create(paste0("/Users/isabel/Documents/work/cu_boulder/GRA/earthlab/data/sdms_niches/wallacev2/20220725/","/",spec[i])) #create species subfolders
  path<-paste0("/Users/isabel/Documents/work/cu_boulder/GRA/earthlab/data/sdms_niches/wallacev2/20220725/","/",spec[i]) #change wd for each species
  setwd(path)


### Write Output Folders ----

wd<-getwd()#save working directory path

#for occurrence outputs
dir.create(paste0(wd,"/","occurrences")) #create occurrences subfolder

#for climate outputs
dir.create(paste0(wd,"/","climate")) #create climate subfolder
dir.create(paste0(wd,"/","climate","/","current")) #create current climate subfolder
dir.create(paste0(wd,"/","climate","/","future")) #create future climate subfolder
dir.create(paste0(wd,"/","climate","/","future","/","2050")) #create future 2050 climate subfolder
dir.create(paste0(wd,"/","climate","/","future","/","2070")) #create future 2070 climate subfolder


#for maxent outputs
dir.create(paste0(wd,"/","maxent")) #create maxent subfolder
dir.create(paste0(wd,"/","maxent","/","plots")) #create maxent plots subfolder
dir.create(paste0(wd,"/","maxent","/","spatial")) #create maxent spatial subfolder
dir.create(paste0(wd,"/","maxent","/","results")) #create maxent results subfolder



### Obtain Occurrence Data ----

# Search the gbif database, limited to 100000 records per species.
# You decided to remove occurrences without uncertainty
# information? FALSE

# Query selected database for occurrence records
queryDb_Ld <- occs_queryDb(
  spNames = spec[i],
  occDb = "gbif",
  occNum = 100000,
  RmUncertain = FALSE)

occs_Ld <-queryDb_Ld[[paste0(word(spec[i],1),"_",word(spec[i],2))]]$cleaned

### Obtain environmental data ---

# Using WorldClim (<http://www.worldclim.org/>) bioclimatic dataset at
# resolution of 2.5 arcmin.

envs_Ld <- envs_worldclim(
  bcRes = 2.5,
  bcSel = c('bio01', 'bio02', 'bio03', 'bio04', 'bio05', 'bio06', 'bio07', 'bio08', 'bio09', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19'),
  mapCntr = c(-115.598, 45.814), # Mandatory for 30 arcsec resolution
  doBrick = TRUE)
occs_xy_Ld <- occs_Ld[c('longitude', 'latitude')]
occs_vals_Ld <- as.data.frame(raster::extract(envs_Ld, occs_xy_Ld, cellnumbers = TRUE))
# Remove duplicated same cell values
occs_Ld <- occs_Ld[!duplicated(occs_vals_Ld[, 1]), ]
occs_vals_Ld <- occs_vals_Ld[!duplicated(occs_vals_Ld[, 1]), -1]
# remove occurrence records with NA environmental values
occs_Ld <- occs_Ld[!(rowSums(is.na(occs_vals_Ld)) >= 1), ]
# also remove variable value rows with NA environmental values
occs_vals_Ld <- na.omit(occs_vals_Ld)
# add columns for env variable values for each occurrence record
occs_Ld <- cbind(occs_Ld, occs_vals_Ld)

#make a spatial points object just to have for plotting
occssp<-SpatialPoints(occs_Ld[,c(3,4)])
#plot(occssp)

### Process Occurrence Data ---

#Remove occurrences outside of user drawn polygon (rectangular bounding box of US)
occs_Ld <- poccs_selectOccs(
  occs = occs_Ld,
  polySelXY = matrix(c(-125.027847, -66.492691, -66.31691, -125.643082, -125.027847, 49.080613, 49.539023, 24.666362, 25.858606, 49.080613), ncol = 2, byrow = FALSE),
  polySelID = 1165)

#Thin the occurrences to 1 km
occs_Ld <- poccs_thinOccs(
  occs = occs_Ld,
  thinDist = 1)

#check on location of occurrences plotted with counterminuous US
#read in and transform counterminuous US
conus<-readOGR('/Users/isabel/Documents/work/cu_boulder/GRA/earthlab/data/USA/conterminous/conterminuous.shp')
#plot(conus)
conust<-spTransform(conus,proj4string(envs_Ld))
#plot(conust)
#make spatial points of cleaned occurences and overlay on counterminuous US
occspp<-SpatialPoints(occs_Ld[,c(3,4)])
#plot(occspp,add=T)



# Background points ----
# Sampling of 100000 background points and corresponding environmental data
# using a user drawn background extent with a 0 degree buffer.

# Create a background extent based on user drawn polygon
bgExt_Ld <- penvs_drawBgExtent(
  polyExtXY = matrix(c(-125.027847, -66.492691, -66.31691, -125.643082, -125.027847, 49.080613, 49.539023, 24.666362, 25.858606, 49.080613), ncol = 2, byrow = FALSE),
  polyExtID = 4903,
  drawBgBuf = 0,
  occs = occs_Ld)


### Process environmental data ---
# Mask environmental data to provided extent
bgMask_Ld <- penvs_bgMask(
  occs = occs_Ld,
  envs = envs_Ld,
  bgExt = bgExt_Ld)
# Sample background points from the provided area
bgSample_Ld <- penvs_bgSample(
  occs = occs_Ld,
  bgMask =  bgMask_Ld,
  bgPtsNum = 100000)
# Extract values of environmental layers for each background point
bgEnvsVals_Ld <- as.data.frame(raster::extract(bgMask_Ld,  bgSample_Ld))
##Add extracted values to background points table
bgEnvsVals_Ld <- cbind(scientific_name = paste0("bg_", spec[i]), bgSample_Ld,
                       occID = NA, year = NA, institution_code = NA, country = NA,
                       state_province = NA, locality = NA, elevation = NA,
                       record_type = NA, bgEnvsVals_Ld)


### Partition occurrence data ----

# Partition occurrences and background points for model training and
# validation using random k-fold, a non-spatial partition method.
groups_Ld <- part_partitionOccs(
  occs = occs_Ld ,
  bg =  bgSample_Ld,
  method = "rand",
  kfolds = 4)

#Export occurrence and climate data ----

#raw occurrences
#write.csv(occs_xy_Ld,paste0(wd,"/","occurrences","/",'occs_xy_Ld.csv'),row.names=F)
save(occs_xy_Ld,file = paste0(wd,"/","occurrences","/",'occs_xy_Ld.Rdata'),row.names=F)

#cleaned occurrences with Worldclim data
#write.csv(occs_Ld,paste0(wd,"/","occurrences","/",'occs_Ld.csv'),row.names=F)
save(occs_Ld,file = paste0(wd,"/","occurrences","/",'occs_Ld.Rdata'),row.names=F)

#raster of Wordclim present data
writeRaster(envs_Ld,paste0(wd,"/","climate","/","current","/",'envs_Ld.grd'),overwrite=T,format="raster") #.grd version preserves layer names as opposed to .tif


### Build and Evaluate Niche Model ----

# Generating a species distribution model using the maxnet algorithm as
# implemented in ENMeval V2.0 (with clamping = TRUE). For tuning using L,
# LQ, H, LQH, LQHP feature classes and regularization multipliers in the
# 0.5, 10 range increasing by 1. Not using any categorical predictor
# variables.

# Run maxent model for the selected species
model_Ld <- model_maxent(
  occs = occs_Ld,
  bg = bgEnvsVals_Ld,
  user.grp = groups_Ld,
  bgMsk = bgMask_Ld,
  rms = c(0.5, 10),
  rmsStep =  1,
  fcs = c('L', 'LQ', 'H', 'LQH', 'LQHP'),
  clampSel = TRUE,
  algMaxent = "maxnet",
  parallel = TRUE,
  numCores = 7)

# Maxent model section ----
#get a sense of all models
model_Ld@results #use model results for model selection
model_Ld@results.partitions
model_Ld@models

#note distinct models
n_distinct(model_Ld@results$tune.args)
model_Ld@results$tune.args

#on choosing the "best model" from: https://wallaceecomod.github.io/wallace/articles/tutorial-v2.html
# There is a mountain of literature about this, and there is really no single answer for all datasets.
# The model performance statistics AUC (Area Under the Curve), OR (Omission Rate), and CBI (Continuous
# Boyce Index) were calculated and averaged across our partitions, and AICc (corrected Akaike information
# criterion) was instead calculated using the model prediction of the full background extent
# (and all of the thinned occurrence points). Although AICc does not incorporate the cross-validation
# results, it does explicitly penalize model complexity—hence, models with more parameters tend to have
# a worse AICc score. It’s really up to the user to decide, and the guidance text has some references
# which should help you learn more.
#
# The evaluation metrics table can be sorted.
# First, we will prioritize models that omitted few occurrence points in the predicted area during cross-validation.
# Sort the results table in ascending order by “or.10p.avg”, or the average omission rate when applying a 10-percentile
# training presence threshold to the (withheld) validation data (see guidance text for details). As we would prefer
# a model that does not omit many withheld occurrences when it makes a range prediction, we are prioritizing low values of “or.10p.avg”.
minomitocc<-which.min(model_Ld@results$or.10p.avg)
minomitocc

#Let’s also look at average validation AUC values (where higher values are better)
maxAUC<-which.max(model_Ld@results$auc.val.avg)
maxAUC

#And AICc (where lower values are better)
minAIC<-which.min(model_Ld@results$delta.AICc)
minAIC

#proceed with minAIC to select best model
bm<-model_Ld@results[minAIC,]
b<-bm$tune.args
bm$tune.args
b<-as.character(b[1])
b #model name


bestmodel<-model_Ld@models[[minAIC]]
#note that this is essentially equivalent to doing this manually, given model name, e.g.:
#bestmodel<-model_Ld@models$fc.LQHP_rm.1.5

#plot(bestmodel)
#note that this is essentially equivalent to doing this manually, given model name, e.g.:
#plot(model_Ld@models$fc.LQHP_rm.1.5)


#extract means, mins and maxes for bioclim variables for best model
sm<-as.data.frame(bestmodel$samplemeans)
names(sm)[1]<-"samplemeans"
vma<-bestmodel$varmax
names(vma)[1]<-"varmax"
vmi<-bestmodel$varmin
names(vmi)[1]<-"varmin"

bestmodelbioclim<-cbind(sm,vma,vmi) #make a dataframe

#export best model results
#write.csv(bestmodelbioclim,paste0(wd,"/","maxent","/","results","/",'bestmodelbioclim.csv'),row.names = T)
save(bestmodelbioclim,file = paste0(wd,"/","maxent","/","results","/",'bestmodelbioclim.Rdata'),row.names=F)

#export all model results
allmods<-as.data.frame(model_Ld@results)
#write.csv(allmods,paste0(wd,"/","maxent","/","results","/",'allmods.csv'),row.names = T)
save(allmods,file = paste0(wd,"/","maxent","/","results","/",'allmods.Rdata'),row.names=F)


# Visualize ----

#Generate a map of the maxnet generated model with no threshold
# Select current model and obtain raster prediction
m_Ld <- model_Ld@models[["fc.L_rm.0.5"]]
predSel_Ld <- predictMaxnet(m_Ld, bgMask_Ld,
                            type = "cloglog",
                            clamp = FALSE)
#Get values of prediction
mapPredVals_Ld <- getRasterVals(predSel_Ld, "cloglog")
#Define colors and legend
rasCols <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
legendPal <- colorNumeric(rev(rasCols), mapPredVals_Ld, na.color = 'transparent')
rasPal <- colorNumeric(rasCols, mapPredVals_Ld, na.color = 'transparent')
#Generate map
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m  %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Training)",
                     values = mapPredVals_Ld, layerId = "train",
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  #add occurrence data
  addCircleMarkers(data = occs_Ld, lat = ~latitude, lng = ~longitude,
                   radius = 5, color = 'red', fill = TRUE, fillColor = "red",
                   fillOpacity = 0.2, weight = 2, popup = ~pop) %>%
  ##Add model prediction
  addRasterImage(predSel_Ld, colors = rasPal, opacity = 0.7,
                 group = 'vis', layerId = 'mapPred', method = "ngb") %>%
  ##add background polygons
  addPolygons(data = bgExt_Ld,fill = FALSE,
              weight = 4, color = "blue", group = 'proj')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'predictedsuitability.pdf'))
plot(predSel_Ld)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'predictedsuitability_hist.pdf'))
hist(predSel_Ld)
dev.off()

#write raster of maxent model
writeRaster(predSel_Ld,paste0(wd,"/","maxent","/","spatial","/",'predSel_Ld.grd'),overwrite=T)


### Transfer model: I : 2050, RCP 8.5 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2050 using a “CC” (CCSM4)
# GCM and an RCP of *8.5*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 60)[1]),
  rcp = 85,
  model = "CC",
  year = 50)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP8.5_2050.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP8.5_2050_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP8.5_2050.grd'),overwrite=T)


### Transfer model: II : 2050, RCP 6.0 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2050 using a “CC” (CCSM4)
# GCM and an RCP of *6.0*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 60)[1]),
  rcp = 60,
  model = "CC",
  year = 50)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP6.0_2050.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP6.0_2050_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP6.0_2050.grd'),overwrite=T)


### Transfer model: III : 2050, RCP 4.5 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2050 using a “CC” (CCSM4)
# GCM and an RCP of *4.5*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 45)[1]),
  rcp = 45,
  model = "CC",
  year = 50)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP4.5_2050.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP4.5_2050_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP4.5_2050.grd'),overwrite=T)


### Transfer model: IV : 2050, RCP 2.6 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2050 using a “CC” (CCSM4)
# GCM and an RCP of *2.6*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 26)[1]),
  rcp = 26,
  model = "CC",
  year = 50)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP2.6_2050.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP2.6_2050_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP2.6_2050.grd'),overwrite=T)


### Transfer model: V : 2070, RCP 8.5 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2070 using a “CC” (CCSM4)
# GCM and an RCP of *8.5*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 60)[1]),
  rcp = 85,
  model = "CC",
  year = 70)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP8.5_2070.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP8.5_2070_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP8.5_2070.grd'),overwrite=T)


### Transfer model: VI : 2070, RCP 6.0 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2070 using a “CC” (CCSM4)
# GCM and an RCP of *6.0*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 60)[1]),
  rcp = 60,
  model = "CC",
  year = 70)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP6.0_2070.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP6.0_2070_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP6.0_2070.grd'),overwrite=T)


### Transfer model: VII : 2070, RCP 4.5 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2070 using a “CC” (CCSM4)
# GCM and an RCP of *4.5*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 45)[1]),
  rcp = 45,
  model = "CC",
  year = 70)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP4.5_2070.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP4.5_2070_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP4.5_2070.grd'),overwrite=T)


### Transfer model: VIII : 2070, RCP 2.6 ----

# Transferring the model to the same modelling area with no threshold rule.
# Note no threshold gives a continuous suitability output as opposed to a binary one.
# New time based on “WorldClim 1.4” variables for 2070 using a “CC” (CCSM4)
# GCM and an RCP of *2.6*.

#Download variables for transferring
xferTimeEnvs_Ld <- raster::getData(
  'CMIP5',
  var = "bio",
  res=2.5,
  #res = as.integer((raster::res(bgMask_Ld) * 26)[1]),
  rcp = 26,
  model = "CC",
  year = 70)

names(xferTimeEnvs_Ld) <- paste0('bio', c(paste0('0',1:9), 10:19))
# Select variables for transferring to match variables used for modelling
xferTimeEnvs_Ld <- xferTimeEnvs_Ld[[names(bgMask_Ld)]]

# Generate a transfer of the model to the desired area and time
xfer_time_Ld <-xfer_time(
  evalOut = model_Ld,
  curModel = "fc.L_rm.0.5",
  envs = xferTimeEnvs_Ld,
  xfExt = bgExt_Ld,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = FALSE
)
# store the cropped variables of transfer
xferExt_Ld <- xfer_time_Ld$xferExt


###Make map of transfer
bb_Ld <-  bgExt_Ld@bbox
bbZoom <- polyZoom(bb_Ld[1, 1], bb_Ld[2, 1], bb_Ld[1, 2],
                   bb_Ld[2, 2], fraction = 0.05)
mapXferVals_Ld <- getRasterVals(xfer_time_Ld$xferTime,"cloglog")
rasCols_Ld <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Ld), mapXferVals_Ld, na.color = 'transparent')
rasPal_Ld <- colorNumeric(rasCols_Ld, mapXferVals_Ld, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap)
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
                     title = "Predicted Suitability<br>(Transferred)",
                     values = mapXferVals_Ld, layerId = 'xfer',
                     labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
  # map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Ld$xferTime, colors = rasPal_Ld, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
  ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Ld, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')

#write map to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP2.6_2070.pdf'))
plot(xfer_time_Ld$xferTime)
dev.off()

#write histogram of cell counts to pdf
pdf(paste0(wd,"/","maxent","/","plots","/",'transfer_CCSM4_RCP2.6_2070_hist.pdf'))
hist(xfer_time_Ld$xferTime)
dev.off()

#export raster of transfer to new time - both Wordclim data (as RasterBrick)
#and maxent suitability (as RasterLayer)
writeRaster(xfer_time_Ld$xferTime,paste0(wd,"/","maxent","/","spatial","/",'transfer_CCSM4_RCP2.6_2070.grd'),overwrite=T)

}

