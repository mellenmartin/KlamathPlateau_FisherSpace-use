#### Klamath Plateau fisher space-use ####
#### Step 2 -- dBBMMs and motion variance ####
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####
### code derived from code written by Jerod Merkle 
### shared during 2019 British Ecological Society movement ecology workshop


rm(list=ls());gc() 

# set working drive
setwd("./KlamathPlateau/")


library(dplyr)
library(tidyr)
library(lubridate)
#load the data
fisher <- read.csv("FisherData_SoOregon_Moriarty_190220_cleaned.csv",header=T)
head(fisher)
fisher <- fisher %>% filter(case_ == "TRUE")

# packages
library(rgdal)
library(move)
library(adehabitatHR)
library(mapview)
library(stringr)
library(snowfall)

#specify the format of the date/time label and the time zone, and convert to POSIX
fisher$datetime <- as_datetime(fisher$datetimenew)
#fisher <- fisher %>% drop_na(start_easting, start_northing, datetime)

#convert to a spatial points data frame
coordinates(fisher) <- c("easting","northing")
proj4string(fisher) <- "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +init=epsg:26910"

#-------------------------------------------#
# create Brownian Bridge movement models ####
#-------------------------------------------#

# We need to mark whether there are locations that are greater than 8 or so hours apart (BB analyses choke on this)
fisher <- fisher[order(fisher$uniqueID, fisher$datetimenew),]  #order the database by id and date
dif <- c(as.numeric(diff(as.numeric(fisher$datetime))),0)  #change in time between points (in hrs)
dif[c(diff(as.numeric(as.factor(paste(fisher$uniqueID)))),0) != 0] <- 0 #replace the difs when it switches individual/season with a 0
fisher$dif <- dif

hist(dif)
hist(dif[dif < 10000])
table(dif > 14400)    #how often are your points > 8 hour apart?

#add a burst column for when there are breaks larger than 4 hours in the data 
MaxFixInterval <- 14400     # identify your max interval as an object because we'll use it later
fisher$connect <- ifelse(fisher$dif > MaxFixInterval, "no", "yes")
table(fisher$connect)
rm(dif)


#---------------#
# Dynamic BB ####
#---------------#

#first, create raster/grid to calculate UDs over
ext <- extent(fisher)
multiplyers <- c((ext[2]-ext[1])*0.2, (ext[4]-ext[3])*0.2)   # add about 20% around the edges of your extent (you can adjust this if necessary)
ext <- extend(ext, multiplyers)
grd <- raster(ext)
res(grd) <- 100
#i'm using a 250m resolution here. Might want to increase this if things are going slowly. Or decrease it if you want more precision
projection(grd) <- proj4string(fisher)
rm(multiplyers, ext)
#plot your grid
plot(extent(grd))   # this is the bounding box of your grid
plot(fisher, add=TRUE)

library(snowfall)

#create a folder to write your BBs to
fldr <- "./KlamathPlateau/dBBMMs"   
if(dir.exists(fldr)==FALSE){   # if it doesn't exist, it'll create it for you
  dir.create(fldr)
}
# check if there any othr files in it
length(dir(fldr))   #this should be 0. Otherwise you may want to remove files from this folder.

u <- unique(fisher$uniqueID) #loop over ids
# Create list to store the dataframes within the loop

# run it on multiple processors, as it can take considerable amount of time!
snowfall::sfInit(parallel = T, cpus = 2)   #must change the cpus
snowfall::sfExport("fisher", "grd", "u", "fldr", "bbmm.list")
snowfall::sfLibrary(move)

DynBB <- do.call(rbind, snowfall::sfClusterApplyLB(1:length(u), function(i){
  temp <- fisher[fisher$uniqueID==u[i],]
  # this is the function to calculate the dynamics BB
  mov <- move(x=coordinates(temp)[,1], y=coordinates(temp)[,2], time=temp$datetime,
              animal=u[i],proj=CRS(proj4string(temp)))   #create mov object
  mov <- move::burst(x=mov, f=temp$connect[1:(nrow(temp)-1)])   #this identifies the bad points (ie > than your MaxFixInterval)
  Dbb <- try(brownian.bridge.dyn(mov,
                                 location.error=39, #this is the location error of your collars
                                 raster=grd,
                                 margin=11,
                                 window.size=31,
                                 burstType="yes"), 
             silent=TRUE)

  #mov$var <- getMotionVariance(Dbb)
  #temp.mov <- as.data.frame(mov)
  #temp.mov$uniqueID <- u[1]
  #write out results to file too, so if there is an error you don't loose all your work!
  if(class(Dbb)=="try-error"){
    return(data.frame(ID=u[i], Failed=TRUE))   #return an NA so you know there was an error
  }else{
    if(length(Dbb@layers)>1){   #check to see if it was a multi-part DBB
      rast <- sum(Dbb)
    }else{
      rast <- Dbb[[1]]
    }
    writeRaster(rast, filename = paste(fldr, "/dynBB_", u[i],".shp", sep=""), format="HFA", overwrite=T)   # if you want to save just the imagine file of the UD
    save(Dbb, file=paste(fldr, "/dynBB_", u[i],".RData", sep=""))   #save the BB you made as a Rdata file so you can reload it later
    return(data.frame(ID=u[i], Failed=FALSE))    # have it return the motion variance
    
  }
}))
sfStop()   #this stops the parallelization loop

head(DynBB)
# were there any BBs that didn't work?
table(DynBB$Failed)  #if any of these are TRUE, then YES

# some code to look at your dynamic BBs
whichID <- u[4]    #choose an ID to look at
rast <- raster(paste(fldr, "/dynBB_", whichID,".img", sep=""))
sum(values(rast))   #this should sum to 1
rast <- getVolumeUD(as(rast, "DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)

plot(rast, xlim=extent(fisher[fisher$uniqueID == whichID,])[1:2]+c(-1000,1000),
     ylim=extent(fisher[fisher$uniqueID == whichID,])[3:4]+c(-1000,1000))
move::contour(rast, levels=c(.5,0.95,.99), add=T)
plot(fisher[fisher$uniqueID == whichID,"uniqueID"], add=T, pch=".", col="black")

# remove excess objects
rm(rast, whichID, MaxFixInterval)


#----------------------#
# HR size for dynBB ####
#----------------------#
# run it on multiple processors


sfInit(parallel = T, cpus = 2)   #must change the cpus
snowfall::sfExport("fldr","u")
snowfall::sfLibrary(move)
dynBB_HRsizes <- do.call(cbind, snowfall::sfClusterApplyLB(1:length(u), function(i){
  rast2 <- raster::raster(paste(fldr, "/dynBB_", u[i],".img", sep=""))
  # load(paste(fldr, "/dynBB_", u[i],".RData", sep=""))
  rast2 <- move::getVolumeUD(as(rast2, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
  hr95 <- reclassify(rast2, rcl=matrix(c(0,0.95,1,0.95,1,0),2,3, byrow=T))
  hr95 <- (res(hr95)[1]*res(hr95)[2]/1000000)*table(values(hr95)==0)[1]
  hr99 <- reclassify(rast2, rcl=matrix(c(0,0.99,1,0.99,1,0),2,3, byrow=T))
  hr99 <- (res(hr99)[1]*res(hr99)[2]/1000000)*table(values(hr99)==0)[1]
  hr50 <- reclassify(rast2, rcl=matrix(c(0,0.5,1,0.5,1,0),2,3, byrow=T))
  hr50 <- (res(hr50)[1]*res(hr50)[2]/1000000)*table(values(hr50)==0)[1]
  return(c(hr50,hr95,hr99))
}))
sfStop()
round(dynBB_HRsizes,1)  #in KM^2

dyn <- as.data.frame(dynBB_HRsizes)
dyn <- t(dyn)
dyn <- as.data.frame(dyn)
dyn$ID <- u

names(dyn)[names(dyn) == "u"] <- "AnimalID"

names(dyn)[names(dyn) == "FALSE."] <- "dBBMM_50"
names(dyn)[names(dyn) == "FALSE..1"] <- "dBBMM_95"
names(dyn)[names(dyn) == "FALSE..2"] <- "dBBMM_99"

write.csv(dyn, file = "KlamathPlateau_dBBMM.csv", row.names = FALSE)

#---------------------------------------------------#
# Create polygon of contour for dynamic BB model ####
#---------------------------------------------------#
plot(DynBB)
contour(dbbmm, add=T, levels=c(.5,.95))
raster2contour(dbbmm)
show(dbbmm)

bbmm.contour = data.frame(x = dbbmm$x, y = dbbmm$y, probability = dbbmm$probability)


#' compare dbmm plot to tracks
#' We can map the data
#' turn back to lat long
trk_map <-
  mk_track(
    mend_data,
    .x = long,
    .y = lat,
    .t = time,
    id = id,
    species = species,
    crs = CRS("+init=epsg:4326")
  )

trk_map <- filter(trk_map, id == "WBV1__44782")

#' plot all of the data on the one graph
qmplot(x_,
       y_,
       data = trk_map,
       maptype = "toner-lite",
       colour = id, 
       legend = "none",
       xlab = "longitude", 
       ylab = "latitude")


bb.95 <- getverticeshr(tata, percent = 95)
bb.95


#---------------------------------------------------#
# export contour polygons for dynamic BB model ####
#---------------------------------------------------#

    #choose an ID to look at
percentile99 <- 0.99   #what contour percentile are you interested in?
percentile95 <- 0.95
percentile50 <- 0.50
# set output directory for KDE loop to write files to #
dirout <- "./KlamathPlateau/dBBMMs/"


for(i in 1:length(u)){
  rast <- raster(paste(fldr, "/dynBB_", u[i],".img", sep=""))
  rastD <- getVolumeUD(as(rast, Class="DBBMM"))    #turn it into a volume, where it starts with the highest use areas (uses move package)
  poly99 <- reclassify(rastD, rcl=matrix(c(0,percentile99,1,percentile99,1,0),2,3, byrow=T))
  poly99[values(poly99)==0] <- NA
  poly99 <- rasterToPolygons(poly99, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour

  poly95 <- reclassify(rastD, rcl=matrix(c(0,percentile95,1,percentile95,1,0),2,3, byrow=T))
  poly95[values(poly95)==0] <- NA
  poly95 <- rasterToPolygons(poly95, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour
  
  poly50 <- reclassify(rastD, rcl=matrix(c(0,percentile50,1,percentile50,1,0),2,3, byrow=T))
  poly50[values(poly50)==0] <- NA
  poly50 <- rasterToPolygons(poly50, dissolve=TRUE) #this is the spatialPolygonsDataFrame representing the contour
  
  rgdal::writeOGR(poly99, dirout, paste(u[i], "_dynBB99", sep=""), "ESRI Shapefile")
  rgdal::writeOGR(poly95, dirout, paste(u[i], "_dynBB95", sep=""), "ESRI Shapefile")
  rgdal::writeOGR(poly50, dirout, paste(u[i], "_dynBB50", sep=""), "ESRI Shapefile")
  
  }
  #plot it

#---------------------------------------------------#
# get motion variance  ####
#---------------------------------------------------#
bbmm.list <- list()

for(i in 1:length(u)){
  temp <- fisher[fisher$uniqueID==u[i],]
  # this is the function to calculate the dynamics BB
  mov <- move(x=coordinates(temp)[,1], y=coordinates(temp)[,2], time=temp$datetime,
              animal=u[i],proj=CRS(proj4string(temp)))   #create mov object
  mov <- move::burst(x=mov, f=temp$connect[1:(nrow(temp)-1)])   #this identifies the bad points (ie > than your MaxFixInterval)
  Dbb <- try(brownian.bridge.dyn(mov,
                                 location.error=39, #this is the location error of your collars
                                 raster=grd,
                                 margin=11,
                                 window.size=31,
                                 burstType="yes"), 
             silent=TRUE)
  
  #### create data frame of the locs and motion variance ####
  mov$var <- getMotionVariance(Dbb)
  temp.mov <- as.data.frame(mov)
  temp.mov$uniqueID <- u[i]
  #### add to list ####
  tmp <- list(temp.mov)
  # add to list
  bbmm.list[[i]] <- temp.mov
}

#### create a table with the coords and the motion variance ###
library(data.table)
combinemot <- do.call(what = rbind, args = bbmm.list)
combinemot <- as.data.frame(combinemot)
write.csv(combinemot, file = "KlamathPlateau_dBBMMMotionVariance.csv", row.names = FALSE)
