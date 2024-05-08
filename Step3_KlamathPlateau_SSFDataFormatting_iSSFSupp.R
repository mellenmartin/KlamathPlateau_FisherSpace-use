#### Klamath Plateau fisher space-use ####
#### Step 3 -- extracting covariates, prepping data, running iSSFs for supplementary ####
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####

rm(list=ls());gc() #clear the memory

#load required packages
library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(mapview)
library(stars)
library(rasterVis)
library(tidyr)
library(ggplot2)
library(raster)
library(lubridate)
library(amt)
library(tidyverse)

# working directory
setwd("./KlamathPlateau/")

dat <- read.csv("FisherData_SoOregon_Moriarty_190220_cleaned.csv",header=T)

dat$datetimenew <- as.POSIXct(strptime(dat$datetime, "%m/%d/%Y %H:%M"),tz="PST8PDT")
dat$timenew <- format(as.POSIXct(dat$datetimenew), format = "%H:%M:%S")
dat$datenew <- format(as.POSIXct(dat$datetimenew), format = "%m/%d/%Y")
dat$timenew <- format(as.POSIXct(dat$datetimenew), format = "%H:%M:%S")
dat$datenew <- format(as.POSIXct(dat$datetimenew), format = "%m/%d/%Y")

dat <- dat %>% dplyr::filter(collartype == "ATSWildlink" & activity > 20|
                               collartype == "Lotek_IridiumGPS"|
                               collartype == "Lotek_IridiumSwift"|
                               collartype == "Lotek_Swift")
dat <- dat[!is.na(dat$datetimenew),]

#once more, change RawDatetimeGPS to a POSIXct object

head(dat$datetimenew)

#create new columns for the consecutive x and y coordinates
#add NA to the last because there is no N+1 row
dat$X.next<-dat$easting[c(2:length(dat$easting),NA)]
dat$Y.next<-dat$northing[c(2:length(dat$northing),NA)]

#calculate step length as Euclidean distance using Pythagorean Theorem
dat$step<-sqrt((dat$X.next-dat$easting)^2 + (dat$Y.next-dat$northing)^2)

#make a new column for the previous X and Y locations
dat$X.prev<-dat$easting[c(NA,1:(length(dat$easting)-1))]
dat$Y.prev<-dat$northing[c(NA,1:(length(dat$northing)-1))]

#make a new column for the previous step length
dat$step.prev<-dat$step[c(NA,1:(length(dat$step)-1))]
hist(dat$step[dat$step<2000])
hist(dat$step[dat$step<2])
dat <- dat %>% dplyr::filter(!step <2)

# stand polygons
ecognition <- sf::st_read("./shapefiles/Ecognition_CoreData_220731.shp", quiet=TRUE)
ecognition <- sf::st_transform(ecognition, crs = 26910)

# read rasters in
width_pxl <- raster("./rasters/width_pxl.tif")
asymmetry <- raster("./rasters/asymmetry.tif")
mnelevp9 <- raster("./rasters/mnelevp9.tif")
sdelevp95 <- raster("./rasters/sdelevp95.tif")
mndiffp95 <- raster("./rasters/mndiffp95.tif")
rumple <- raster("./rasters/rumple.tif")
cc2m <- raster("./rasters/1st_cov_abv_2m_RV_K_UR_SS_KL_Null0.tif")
cc <- raster("./rasters/1st_cov_abv_mean_RV_K_UR_SS_KL_Null0.tif")
height <- raster("./rasters/elev_p95_RV_K_UR_SS_KL_Null0.tif") 
rumple <- raster("./rasters/canopy_rumple_RV_K_UR_SS_KL.tif") 
youngdist <- raster("./rasters/distyoungforest.tif")

cc2m_mask <- mask(x = cc2m, mask = shp2)
cc_mask<- mask(x = cc, mask = shp2)
height_mask <- mask(x = height, mask = shp2)

cc2m_sd <- focal(cc2m_mask, w=matrix(1,3,3), fun=sd, na.rm = TRUE)
cc_sd <- focal(cc, w=matrix(1,3,3), fun=sd, na.rm = TRUE)
height_sd <- focal(height, w=matrix(1,3,3), fun=sd, na.rm = TRUE)

cc2m_mean <- focal(cc2m, w=matrix(1,3,3), fun=mean, na.rm = TRUE)
cc_mean <- focal(cc, w=matrix(1,3,3), fun=mean, na.rm = TRUE)
height_mean <- focal(height, w=matrix(1,3,3), fun=mean, na.rm = TRUE)

cc2m_meandiff <- cc2m_mask - cc2m_mean
cc_meandiff <- cc_mask - cc_mean
height_meandiff <- height_mask - height_mean

##### reformat data for ssf with all animals ######

dat_all <- dat %>% 
  nest(-"uniqueID")

set.seed(1234)

dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
    make_track(d, easting, northing, datetimenew, sex = sex, crs = 26910) %>% 
      track_resample(rate = minutes(15), tolerance = minutes(30)) %>% 
      amt::filter_min_n_burst(min_n = 3) %>%
      amt::steps_by_burst() %>% amt::random_steps(n = 20) %>%
      amt::time_of_day(include.crepuscule = FALSE) %>%
      amt::extract_covariates(youngdist, where = "start") %>% 
      rename(youngdist_start = "distyoungforest") %>% 
      amt::extract_covariates(youngdist, where = "end") %>% 
      rename(youngdist_end = "distyoungforest") %>%
      amt::extract_covariates(asymmetry, where = "start") %>% 
      rename(asym_start = "asymmetry") %>% 
      amt::extract_covariates(asymmetry, where = "end") %>% 
      rename(asym_end = "asymmetry") %>% 
      amt::extract_covariates(width_pxl, where = "start") %>% 
      rename(width_start = "width_pxl") %>% 
      amt::extract_covariates(width_pxl, where = "end") %>% 
      rename(width_end = "width_pxl") %>%  
      amt::extract_covariates(height_meandiff, where = "start") %>% 
      rename(hgtdiff_start = "layer") %>% 
      amt::extract_covariates(height_meandiff, where = "end") %>% 
      rename(hgtdiff_end = "layer") %>% 
      amt::extract_covariates(cc_meandiff, where = "start") %>% 
      rename(ccdiff_start = "layer") %>% 
      amt::extract_covariates(cc_meandiff, where = "end") %>% 
      rename(ccdiff_end = "layer") %>% 
      amt::extract_covariates(height_sd, where = "start") %>% 
      rename(heightsd_start = "layer") %>% 
      amt::extract_covariates(height_sd, where = "end") %>% 
      rename(heightsd_end = "layer") %>% 
      amt::extract_covariates(cc_sd, where = "start") %>% 
      rename(ccsd_start = "layer") %>% 
      amt::extract_covariates(cc_sd, where = "end") %>% 
      rename(ccsd_end = "layer") %>% 
      amt::extract_covariates(cc, where = "start") %>% 
      rename(cc_start = "X1st_cov_abv_mean_RV_K_UR_SS_KL_Null0") %>% 
      amt::extract_covariates(cc, where = "end") %>% 
      rename(cc_end = "X1st_cov_abv_mean_RV_K_UR_SS_KL_Null0") %>% 
      amt::extract_covariates(height, where = "start") %>% 
      rename(height_start = "elev_p95_RV_K_UR_SS_KL_Null0") %>% 
      amt::extract_covariates(height, where = "end") %>% 
      rename(height_end = "elev_p95_RV_K_UR_SS_KL_Null0") %>% 
      #st_distance(.,ecognition) %>% 
      mutate(log_sl_=log(sl_),
             cos_ta_ = cos(ta_))%>%
      filter(is.infinite(log_sl_)!=TRUE) %>% 
      na.omit()
  }))

m1 <- dat_all %>% mutate(ssf = lapply(trk, function(x) {
  x %>% amt::fit_issf(case_ ~ scale(cc_end) + I(scale(cc_end)^2)+ scale(ccdiff_end) + scale(youngdist_end) +
                        scale(height_end) +I(scale(height_end)^2) + scale(hgtdiff_end) +
                        sl_ + log_sl_ + cos_ta_ + 
                        strata(step_id_), model = TRUE)
}))

#### unnest data and write for later use ####
dat_ssf <- m1 %>% 
  select(uniqueID, trk) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(uniqueID)), 
    step_id = paste0(id,"-", burst_))

dat_ssf$step_id2 <- paste0(as.numeric(factor(dat_ssf$uniqueID)),"-", dat_ssf$burst_)

dat_ssf

write.csv(dat_ssf, file = "KlamathPlateau_FisherSSFLocs_UsedandRandom.csv" , row.names = FALSE)
