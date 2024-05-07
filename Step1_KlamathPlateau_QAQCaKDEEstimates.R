#### Klamath Plateau fisher space-use ####
#### Step 1 -- data grooming & aKDE summaries ####
#### data collected by Moriarty, Matthews, et al. ####

# packages
library(foreign)
library(tidyverse)
library(dplyr)
library(lubridate)
library(raster)
library(mapview)
library(sf)
library(amt)
library(ggplot2) 
library(ctmm)

# directory
setwd("./KlamathPlateau/")
dat <- read.csv("FisherData_SoOregon_Moriarty_190220_memedit.csv",header=T)
names(dat) = tolower(names(dat)) 
# create unique ID for individual and deployment
dat$uniqueID <- paste0(dat$fid_,"-",dat$deployment)

# drop outliers, unsuccessful, and low precision fixes
dat <- dat %>% drop_na(easting, northing, date, time, datetime)
dat <- dat %>% filter(easting > 1,
                      northing > 1,
                      !easting<450000,
                      !easting>578000,
                      collartype == "ATSWildlink" & hdop < 5|
                        collartype == "Lotek_IridiumGPS" & hdop < 5|
                        collartype == "Lotek_IridiumSwift" & hdop < 5|
                        collartype == "Lotek_Swift" & hdop < 5)

####
# take a look at the points in mapview to find and drop outliers
fish <- dat
coordinates(fish) <- c("easting","northing")
proj4string(fish) <- "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +init=epsg:26910"
mapview(fish, zcol = "uniqueID", burst = TRUE)

# drop duplicates and outliers for each individual 
dat <- dat %>%  filter(!between(id,7673,7681),
                       !id == 26089,
                       !id == 33728,
                       !id == 59248,
                       !id == 57121,
                       !id == 83436,
                       !id == 83437,
                       !id == 90919,
                       !id == 90921,
                       !id == 91515,
                       !id == 91516,
                       !id == 102093,
                       !id == 102095,
                       !id == 102097,
                       !id == 102098,
                       !id == 102099,
                       !between(id,102101,102103))

####
fish <- dat
coordinates(fish) <- c("easting","northing")
proj4string(fish) <- "+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +init=epsg:26910"
mapview(fish, zcol = "uniqueID", burst = TRUE)

#specify the format of the date/time label and the time zone, and convert to POSIX
dat$datetimenew <- as.POSIXct(paste(dat$date, dat$time), format = "%m/%d/%Y %H:%M:%S")
dat$datetimenew <- as.POSIXct(strptime(dat$datetime, "%m/%d/%Y %H:%M"),tz="PST8PDT")
#dat$datetimenew <- as.POSIXct(dat$datetime, tz = "America/Los_Angeles")
#dat$dattimenew <- dat$datetimenew %m+% years(2000)
dat$timenew <- format(as.POSIXct(dat$datetimenew), format = "%H:%M:%S")
dat$datenew <- format(as.POSIXct(dat$datetimenew), format = "%m/%d/%Y")
dat$julian <- yday(as.POSIXlt(dat$datetimenew))
dat$season <- ifelse(dat$julian > 32 & dat$julian < 121, "breeding", "non-breeding")
dat <- dat %>% drop_na(easting, northing, datetimenew)
dat <- dat %>% filter(easting > 1,
                      northing > 1)

dat$sex <- substr(dat$uniqueID,1,1)
dat$sex <- ifelse(dat$sex == "F", "Female", "Male")
dat <- dat[order(dat$fid_, dat$datetimenew),]  #order the database by id and date

#### autocorrelated kernel density estimates for supplementary/BLM/USFWS summary
library(foreign)
library(tidyverse)
library(dplyr)
library(lubridate)
library(raster)
library(mapview)
library(amt)
library(sf)
library(ggplot2) 
library(ctmm)

# nest cleaned data by "uniqueID" -- animal and season -- to nest data for each animal and each collar deployment
fisher_nest <- dat %>% nest(data=-"uniqueID")

fisher_trk <- fisher_nest %>%
  mutate(trk = map(data, function(d) {
    make_track(d, easting, northing, datetimenew, sex = sex, crs = 26910) 
  }))

#### aKDE #####
#fisher_trk <- fisher_nest %>%
#  mutate(trk = map(data, function(d) {
#    amt::make_track(d, easting, northing, datetimenew, sex = sex, 
#                    crs = 4326) %>% 
#      amt::transform_coords(26910)
#  }))

system.time(
  hr_akde_ouf <- fisher_trk %>%
    mutate(
      hr_akde = map(trk, ~ hr_akde(., fit_ctmm(., "ouf")))
    )
)

iso_hr_akde_ouf <- hr_akde_ouf %>%
  mutate(
    iso_hr_akde = map(hr_akde, ~ hr_isopleths(.))
  )

iso_hr_akde_ouf_unn <- iso_hr_akde_ouf %>% dplyr::select(c(uniqueID,iso_hr_akde)) %>%
  unnest(cols = iso_hr_akde) %>%
  filter(what=="estimate") %>%
  mutate(area_km2 = as.numeric(area)/1000000)

# export the datasets
# export as shapefiles
st_write(iso_hr_akde_ouf_unn, "~/aKDEs/KlamathPlateau_aKDEEstimates.shp", append = TRUE)
#export as csv file
fisher_hr_table <- iso_hr_akde_ouf_unn %>% dplyr::select(-geometry)
write.csv(fisher_hr_table, file = "KlamathPlateau_aKDEEstimates.csv", row.names = FALSE)

##### write the cleaned dat for dBBMMs and further processing for the iSSF, etc
write.csv(dat, file = "FisherData_SoOregon_Moriarty_190220_cleaned230525.csv", row.names = FALSE)
