#### Klamath Plateau fisher space-use ####
#### Step 3 -- extract stand-level covariates for used and random polygons ####
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####

rm(list=ls());gc() #clear the memory

#load required packages
library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(stars)
library(rasterVis)
library(nngeo)
library(landscapemetrics)
library(landscapetools)
 
#
##### load in and bind BBMM shapefiles #####
library(sf)
library(dplyr)
wd <- "./KlamathPlateau/dBBMMs/dBBMMContours"
setwd(wd)

# load in HR contours
filenames <- list.files(wd, pattern = "*.shp$")

# Create list to store the contours within the loop
bbmm.list <- list()

## Import all bbmm shapefiles, add individual identifiers, calculate area, and calculate stream density
for(i in 1:length(filenames)){
  wd <- "./KlamathPlateau/dBBMMs/dBBMMContours"
  setwd(wd)
  
  # set CRS for the shapefiles
  crs <-  "+proj=utm +zone=10 +ellps=WGS84 +datum=NAD83 units=m"
  
  library(stringr)
  
  temp.sf <- st_read(filenames[i])
  polyname <- sub("*.shp", "", filenames[i]) 
  if(nrow(temp.sf) > 1){
    temp.sf <- temp.sf %>% slice(2)
  }
  temp.sf <- st_transform(temp.sf, crs = crs)
  temp.sf <- temp.sf %>%  rename("ID" = 1)
  temp.sf[c(1)] <- c(i)
  temp.sf$ID <- paste0(polyname)
  
  # add identifiers 
  temp.sf$HR_ID <-  sub('.*-',"", temp.sf$ID)
  temp.sf$contour <- str_sub(temp.sf$ID, - 2)
  temp.sf$Year <-  sub(".*-(.*)_.+", "\\1", temp.sf$ID)
  
  # calculate area
  temp.sf$areakm <- as.numeric(st_area(temp.sf)/1000000)
  
  library(maptools)
  library(raster)
  
  #### extract ecognition data and summarize ####
  ecog <- st_read("./KlamathPlateau/rasters/Ecognition_CoreData_220731.shp")
  ecog <- st_buffer(ecog, dist=0)
  temp.sf_ecog <- st_intersection(ecog, temp.sf) 
  temp.sf$patchcount <- length(temp.sf_ecog$OBJECTID)
  temp.sf$medianpatchsize <- median(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$minpatchsize <- min(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$maxpatchsize <- max(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$medianpatchlength <- median(temp.sf_ecog$Length_Pxl)*30
  temp.sf$minpatchlength <- min(temp.sf_ecog$Length_Pxl)*30
  temp.sf$maxpatchlength <- max(temp.sf_ecog$Length_Pxl)*30
  temp.sf$medianpatchwidth <- median(temp.sf_ecog$Width_Pxl)*30
  temp.sf$minpatchwidth <- min(temp.sf_ecog$Width_Pxl)*30
  temp.sf$maxpatchwidth <- max(temp.sf_ecog$Width_Pxl)*30
  temp.sf$medianasym <- median(temp.sf_ecog$Asymmetry)
  temp.sf$minasym <- min(temp.sf_ecog$Asymmetry)
  temp.sf$maxasym <- max(temp.sf_ecog$Asymmetry)
  temp.sf$mediancompact <- median(temp.sf_ecog$Compactnes)
  temp.sf$mincompact <- min(temp.sf_ecog$Compactnes)
  temp.sf$maxcompact <- max(temp.sf_ecog$Compactnes)
  temp.sf$medneigh <- median(temp.sf_ecog$NumNeighbo)
  temp.sf$minneigh <- min(temp.sf_ecog$NumNeighbo)
  temp.sf$maxneigh <- max(temp.sf_ecog$NumNeighbo)
  temp.sf$medrump <- median(temp.sf_ecog$Mn_CanRump)
  temp.sf$minrump <- min(temp.sf_ecog$Mn_CanRump)
  temp.sf$maxrump <- max(temp.sf_ecog$Mn_CanRump)
  temp.sf$mnheight <- median(temp.sf_ecog$Mn_elev_p9)
  temp.sf$minheight <- min(temp.sf_ecog$Mn_elev_p9)
  temp.sf$maxheight <- max(temp.sf_ecog$Mn_elev_p9)
  temp.sf$medcc <- median(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$mincc <- min(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$maxcc <- max(temp.sf_ecog$Mn_1stCo_1)
  ##
  # write the combined HR polygons w/ covariate info
  setwd("./KlamathPlateau/HRSummaries")
  st_write(temp.sf, paste0(polyname,"KlamathPlateau_BBMMCombined_HRCovs.shp"))
  
  #### add to list ####
  temp <- list(temp.sf)
  # add to list
  bbmm.list[[i]] <- temp.sf
}

library(sf)
library(data.table)
combineShp <- do.call(what = sf:::rbind.sf, args = bbmm.list)
combineShp <- st_transform(combineShp, crs = crs)

# write the combined HR polygons w/ stream info for future use
st_write(combineShp, "./KlamathPlateau/HRSummaries/KlamathPlateau_BBMMCombined_HRCovs.shp")

##### combine polygons for overall estimates ####
library(sf)
library(dplyr)

combineShp <- st_read("./KlamathPlateau/HRSummaries/KlamathPlateau_BBMMCombined.shp")
combineShp <- st_transform(combineShp, crs = crs)

# combine polygons buffer polygons by 3km, ~2x sigma for males (Green et al.)
combineShp_dissolve <- st_combine(combineShp)
combineShp_dissolve <- st_union(combineShp_dissolve, by_feature = FALSE)
combineShp_buffer <- st_buffer(combineShp_dissolve, dist = 12000)
st_write(combineShp_buffer, "./KlamathPlateau/HRSummaries/KlamathPlateau_BBMMDissolvedBuffered12km.shp")

##### summarize the patch conditions in the study area #####
combineShp_df <- combineShp

# calculate area
combineShp_df$areakm <- as.numeric(st_area(combineShp)/1000000)

##### use voronoi tesselation to create random polygons throughout the study area #####
library(tidyverse)
library(sf)

# perform first for males
malefisher_randbbmmpts <- st_sample(combineShp_buffer, 1000) %>% # random points, as a list ...
  st_sf() %>%  # ... to data frame ...
  st_transform(crs = crs)  # ... and a metric CRS

i <- 1 # iterator start

buffer_size <- 6000 # minimal distance to be enforced (in meters)

repeat( {
  #  create buffer around i-th point
  buffer <- st_buffer(malefisher_randbbmmpts[i,], buffer_size ) 
  
  offending <- malefisher_randbbmmpts %>%  # start with the intersection of master points... 
    st_intersects(buffer, sparse = F) # ... and the buffer, as a vector
  
  # i-th point is not really offending - it is the origin (not to be excluded)
  offending[i] <- FALSE
  
  # if there are any offending points left - re-assign the master points, 
  # with the offending ones excluded / this is the main pruning part :)
  malefisher_randbbmmpts <- malefisher_randbbmmpts[!offending,] 
  
  if ( i >= nrow(malefisher_randbbmmpts)) {
    # the end was reached; no more points to process
    break 
  } else {
    # rinse & repeat
    i <- i + 1 
  }
  
} )

malefisher_randbbmmpolys <- malefisher_randbbmmpts %>%  # consider the master points
  st_geometry() %>% # ... as geometry only (= throw away the data items)
  st_union() %>% # unite them ...
  st_voronoi() %>% # ... and perform the voronoi tessellation
  st_collection_extract(type = "POLYGON") %>% # select the polygons
  st_sf(crs = crs) %>% # set metric crs
  st_intersection(combineShp_buffer) %>% # 
  st_join(malefisher_randbbmmpts) %>% # & re-connect the data items
  st_set_agr("aggregate") # clean up

malefisher_randbbmmpolys_covs <- malefisher_randbbmmpolys 
# Create list to store the contours within the loop
malerand.list <- list()

for (i in 1:nrow(malefisher_randbbmmpolys_covs)){
  
  temp.sf <- st_sf(malefisher_randbbmmpolys_covs[i,])
  temp.sf$ID <- paste0(i)
  
  # add identifiers 
  temp.sf$HR_ID <-  paste0("Male_",i)
  temp.sf$contour <- paste0("NA")
  temp.sf$Year <-  paste0("NA")
  # calculate area
  temp.sf$areakm <- as.numeric(st_area(temp.sf)/1000000)
  
  #### extract ecognition data and summarize ####
  ecog <- st_read("./KlamathPlateau/rasters/Ecognition_CoreData_220731.shp")
  ecog <- st_buffer(ecog, dist=0)
  temp.sf_ecog <- st_intersection(ecog, temp.sf) 
  temp.sf$patchcount <- length(temp.sf_ecog$OBJECTID)
  temp.sf$medianpatchsize <- median(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$minpatchsize <- min(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$maxpatchsize <- max(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$medianpatchlength <- median(temp.sf_ecog$Length_Pxl)*30
  temp.sf$minpatchlength <- min(temp.sf_ecog$Length_Pxl)*30
  temp.sf$maxpatchlength <- max(temp.sf_ecog$Length_Pxl)*30
  temp.sf$medianpatchwidth <- median(temp.sf_ecog$Width_Pxl)*30
  temp.sf$minpatchwidth <- min(temp.sf_ecog$Width_Pxl)*30
  temp.sf$maxpatchwidth <- max(temp.sf_ecog$Width_Pxl)*30
  temp.sf$medianasym <- median(temp.sf_ecog$Asymmetry)
  temp.sf$minasym <- min(temp.sf_ecog$Asymmetry)
  temp.sf$maxasym <- max(temp.sf_ecog$Asymmetry)
  temp.sf$mediancompact <- median(temp.sf_ecog$Compactnes)
  temp.sf$mincompact <- min(temp.sf_ecog$Compactnes)
  temp.sf$maxcompact <- max(temp.sf_ecog$Compactnes)
  temp.sf$medneigh <- median(temp.sf_ecog$NumNeighbo)
  temp.sf$minneigh <- min(temp.sf_ecog$NumNeighbo)
  temp.sf$maxneigh <- max(temp.sf_ecog$NumNeighbo)
  temp.sf$medrump <- median(temp.sf_ecog$Mn_CanRump)
  temp.sf$minrump <- min(temp.sf_ecog$Mn_CanRump)
  temp.sf$maxrump <- max(temp.sf_ecog$Mn_CanRump)
  temp.sf$mnheight <- median(temp.sf_ecog$Mn_elev_p9)
  temp.sf$minheight <- min(temp.sf_ecog$Mn_elev_p9)
  temp.sf$maxheight <- max(temp.sf_ecog$Mn_elev_p9)
  temp.sf$medcc <- median(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$mincc <- min(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$maxcc <- max(temp.sf_ecog$Mn_1stCo_1)

  #### add to list ####
  temp <- list(temp.sf)
  # add to list
  malerand.list[[i]] <- temp.sf
  
  }

library(sf)
library(data.table)
malefisher_randbbmmpolys_covs <- do.call(what = sf:::rbind.sf, args = malerand.list)
malefisher_randbbmmpolys_covs <- st_transform(malefisher_randbbmmpolys_covs, crs = crs)

#### then repeat for females ####
femalefisher_randbbmmpts <- st_sample(combineShp_buffer, 1000) %>% # random points, as a list ...
  st_sf() %>%  # ... to data frame ...
  st_transform(crs = crs)  # ... and a metric CRS

i <- 1 # iterator start

buffer_size <- 4000 # minimal distance to be enforced (in meters)

repeat( {
  #  create buffer around i-th point
  buffer <- st_buffer(femalefisher_randbbmmpts[i,], buffer_size ) 
  
  offending <- femalefisher_randbbmmpts %>%  # start with the intersection of master points... 
    st_intersects(buffer, sparse = F) # ... and the buffer, as a vector
  
  # i-th point is not really offending - it is the origin (not to be excluded)
  offending[i] <- FALSE
  
  # if there are any offending points left - re-assign the master points, 
  # with the offending ones excluded / this is the main pruning part :)
  femalefisher_randbbmmpts <- femalefisher_randbbmmpts[!offending,] 
  
  if ( i >= nrow(femalefisher_randbbmmpts)) {
    # the end was reached; no more points to process
    break 
  } else {
    # rinse & repeat
    i <- i + 1 
  }
  
} )

femalefisher_randbbmmpolys <- femalefisher_randbbmmpts %>%  # consider the master points
  st_geometry() %>% # ... as geometry only (= throw away the data items)
  st_union() %>% # unite them ...
  st_voronoi() %>% # ... and perform the voronoi tessellation
  st_collection_extract(type = "POLYGON") %>% # select the polygons
  st_sf(crs = crs) %>% # set metric crs
  st_intersection(combineShp_buffer) %>% # limit to Prague city boundaries
  st_join(femalefisher_randbbmmpts) %>% # & re-connect the data items
  st_set_agr("aggregate") # clean up

plot(femalefisher_randbbmmpolys)

femalefisher_randbbmmpolys_covs <- femalefisher_randbbmmpolys 
# Create list to store the contours within the loop
femalerand.list <- list()

for (i in 1:nrow(femalefisher_randbbmmpolys_covs)){
  
  temp.sf <- st_sf(femalefisher_randbbmmpolys_covs[i,])
  temp.sf$ID <- paste0(i)
  
  # add identifiers 
  temp.sf$HR_ID <-  paste0("Female_",i)
  temp.sf$contour <- paste0("NA")
  temp.sf$Year <-  paste0("NA")
  # calculate area
  temp.sf$areakm <- as.numeric(st_area(temp.sf)/1000000)
  
  #### extract ecognition data and summarize ####
  ecog <- st_read("./KlamathPlateau/rasters/Ecognition_CoreData_220731.shp")
  ecog <- st_buffer(ecog, dist=0)
  temp.sf_ecog <- st_intersection(ecog, temp.sf) 
  temp.sf$patchcount <- length(temp.sf_ecog$OBJECTID)
  temp.sf$medianpatchsize <- median(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$minpatchsize <- min(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$maxpatchsize <- max(temp.sf_ecog$Shape_Area)/1000000
  temp.sf$medianpatchlength <- median(temp.sf_ecog$Length_Pxl)*30
  temp.sf$minpatchlength <- min(temp.sf_ecog$Length_Pxl)*30
  temp.sf$maxpatchlength <- max(temp.sf_ecog$Length_Pxl)*30
  temp.sf$medianpatchwidth <- median(temp.sf_ecog$Width_Pxl)*30
  temp.sf$minpatchwidth <- min(temp.sf_ecog$Width_Pxl)*30
  temp.sf$maxpatchwidth <- max(temp.sf_ecog$Width_Pxl)*30
  temp.sf$medianasym <- median(temp.sf_ecog$Asymmetry)
  temp.sf$minasym <- min(temp.sf_ecog$Asymmetry)
  temp.sf$maxasym <- max(temp.sf_ecog$Asymmetry)
  temp.sf$mediancompact <- median(temp.sf_ecog$Compactnes)
  temp.sf$mincompact <- min(temp.sf_ecog$Compactnes)
  temp.sf$maxcompact <- max(temp.sf_ecog$Compactnes)
  temp.sf$medneigh <- median(temp.sf_ecog$NumNeighbo)
  temp.sf$minneigh <- min(temp.sf_ecog$NumNeighbo)
  temp.sf$maxneigh <- max(temp.sf_ecog$NumNeighbo)
  temp.sf$medrump <- median(temp.sf_ecog$Mn_CanRump)
  temp.sf$minrump <- min(temp.sf_ecog$Mn_CanRump)
  temp.sf$maxrump <- max(temp.sf_ecog$Mn_CanRump)
  temp.sf$mnheight <- median(temp.sf_ecog$Mn_elev_p9)
  temp.sf$minheight <- min(temp.sf_ecog$Mn_elev_p9)
  temp.sf$maxheight <- max(temp.sf_ecog$Mn_elev_p9)
  temp.sf$medcc <- median(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$mincc <- min(temp.sf_ecog$Mn_1stCo_1)
  temp.sf$maxcc <- max(temp.sf_ecog$Mn_1stCo_1)
  
  #### add to list ####
  temp <- list(temp.sf)
  # add to list
  femalerand.list[[i]] <- temp.sf
  
}

library(sf)
library(data.table)
femalefisher_randbbmmpolys_covs <- do.call(what = sf:::rbind.sf, args = femalerand.list)
femalefisher_randbbmmpolys_covs <- st_transform(femalefisher_randbbmmpolys_covs, crs = crs)

##### combine male and female rand polys #####
st_geometry(femalefisher_randbbmmpolys_covs) <- NULL
st_geometry(malefisher_randbbmmpolys_covs) <- NULL
femalefisher_randbbmmpolys_covs$Sex <- "Female"
malefisher_randbbmmpolys_covs$Sex <- "Male"

allrandpolys <- rbind(femalefisher_randbbmmpolys_covs,
                      malefisher_randbbmmpolys_covs)

allrandpolys$type <- "Random"

st_geometry(combineShp) <- NULL
combineShp$Sex <- c(rep("Female",42),
                    rep("Male",33))
combineShp95 <- combineShp %>% dplyr::filter(contour == 95)
combineShp95$type <- "Used"

randandusedpolys <- rbind(combineShp95, allrandpolys)

##### home range summaries #####
st_geometry(combineShp) <- NULL
combineShp$case <- rep(1, length(combineShp$ID))

combineShp$fisher <- substr(combineShp$ID, 1, 4)
combineShp$HR_ID <- substr(combineShp$ID,1,nchar(combineShp$ID)-8)
combineShp$contour <- substr(combineShp$ID, nchar(combineShp$ID)-1, nchar(combineShp$ID))
combineShp$Sex <- substr(combineShp$fisher,1,1)
combineShp$Sex <- ifelse(combineShp$Sex == "F", "Female", "Male")
class(combineShp)

write.csv(combineShp, "./KlamathPlateau/HRSummaries/AnnualHomeRange_IndividualHRSizeLandcover.csv", row.names = FALSE)
indsums <- combineShp

options(digits=2)
annualsummaries <- randandusedpolys %>%
  group_by(Sex, type) %>%
  summarize(area_mean = mean(areakm),
            area_median = median(areakm),
            area_sd = sd(areakm),
            area_min = min(areakm),
            area_max = max(areakm),
            patctmed = median(patchcount),
            patctsd = sd(patchcount),
            patctmin = min(patchcount),
            patctmax = max(patchcount),
            patszmed = median(medianpatchsize, na.rm = TRUE),
            patszsd = sd(medianpatchsize, na.rm = TRUE),
            patszmin = min(medianpatchsize, na.rm = TRUE),
            patszmax = max(medianpatchsize, na.rm = TRUE),
            patlenmed = median(medianpatchlength, na.rm = TRUE),
            patlensd = sd(medianpatchlength, na.rm = TRUE),
            patlenmin = min(medianpatchlength, na.rm = TRUE),
            patlenmax = max(medianpatchlength, na.rm = TRUE),
            patwidmed = median(medianpatchwidth, na.rm = TRUE),
            patwidsd = sd(medianpatchwidth, na.rm = TRUE),
            patwidmin = min(medianpatchwidth, na.rm = TRUE),
            patwidmax = max(medianpatchwidth, na.rm = TRUE),
            asym_median = median(medianasym, na.rm = TRUE),
            asym_sd = sd(medianasym, na.rm = TRUE),
            asym_min = min(medianasym, na.rm = TRUE),
            asym_max = max(medianasym, na.rm = TRUE),
            commed = median(mediancompact, na.rm = TRUE),
            comsd = sd(mediancompact, na.rm = TRUE),
            commin = min(mediancompact, na.rm = TRUE),
            commax = max(mediancompact, na.rm = TRUE),
            rum = median(medrump, na.rm = TRUE),
            rumsd = sd(medrump, na.rm = TRUE),
            rummin = min(medrump, na.rm = TRUE),
            rummax = max(medrump, na.rm = TRUE),
            cc = median(medcc, na.rm = TRUE),
            ccsd = sd(medcc, na.rm = TRUE),
            ccmin = min(medcc, na.rm = TRUE),
            ccmax = max(medcc, na.rm = TRUE),
            mnhei = median(mnheight, na.rm = TRUE),
            mnheisd = sd(mnheight, na.rm = TRUE),
            mnheimin = min(mnheight, na.rm = TRUE),
            mnheimax = max(mnheight, na.rm = TRUE),
            mnneigh = median(medneigh, na.rm = TRUE),
            mnneighsd = sd(medneigh, na.rm = TRUE),
            mnneighmin = min(medneigh, na.rm = TRUE),
            mnneighmax = max(medneigh, na.rm = TRUE))
annualsummaries

annualsummarieslong <- annualsummaries %>% tidyr::pivot_longer(cols = c(area_mean:mnneighmax), names_to = "covariate", values_to = "value")
write.csv(annualsummarieslong, "./KlamathPlateau/HRSummaries/AnnualHomeRange_IndividualHRSizeLandcoverlong.csv", row.names = FALSE)

##### rand vs used summaries #####
options(digits=2)
annualsummaries <- indsums %>%
  group_by(Sex, contour) %>%
  summarize(fishers = n_distinct(fisher),
            Deployments = n_distinct(HR_ID),
            area_mean = mean(areakm),
            area_median = median(areakm),
            area_min = min(areakm),
            area_max = max(areakm),
            area_sd = sd(areakm),
            patctmed = median(ptchcnt),
            patctmin = min(ptchcnt),
            patctmax = max(ptchcnt),
            patctsd = sd(ptchcnt),
            medpatszmed = median(mdnptchs),
            medpatszmin = min(mdnptchs),
            medpatszmax = max(mdnptchs),
            medpatszsd = sd(mdnptchs),
            medpatlenmed = median(mdnptchl),
            medpatlenmin = min(mdnptchl),
            medpatlenmax = max(mdnptchl),
            medpatlensd = sd(mdnptchl),
            medasym_median = median(mednsym),
            medasym_min = min(mednsym),
            medasym_max = max(mednsym),
            medasym_sd = sd(mednsym),
            medcommed = median(mdncmpc),
            medcommin = min(mdncmpc),
            medcommax = max(mdncmpc),
            medcomsd = sd(mdncmpc),
            medrum = median(medrump),
            medrummin = min(medrump),
            medrummax = max(medrump),
            medrumsd = sd(medrump),
            medrum = median(medrump),
            medrummin = min(medrump),
            medrummax = max(medrump),
            medrumsd = sd(medrump),
            medmnhei = median(mnheght),
            medmnheimin = min(mnheght),
            medmnheimax = max(mnheght),
            medmnheisd = sd(mnheght))

write.csv(annualsummaries, "./KlamathPlateau/HRSummaries/AnnualHomeRange_SexdBBMMSizeLandcover.csv", row.names = FALSE)
