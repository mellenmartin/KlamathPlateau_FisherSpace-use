#### Klamath Plateau fisher space-use ####
#### Step 6 -- use formatted movement data to delineate corridor vs. non-corridor tracks ####
#### then buffer and extract veg data in corridor vs. non-corridor path segments
#### code based on vignettes from Movebank and Marco Smolla, Bart Kranstauber & Anne Scharf
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####

rm(list=ls());gc() #clear the memory

# set working drive
setwd("./KlamathPlateau/")

library(dplyr)
library(tidyr)
library(move)
library(lubridate)

#load the data
fisher <- read.csv("KlamathPlateau_FisherSSFLocs_UsedandRandom.csv",header=T)
head(fisher)
fisher <- fisher %>% filter(case_ == "TRUE")
fisher$fisher = substr(fisher$uniqueID, 1, 4)

#specify the format of the date/time label and the time zone, and convert to POSIX
fisher$datetime <- as_datetime(fisher$t1_)

#specify the format of the date/time label and the time zone, and convert to POSIX
fisher <- fisher %>% dplyr::select(uniqueID,
                            fisher,
                            burst_,
                            x1_,
                            x2_,
                            y1_,
                            y2_,
                            sl_,
                            ta_,
                            t1_,
                            t2_,
                            dt_,
                            case_,
                            log_sl_,
                            cos_ta_,
                            y,
                            step_id,
                            step_id2,
                            datetime) %>% drop_na(x1_|y1_|x2_|y2_|datetime)

fisher$group <- c(1:length(fisher$burst_))
fisher$acID <- paste(fisher$uniqueID, fisher$group, sep = "_")

fisher_move <- move::move(x=fisher$x1_, y=fisher$y1_, 
              time=fisher$datetime, 
              proj=CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84"), 
              data=fisher, animal=fisher$uniqueID, acID = acID)

fisher_move_latlong <- spTransform(fisher_move, CRS("+proj=longlat +datum=WGS84"))

library(moveVis)
library(move)
##### corridors #####

F01TF15_corr <- move::corridor(fisher_move_latlong[[1]], speedProp = 0.75, circProp = 0.25)
plot(F01TF15_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")


F01TS16a_corr <- corridor(fisher_move_latlong[[2]], speedProp = 0.75, circProp = 0.25)
plot(F01TS16a_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F01TS16b_corr <- corridor(fisher_move_latlong[[3]], speedProp = 0.75, circProp = 0.25)
plot(F01TS16b_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F01TW17_corr <- corridor(fisher_move_latlong[[4]], speedProp = 0.75, circProp = 0.25)
plot(F01TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F01TS18_corr <- corridor(fisher_move_latlong[[5]], speedProp = 0.75, circProp = 0.25)
plot(F01TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F02TF15_corr <- corridor(fisher_move_latlong[[6]], speedProp = 0.75, circProp = 0.25)
plot(F02TF15_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F02TS16b_corr <- corridor(fisher_move_latlong[[7]], speedProp = 0.75, circProp = 0.25)
plot(F02TS16b_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F03TF15_corr <- corridor(fisher_move_latlong[[8]], speedProp = 0.75, circProp = 0.25)
plot(F03TF15_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F03TS16a_corr <- corridor(fisher_move_latlong[[9]], speedProp = 0.75, circProp = 0.25)
plot(F03TS16a_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F03TS16b_corr <- corridor(fisher_move_latlong[[10]], speedProp = 0.75, circProp = 0.25)
plot(F03TS16b_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F03TW17_corr <- corridor(fisher_move_latlong[[11]], speedProp = 0.75, circProp = 0.25)
plot(F03TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F07TF16_corr <- corridor(fisher_move_latlong[[12]], speedProp = 0.75, circProp = 0.25)
plot(F07TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F07TF17_corr <- corridor(fisher_move_latlong[[13]], speedProp = 0.75, circProp = 0.25)
plot(F07TF17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

F07TS18_corr <- corridor(fisher_move_latlong[[14]], speedProp = 0.75, circProp = 0.25)
plot(F07TS18_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M03TF16_corr <- corridor(fisher_move_latlong[[15]], speedProp = 0.75, circProp = 0.25)
plot(M03TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M03TW17_corr <- corridor(fisher_move_latlong[[16]], speedProp = 0.75, circProp = 0.25)
plot(M03TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M03TF17_corr <- corridor(fisher_move_latlong[[17]], speedProp = 0.75, circProp = 0.25)
plot(M03TF17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M05TF16_corr <- corridor(fisher_move_latlong[[18]], speedProp = 0.75, circProp = 0.25)
plot(M05TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M07TF16_corr <- corridor(fisher_move_latlong[[19]], speedProp = 0.75, circProp = 0.25)
plot(M07TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M08TF16_corr <- corridor(fisher_move_latlong[[20]], speedProp = 0.75, circProp = 0.25)
plot(M08TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M08TW17_corr <- corridor(fisher_move_latlong[[21]], speedProp = 0.75, circProp = 0.25)
plot(M08TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M08TF17_corr <- corridor(fisher_move_latlong[[22]], speedProp = 0.75, circProp = 0.25)
plot(M08TF17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M08TS18_corr <- corridor(fisher_move_latlong[[23]], speedProp = 0.75, circProp = 0.25)
plot(M08TF16_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M09TW17_corr <- corridor(fisher_move_latlong[[24]], speedProp = 0.75, circProp = 0.25)
plot(M09TW17_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")

M09TS18_corr <- corridor(fisher_move_latlong[[25]], speedProp = 0.75, circProp = 0.25)
plot(M09TS18_corr, type="l", xlab="Easting", ylab="Northing", col=c("black", "grey"), lwd=c(1,2))
legend("bottomleft", c("Corridor", "Non-corridor"), col=c("black", "grey"), lty=c(1,1), bty="n")


#### write and export the line data
library(readr)
library(sf)
library(tidyverse)

F01TF15_data <- as.data.frame(fisher_move_latlong[[1]])
F01TF15_data$corridor <- c(2,F01TF15_corr@burstId)
F01TF15_data$corridortext <- ifelse(F01TF15_data$corridor == 1, "corridor", "non-corridor")

F01TS16a_data <- as.data.frame(fisher_move_latlong[[2]])
F01TS16a_data$corridor <- c(2,F01TS16a_corr@burstId)
F01TS16a_data$corridortext <- ifelse(F01TS16a_data$corridor == 1, "corridor", "non-corridor")

F01TS16b_data <- as.data.frame(fisher_move_latlong[[3]])
F01TS16b_data$corridor <- c(2,F01TS16b_corr@burstId)
F01TS16b_data$corridortext <- ifelse(F01TS16b_data$corridor == 1, "corridor", "non-corridor")

F01TW17_data <- as.data.frame(fisher_move_latlong[[4]])
F01TW17_data$corridor <- c(2,F01TW17_corr@burstId)
F01TW17_data$corridortext <- ifelse(F01TW17_data$corridor == 1, "corridor", "non-corridor")

F01TS18_data <- as.data.frame(fisher_move_latlong[[5]])
F01TS18_data$corridor <- c(2,F01TS18_corr@burstId)
F01TS18_data$corridortext <- ifelse(F01TS18_data$corridor == 1, "corridor", "non-corridor")

F02TF15_data <- as.data.frame(fisher_move_latlong[[6]])
F02TF15_data$corridor <- c(2,F02TF15_corr@burstId)
F02TF15_data$corridortext <- ifelse(F02TF15_data$corridor == 1, "corridor", "non-corridor")

F02TS16b_data <- as.data.frame(fisher_move_latlong[[7]])
F02TS16b_data$corridor <- c(2,F02TS16b_corr@burstId)
F02TS16b_data$corridortext <- ifelse(F02TS16b_data$corridor == 1, "corridor", "non-corridor")

F03TF15_data <- as.data.frame(fisher_move_latlong[[8]])
F03TF15_data$corridor <- c(2,F03TF15_corr@burstId)
F03TF15_data$corridortext <- ifelse(F03TF15_data$corridor == 1, "corridor", "non-corridor")

F03TS16a_data <- as.data.frame(fisher_move_latlong[[9]])
F03TS16a_data$corridor <- c(2,F03TS16a_corr@burstId)
F03TS16a_data$corridortext <- ifelse(F03TS16a_data$corridor == 1, "corridor", "non-corridor")

F03TS16b_data <- as.data.frame(fisher_move_latlong[[10]])
F03TS16b_data$corridor <- c(2,F03TS16b_corr@burstId)
F03TS16b_data$corridortext <- ifelse(F03TS16b_data$corridor == 1, "corridor", "non-corridor")

F03TW17_data <- as.data.frame(fisher_move_latlong[[11]])
F03TW17_data$corridor <- c(2,F03TW17_corr@burstId)
F03TW17_data$corridortext <- ifelse(F03TW17_data$corridor == 1, "corridor", "non-corridor")

F07TF16_data <- as.data.frame(fisher_move_latlong[[12]])
F07TF16_data$corridor <- c(2,F07TF16_corr@burstId)
F07TF16_data$corridortext <- ifelse(F07TF16_data$corridor == 1, "corridor", "non-corridor")

F07TF17_data <- as.data.frame(fisher_move_latlong[[13]])
F07TF17_data$corridor <- c(2,F07TF17_corr@burstId)
F07TF17_data$corridortext <- ifelse(F07TF17_data$corridor == 1, "corridor", "non-corridor")

F07TS18_data <- as.data.frame(fisher_move_latlong[[14]])
F07TS18_data$corridor <- c(2,F07TS18_corr@burstId)
F07TS18_data$corridortext <- ifelse(F07TS18_data$corridor == 1, "corridor", "non-corridor")

M03TF16_data <- as.data.frame(fisher_move_latlong[[15]])
M03TF16_data$corridor <- c(2,M03TF16_corr@burstId)
M03TF16_data$corridortext <- ifelse(M03TF16_data$corridor == 1, "corridor", "non-corridor")

M03TW17_data <- as.data.frame(fisher_move_latlong[[16]])
M03TW17_data$corridor <- c(2,M03TW17_corr@burstId)
M03TW17_data$corridortext <- ifelse(M03TW17_data$corridor == 1, "corridor", "non-corridor")

M03TF17_data <- as.data.frame(fisher_move_latlong[[17]])
M03TF17_data$corridor <- c(2,M03TF17_corr@burstId)
M03TF17_data$corridortext <- ifelse(M03TF17_data$corridor == 1, "corridor", "non-corridor")

M05TF16_data <- as.data.frame(fisher_move_latlong[[18]])
M05TF16_data$corridor <- c(2,M05TF16_corr@burstId)
M05TF16_data$corridortext <- ifelse(M05TF16_data$corridor == 1, "corridor", "non-corridor")

M07TF16_data <- as.data.frame(fisher_move_latlong[[19]])
M07TF16_data$corridor <- c(2,M07TF16_corr@burstId)
M07TF16_data$corridortext <- ifelse(M07TF16_data$corridor == 1, "corridor", "non-corridor")

M08TF16_data <- as.data.frame(fisher_move_latlong[[20]])
M08TF16_data$corridor <- c(2,M08TF16_corr@burstId)
M08TF16_data$corridortext <- ifelse(M08TF16_data$corridor == 1, "corridor", "non-corridor")

M08TW17_data <- as.data.frame(fisher_move_latlong[[21]])
M08TW17_data$corridor <- c(2,M08TW17_corr@burstId)
M08TW17_data$corridortext <- ifelse(M08TW17_data$corridor == 1, "corridor", "non-corridor")

M08TF17_data <- as.data.frame(fisher_move_latlong[[22]])
M08TF17_data$corridor <- c(2,M08TF17_corr@burstId)
M08TF17_data$corridortext <- ifelse(M08TF17_data$corridor == 1, "corridor", "non-corridor")

M08TS18_data <- as.data.frame(fisher_move_latlong[[23]])
M08TS18_data$corridor <- c(2,M08TS18_corr@burstId)
M08TS18_data$corridortext <- ifelse(M08TS18_data$corridor == 1, "corridor", "non-corridor")

M09TW17_data <- as.data.frame(fisher_move_latlong[[24]])
M09TW17_data$corridor <- c(2,M09TW17_corr@burstId)
M09TW17_data$corridortext <- ifelse(M09TW17_data$corridor == 1, "corridor", "non-corridor")

M09TS18_data <- as.data.frame(fisher_move_latlong[[25]])
M09TS18_data$corridor <- c(2,M09TS18_corr@burstId)
M09TS18_data$corridortext <- ifelse(M09TS18_data$corridor == 1, "corridor", "non-corridor")

####
corrdata <- rbind(F01TF15_data,
                      F01TS16a_data,
                      F01TS16b_data,
                      F01TS18_data,
                      F01TW17_data,
                      F02TF15_data,
                      F02TS16b_data,
                      F03TF15_data,
                      F03TS16a_data,
                      F03TS16b_data,
                      F03TW17_data,
                      F07TF16_data,
                      F07TF17_data,
                      F07TS18_data,
                      M03TF16_data,
                      M03TF17_data,
                      M03TW17_data,
                      M05TF16_data,
                      M07TF16_data,
                      M08TF16_data,
                      M08TF17_data,
                      M08TW17_data,
                      M08TS18_data,
                      M09TW17_data,
                      M09TS18_data)


F01TF15_start_end_xy <- F01TF15_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F01TF15_start_end_xy, seq(nrow(F01TF15_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F01TF15_lines <- st_sfc(lines)
F01TF15_lines_sf <- st_sf('acID' = F01TF15_start_end_xy$acID, 'geometry' = lines)
F01TF15_lines_sf <- inner_join(F01TF15_lines_sf, F01TF15_data, keep = FALSE)
st_write(F01TF15_lines_sf, "Corridors/KlamathPlateau_F01TF15_lines_230411.shp")

F01TS16a_start_end_xy <- F01TS16a_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F01TS16a_start_end_xy, seq(nrow(F01TS16a_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F01TS16a_lines <- st_sfc(lines)
F01TS16a_lines_sf <- st_sf('acID' = F01TS16a_start_end_xy$acID, 'geometry' = lines)
F01TS16a_lines_sf <- inner_join(F01TS16a_lines_sf, F01TS16a_data, keep = FALSE)
st_write(F01TS16a_lines_sf, "Corridors/KlamathPlateau_F01TS16a_lines_230411.shp")

F01TS16b_start_end_xy <- F01TS16b_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F01TS16b_start_end_xy, seq(nrow(F01TS16b_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F01TS16b_lines <- st_sfc(lines)
F01TS16b_lines_sf <- st_sf('acID' = F01TS16b_start_end_xy$acID, 'geometry' = lines)
F01TS16b_lines_sf <- inner_join(F01TS16b_lines_sf, F01TS16b_data, keep = FALSE)
st_write(F01TS16b_lines_sf, "Corridors/KlamathPlateau_F01TS16b_lines_230411.shp")

F01TW17_start_end_xy <- F01TW17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F01TW17_start_end_xy, seq(nrow(F01TW17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F01TW17_lines <- st_sfc(lines)
F01TW17_lines_sf <- st_sf('acID' = F01TW17_start_end_xy$acID, 'geometry' = lines)
F01TW17_lines_sf <- inner_join(F01TW17_lines_sf, F01TW17_data, keep = FALSE)
st_write(F01TW17_lines_sf, "Corridors/KlamathPlateau_F01TW17_lines_230411.shp")

F01TS18_start_end_xy <- F01TS18_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F01TS18_start_end_xy, seq(nrow(F01TS18_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F01TS18_lines <- st_sfc(lines)
F01TS18_lines_sf <- st_sf('acID' = F01TS18_start_end_xy$acID, 'geometry' = lines)
F01TS18_lines_sf <- inner_join(F01TS18_lines_sf, F01TS18_data, keep = FALSE)
st_write(F01TS18_lines_sf, "Corridors/KlamathPlateau_F01TS18_lines_230411.shp")

F02TF15_start_end_xy <- F02TF15_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F02TF15_start_end_xy, seq(nrow(F02TF15_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F02TF15_lines <- st_sfc(lines)
F02TF15_lines_sf <- st_sf('acID' = F02TF15_start_end_xy$acID, 'geometry' = lines)
F02TF15_lines_sf <- inner_join(F02TF15_lines_sf, F02TF15_data, keep = FALSE)
st_write(F02TF15_lines_sf, "Corridors/KlamathPlateau_F02TF15_lines_230411.shp")

F02TS16b_start_end_xy <- F02TS16b_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F02TS16b_start_end_xy, seq(nrow(F02TS16b_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F02TS16b_lines <- st_sfc(lines)
F02TS16b_lines_sf <- st_sf('acID' = F02TS16b_start_end_xy$acID, 'geometry' = lines)
F02TS16b_lines_sf <- inner_join(F02TS16b_lines_sf, F02TS16b_data, keep = FALSE)
st_write(F02TS16b_lines_sf, "Corridors/KlamathPlateau_F02TS16b_lines_230411.shp")

F03TF15_start_end_xy <- F03TF15_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F03TF15_start_end_xy, seq(nrow(F03TF15_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F03TF15_lines <- st_sfc(lines)
F03TF15_lines_sf <- st_sf('acID' = F03TF15_start_end_xy$acID, 'geometry' = lines)
F03TF15_lines_sf <- inner_join(F03TF15_lines_sf, F03TF15_data, keep = FALSE)
st_write(F03TF15_lines_sf, "Corridors/KlamathPlateau_F03TF15_lines_230411.shp")

F03TS16a_start_end_xy <- F03TS16a_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F03TS16a_start_end_xy, seq(nrow(F03TS16a_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F03TS16a_lines <- st_sfc(lines)
F03TS16a_lines_sf <- st_sf('acID' = F03TS16a_start_end_xy$acID, 'geometry' = lines)
F03TS16a_lines_sf <- inner_join(F03TS16a_lines_sf, F03TS16a_data, keep = FALSE)
st_write(F03TS16a_lines_sf, "Corridors/KlamathPlateau_F03TS16a_lines_230411.shp")

F03TS16b_start_end_xy <- F03TS16b_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F03TS16b_start_end_xy, seq(nrow(F03TS16b_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F03TS16b_lines <- st_sfc(lines)
F03TS16b_lines_sf <- st_sf('acID' = F03TS16b_start_end_xy$acID, 'geometry' = lines)
F03TS16b_lines_sf <- inner_join(F03TS16b_lines_sf, F03TS16b_data, keep = FALSE)
st_write(F03TS16b_lines_sf, "Corridors/KlamathPlateau_F03TS16b_lines_230411.shp")

F03TW17_start_end_xy <- F03TW17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F03TW17_start_end_xy, seq(nrow(F03TW17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F03TW17_lines <- st_sfc(lines)
F03TW17_lines_sf <- st_sf('acID' = F03TW17_start_end_xy$acID, 'geometry' = lines)
F03TW17_lines_sf <- inner_join(F03TW17_lines_sf, F03TW17_data, keep = FALSE)
st_write(F03TW17_lines_sf, "Corridors/KlamathPlateau_F03TW17_lines_230411.shp")

F07TF16_start_end_xy <- F07TF16_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F07TF16_start_end_xy, seq(nrow(F07TF16_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F07TF16_lines <- st_sfc(lines)
F07TF16_lines_sf <- st_sf('acID' = F07TF16_start_end_xy$acID, 'geometry' = lines)
F07TF16_lines_sf <- inner_join(F07TF16_lines_sf, F07TF16_data, keep = FALSE)
st_write(F07TF16_lines_sf, "Corridors/KlamathPlateau_F07TF16_lines_230411.shp")

F07TF17_start_end_xy <- F07TF17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F07TF17_start_end_xy, seq(nrow(F07TF17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F07TF17_lines <- st_sfc(lines)
F07TF17_lines_sf <- st_sf('acID' = F07TF17_start_end_xy$acID, 'geometry' = lines)
F07TF17_lines_sf <- inner_join(F07TF17_lines_sf, F07TF17_data, keep = FALSE)
st_write(F07TF17_lines_sf, "Corridors/KlamathPlateau_F07TF17_lines_220411.shp")

F07TS18_start_end_xy <- F07TS18_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(F07TS18_start_end_xy, seq(nrow(F07TS18_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
F07TS18_lines <- st_sfc(lines)
F07TS18_lines_sf <- st_sf('acID' = F07TS18_start_end_xy$acID, 'geometry' = lines)
F07TS18_lines_sf <- inner_join(F07TS18_lines_sf, F07TS18_data, keep = FALSE)
st_write(F07TS18_lines_sf, "Corridors/KlamathPlateau_F07TS18_lines_230411.shp")

M03TF16_start_end_xy <- M03TF16_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M03TF16_start_end_xy, seq(nrow(M03TF16_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M03TF16_lines <- st_sfc(lines)
M03TF16_lines_sf <- st_sf('acID' = M03TF16_start_end_xy$acID, 'geometry' = lines)
M03TF16_lines_sf <- inner_join(M03TF16_lines_sf, M03TF16_data, keep = FALSE)
st_write(M03TF16_lines_sf, "Corridors/KlamathPlateau_M03TF16_lines_230411.shp")

M03TF17_start_end_xy <- M03TF17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M03TF17_start_end_xy, seq(nrow(M03TF17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M03TF17_lines <- st_sfc(lines)
M03TF17_lines_sf <- st_sf('acID' = M03TF17_start_end_xy$acID, 'geometry' = lines)
M03TF17_lines_sf <- inner_join(M03TF17_lines_sf, M03TF17_data, keep = FALSE)
st_write(M03TF17_lines_sf, "Corridors/KlamathPlateau_M03TF17_lines_230411.shp")

M03TW17_start_end_xy <- M03TW17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M03TW17_start_end_xy, seq(nrow(M03TW17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M03TW17_lines <- st_sfc(lines)
M03TW17_lines_sf <- st_sf('acID' = M03TW17_start_end_xy$acID, 'geometry' = lines)
M03TW17_lines_sf <- inner_join(M03TW17_lines_sf, M03TW17_data, keep = FALSE)
st_write(M03TW17_lines_sf, "Corridors/KlamathPlateau_M03TW17_lines_230411.shp")

M05TF16_start_end_xy <- M05TF16_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M05TF16_start_end_xy, seq(nrow(M05TF16_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M05TF16_lines <- st_sfc(lines)
M05TF16_lines_sf <- st_sf('acID' = M05TF16_start_end_xy$acID, 'geometry' = lines)
M05TF16_lines_sf <- inner_join(M05TF16_lines_sf, M05TF16_data, keep = FALSE)
st_write(M05TF16_lines_sf, "Corridors/KlamathPlateau_M05TF16_lines_230411.shp")

M07TF16_start_end_xy <- M07TF16_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M07TF16_start_end_xy, seq(nrow(M07TF16_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M07TF16_lines <- st_sfc(lines)
M07TF16_lines_sf <- st_sf('acID' = M07TF16_start_end_xy$acID, 'geometry' = lines)
M07TF16_lines_sf <- inner_join(M07TF16_lines_sf, M07TF16_data, keep = FALSE)
st_write(M07TF16_lines_sf, "Corridors/KlamathPlateau_M07TF16_lines_230411.shp")

M08TF16_start_end_xy <- M08TF16_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M08TF16_start_end_xy, seq(nrow(M08TF16_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M08TF16_lines <- st_sfc(lines)
M08TF16_lines_sf <- st_sf('acID' = M08TF16_start_end_xy$acID, 'geometry' = lines)
M08TF16_lines_sf <- inner_join(M08TF16_lines_sf, M08TF16_data, keep = FALSE)
st_write(M08TF16_lines_sf, "Corridors/KlamathPlateau_M08TF16_lines_230411.shp")

M08TF17_start_end_xy <- M08TF17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M08TF17_start_end_xy, seq(nrow(M08TF17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M08TF17_lines <- st_sfc(lines)
M08TF17_lines_sf <- st_sf('acID' = M08TF17_start_end_xy$acID, 'geometry' = lines)
M08TF17_lines_sf <- inner_join(M08TF17_lines_sf, M08TF17_data, keep = FALSE)
st_write(M08TF17_lines_sf, "Corridors/KlamathPlateau_M08TF17_lines_230411.shp")


M08TW17_start_end_xy <- M08TW17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M08TW17_start_end_xy, seq(nrow(M08TW17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M08TW17_lines <- st_sfc(lines)
M08TW17_lines_sf <- st_sf('acID' = M08TW17_start_end_xy$acID, 'geometry' = lines)
M08TW17_lines_sf <- inner_join(M08TW17_lines_sf, M08TW17_data, keep = FALSE)
st_write(M08TW17_lines_sf, "Corridors/KlamathPlateau_M08TW17_lines_230411.shp")

M08TS18_start_end_xy <- M08TS18_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M08TS18_start_end_xy, seq(nrow(M08TS18_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M08TS18_lines <- st_sfc(lines)
M08TS18_lines_sf <- st_sf('acID' = M08TS18_start_end_xy$acID, 'geometry' = lines)
M08TS18_lines_sf <- inner_join(M08TS18_lines_sf, M08TS18_data, keep = FALSE)
st_write(M08TS18_lines_sf, "Corridors/KlamathPlateau_M08TS18_lines_230411.shp")

M09TW17_start_end_xy <- M09TW17_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M09TW17_start_end_xy, seq(nrow(M09TW17_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M09TW17_lines <- st_sfc(lines)
M09TW17_lines_sf <- st_sf('acID' = M09TW17_start_end_xy$acID, 'geometry' = lines)
M09TW17_lines_sf <- inner_join(M09TW17_lines_sf, M09TW17_data, keep = FALSE)
st_write(M09TW17_lines_sf, "Corridors/KlamathPlateau_M09TW17_lines_230411.shp")

M09TS18_start_end_xy <- M09TS18_data %>% dplyr::select(acID, x1_, y1_, x2_, y2_, corridor) %>% drop_na(x1_, y1_, x2_, y2_)
rows <- split(M09TS18_start_end_xy, seq(nrow(M09TS18_start_end_xy)))
lines <- lapply(rows, function(row) {
  lmat <- matrix(unlist(row[2:5]), ncol = 2, byrow = TRUE)
  st_linestring(lmat)
})
M09TS18_lines <- st_sfc(lines)
M09TS18_lines_sf <- st_sf('acID' = M09TS18_start_end_xy$acID, 'geometry' = lines)
M09TS18_lines_sf <- inner_join(M09TS18_lines_sf, M09TS18_data, keep = FALSE)
st_write(M09TS18_lines_sf, "Corridors/KlamathPlateau_M09TS18_lines_230411.shp")

sf_lines_list <- list(F01TF15_lines_sf,
                      F01TS16a_lines_sf,
                      F01TS16b_lines_sf,
                      F01TS18_lines_sf,
                      F01TW17_lines_sf,
                      F02TF15_lines_sf,
                      F02TS16b_lines_sf,
                      F03TF15_lines_sf,
                      F03TS16a_lines_sf,
                      F03TS16b_lines_sf,
                      F03TW17_lines_sf,
                      F07TF16_lines_sf,
                      F07TF17_lines_sf,
                      F07TS18_lines_sf,
                      M03TF16_lines_sf,
                      M03TF17_lines_sf,
                      M03TW17_lines_sf,
                      M05TF16_lines_sf,
                      M07TF16_lines_sf,
                      M08TF16_lines_sf,
                      M08TF17_lines_sf,
                        M08TW17_lines_sf,
                        M08TS18_lines_sf,
                        M09TW17_lines_sf,
                        M09TS18_lines_sf)

library(sf)
library(data.table)
KlamathPlateau_alllines <- do.call(what = sf:::rbind.sf, args = sf_lines_list)
st_write(KlamathPlateau_alllines, "Corridors/KlamathPlateau_AllFishers_lines_230411.shp")

#### buffer and extract covariate values in paths/path segments ####
library(stars) 
library(ggnewscale) 
library(geobgu) # install from GitHub ("michaeldorman/geobgu")
library(dplyr) 
library(sf) 

wd <- "./KlamathPlateau/Corridors"
setwd(wd)

# load in HR contours
filenames <- list.files(wd, pattern = "*.shp$")

# Create list to store the movement paths within the loop
lines.list <- list()

## Import all bbmm shapefiles, add individual identifiers, calculate area, and calculate stream density
for(i in 2:length(filenames)){
  library(sf)
  wd <- "./KlamathPlateau/Corridors"
  setwd(wd)
  
  # set CRS for the shapefiles
  crs <-  "+proj=utm +zone=10 +ellps=WGS84 +datum=NAD83 units=m"
  
  temp.sf <- st_read(filenames[i])
  buffer <- mean(temp.sf$sl_)
  temp.sf.lines <- st_sf(temp.sf)
  st_crs(temp.sf.lines) <- 26910
  
  polyname <- sub("*.shp", "", filenames[i]) 
  temp.sf.lines <- st_buffer(temp.sf.lines, dist = buffer)
  temp.sf.lines$areakm <- as.numeric(st_area(temp.sf.lines)/1000000)
  
  ecog <- st_read("./KlamathPlateau/rasters/Ecognition_CoreData_220731.shp")
  ecog <- st_buffer(ecog, dist=0)
  temp.sf_ecog <- st_crop(ecog, temp.sf.lines) 
  
  library(stars)
  library(geobgu)
  # read rasters in
  ddi <- read_stars("./KlamathPlateau/rasters/ddi_reproject.tif")
  cc <- read_stars("./KlamathPlateau/rasters/cc_reproject.tif")
  cclayers <- read_stars("./KlamathPlateau/rasters/cclayers_reproject.tif")
  basal <- read_stars("./KlamathPlateau/rasters/basal_reproject.tif")
  height <- read_stars("./KlamathPlateau/rasters/height_reproject.tif") 
  qmd <- read_stars("./KlamathPlateau/rasters/qmd_reproject.tif") 
  
  st_crs(ddi) <- 26910
  st_crs(cc) <- 26910
  st_crs(cclayers) <- 26910
  st_crs(basal) <- 26910
  st_crs(height) <- 26910
  st_crs(qmd)  <- 26910
  
  ecog_covs <- temp.sf.lines %>%
    st_join(temp.sf_ecog) %>%
    group_by(acID) %>% 
    summarize(n_patches = n(),
              medianpatchsize = median(Shape_Area),
              medianpatchlength = median(Shape_Leng),
              medianasym = median(Asymmetry),
              mediancompact = median(Compactnes),
              medianbright = median(Brightness),
              medianmaxdiff = median(Max_diff),
              medianMnDiff_p95 = median(MnDiff_p95),
              medianMnDiff_ele = median(MnDiff_ele),
              medianMn_CanRelR = median(Mn_CanRelR),
              medianMn_CanRump = median(Mn_CanRump),
              medianMn_Kurtosi = median(Mn_Kurtosi),
              medianSD_CanRump = median(SD_CanRump))
  
  gnn_covs <- temp.sf.lines %>% 
    mutate(ddicov = raster_extract(ddi, temp.sf.lines, fun = mean, na.rm = TRUE),
           cccov = raster_extract(cc, temp.sf.lines, fun = mean, na.rm = TRUE),
           cclayerscov = raster_extract(cclayers, temp.sf.lines, fun = mean, na.rm = TRUE),
           basalcov = raster_extract(basal, temp.sf.lines, fun = mean, na.rm = TRUE),
           heightcov = raster_extract(height, temp.sf.lines, fun = mean, na.rm = TRUE),
           qmdcov = raster_extract(qmd, temp.sf.lines, fun = mean, na.rm = TRUE))
  
  #temp.sf.lines <- inner_join(temp.sf.lines %>% as.data.frame(), stream_covs %>% as.data.frame(), by = "acID")
  temp.sf.lines1 <- left_join(temp.sf.lines %>% as.data.frame(), ecog_covs %>% as.data.frame(), by = "acID")
  temp.sf.lines2 <- left_join(temp.sf.lines1 %>% as.data.frame(), gnn_covs %>% as.data.frame(), by = "acID")
  
  temp.sf.lines <- st_sf(temp.sf.lines2)
  # write the combined HR polygons w/ stream info for future use
  setwd("./KlamathPlateau/Corridors/")
  #temp.sf.lines <- inner_join(temp.sf.lines, temp.sf, keep = TRUE)
  sf::st_write(temp.sf.lines, paste0(polyname,"KlamathPlateau_lines_HRcovs_230411.shp"))
  
  #### add to list ####
  temp <- list(temp.sf.lines)
  # add to list
  lines.list[[i]] <- temp
  
}

wd <- "./KlamathPlateau/Corridors/"
setwd(wd)

# load in HR contours
filenames <- list.files(wd, pattern = "*.shp$")

library(sf)
library(data.table)

# load in HR contours
filenames <- list.files(wd, pattern = "*.shp$")
sffiles <- lapply(filenames, st_read)
combineShp <- mapedit:::combine_list_of_sf(sffiles)
combineShp <- do.call(what = sf:::rbind.sf, args = sffiles)
#combineShp <- rbind(lines.list)
st_crs(combineShp) <- 26910

# write the combined HR polygons w/ stream info for future use
st_write(combineShp, "./KlamathPlateau/Corridors/KlamathPlateau_allfisherlines_HRcovs_230411.shp")

# turn shapefile to dataframe 
corridorcovariates <- as.data.frame(combineShp)
alldata <- corridorcovariates
st_geometry(alldata) <- NULL

# format data 
alldata$corridorbin <- ifelse(alldata$crrdr_x == 2 , 0, 1)
alldata$basalcv <- alldata$basalcv/100
alldata$cccov <- alldata$cccov/100
alldata$ddicov <- alldata$ddicov/100
alldata$heghtcv <- alldata$heghtcv/100
alldata$qmdcov <- alldata$qmdcov/10
write.csv(alldata, file = "KlamathPlateau_allfisherlines_HRcovs.csv" , row.names = FALSE)