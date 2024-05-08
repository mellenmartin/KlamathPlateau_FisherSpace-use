#### Klamath Plateau fisher space-use ####
#### Step 5 -- fitting SSFs with GLMMTMB for main manuscript results ####
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####

rm(list=ls());gc() #clear the memory

##### fit SSFs with GLMMTMB ####
library(glmmTMB)
library(INLA)
library(tidyverse)
library(raster)
library(survival)
library(TwoStepCLogit)
library(amt)
library(sf)
dat_ssf <- read.csv("./KlamathPlateau/KlamathPlateau_FisherSSFLocs_UsedandRandom.csv")

dat_ssf$fisher = substr(dat_ssf$uniqueID, 1, 4)
dat_ssf$fishernum <- as.numeric(unique(dat_ssf$fisher))
dat_ssf$julian <- yday(dat_ssf$t2_)
dat_ssf$season <- ifelse(dat_ssf$julian>79 & dat_ssf$julian<264, "Summer", "Winter")
dat_ssf$step_id3 <- paste0(as.numeric(factor(dat_ssf$fishernum)),"-", dat_ssf$burst_)
randlocs <- dat_ssf %>% dplyr::filter(case_ == "FALSE")
usedlocs <- dat_ssf %>% dplyr::filter(case_ != "FALSE")

### w/ glmmTMB
##### now fit the model with multiple covariates #####
dat_ssf <- dat_ssf %>% dplyr::filter(uniqueID != "M08T-Spring18") 
dat_ssfsf <- dat_ssf %>% filter(!x2_ == "") %>% 
  st_as_sf(., coords = c("x2_","y2_"), crs = 26910)

# weight the available locations
dat_ssf$weight <- 10000^(1-dat_ssf$y)
set.seed(516)
TMBStruc <- glmmTMB(y ~ -1 + scc_end 
                    + I(scc_end^2)
                    + sheight_end
                    + I(sheight_end^2)
                    + sdistedge
                    + sccdiff_end
                    + shgtdiff_end
                    + (1|step_id3) 
                    + (0 + scc_end | fisher)
                    + (0 + I(scc_end^2) | fisher)
                    + (0 + sheight_end| fisher)
                    + (0 + I(sheight_end^2)| fisher)
                    + (0 + sdistedge| fisher)
                    + (0 + sccdiff_end | fisher)
                    + (0 + shgtdiff_end| fisher),
                    family=poisson, weights = weight, data = dat_ssf, doFit = FALSE) 
#Set the value of the standard deviation of the first random effect (here (1|step_id)):
TMBStruc$parameters$theta[1] <- log(1e3) 

#Tell glmmTMB not to change the first standard deviation, all other values are freely estimated 
#(and are different from each other)
TMBStruc$mapArg = list(theta=factor(c(NA,1:7)))

#Fit the model and look at the summary:
glmm.TMB.random3 <- glmmTMB:::fitTMB(TMBStruc)

summary(glmm.TMB.random3)
