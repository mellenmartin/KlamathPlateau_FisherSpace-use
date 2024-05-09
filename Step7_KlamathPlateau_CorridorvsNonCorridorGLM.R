#### Klamath Plateau fisher space-use ####
#### Step 7 --   glms to understand effects of stand structure on corridor vs. non-corridor movement
#### Martin et al. in review -- forest structure affects fisher space-use ####
#### data collected by Moriarty, Matthews, et al. ####

rm(list=ls());gc() #clear the memory

library(ResourceSelection)
library(glmmTMB)

### load data ####
alldata <- read.csv(".KlamathPlateau/Corridors/KlamathPlateau_allfisherlines_HRcovs.csv")

with(alldata, prop.table(table(unqID_y, corridorbin), 1))

alldata$weight <- 1000^(1-alldata$corridorbin)
alldata$mdnptchs <- alldata$mdnptchs/1000000 # median patch size to km2
alldata$strata <- paste0(alldata$unqID_y,"-",alldata$brst__y) # create strata

# scale covariates
alldata$snpatches <- scale(alldata$n_ptchs)
alldata$smdnpatchsize <- scale(alldata$mdnptchs)
alldata$smdnpatchlen <- scale(alldata$mdnptchl)
alldata$smdnasym <- scale(alldata$mednsym)
alldata$smdncanrum <- scale(alldata$mdM_CRR)
alldata$smeandiffhgt <- scale(alldata$mdMD_95)

#### corridor rsf -- how do paths classified as corridors vary from non-corridors ####
corridorselection <- glmmTMB(corridorbin ~ 
                               snpatches + 
                               smdnpatchsize +
                               smdnpatchlen +
                               smdnasym +  
                               smdncanrum +
                               smeandiffhgt, 
                               family=poisson(), weights = weight, data = alldata)

summary(corridorselection3)
tabranef(corridorselection)
fixef(corridorselection)