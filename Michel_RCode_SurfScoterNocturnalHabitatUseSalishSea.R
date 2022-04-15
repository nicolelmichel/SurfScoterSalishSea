###############################################################################
## R Code from "Surf Scoters use deeper offshore waters during nocturnal 
##     resting periods in the Salish Sea of Washington and British 
##     Columbia" published in 2022 in Ornithological Applications
##     Authors: Lindsey Hamilton, Nicole Michel, Joseph Evenson, Dina Roberts
##
## Code and analyses by Nicole Michel 
## Contact: Nicole.L.Michel1@gmail.com
###############################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### prep for measuring Euclidean least cost paths over water
# use igraph method from here: https://insileco.github.io/2019/04/08/r-as-a-ruler-how-to-calculate-distances-between-geographical-objects/
library(knitr)
library(sf)
library(sp)
library(igraph)
library(raster)
library(rgeos)

# load the hex grid (created with 100m sides in ArcMap)
ss.hexgrid <- st_read("SalishSeaHexGrid_100mClip.shp")

# create an edgelist (i.e. list to 'from' and 'to' cell id's that represent neighbours)
edgelist.ss <- as.matrix(as.data.frame(st_intersects(ss.hexgrid)))

# convert the edgelist into a matrix
gr.ss <- graph_from_edgelist(edgelist.ss)

# set the weights to the cellsize = distance between centroids (could also incorporate habitat resistance). 
E(gr.ss)$weight <- 173



###############################################################################################
## OBJ. 1: ESTIMATE INDIVIDUAL MAX MOVEMENT DISTANCES BETWEEN DIURNAL AND NOCTURNAL LOCATIONS
###############################################################################################

# read the diurnal/nocturnal locations
dinoc <- read.csv("DiNocLocations_ToPair.csv")

#### PAIR NOCTURNAL AND DIURNAL LOCATIONS FROM THE SAME INDIVIDUAL AND DUTY CYCLE AND IDENTIFY POINTS AT GREATEST DISTANCE
####   THESE WILL BE USED AS CENTROIDS FOR RANDOM POINT GENERATION
# exclude dawn/dusk observations (solar elevation between -6:0)
length(unique(dinoc$tag_local_identifier)) # 34 individuals = correct.  

indlist <- unique(dinoc$tag_local_identifier)

dist.trav <- data.frame()
start.time <- Sys.time()
for (ind in indlist[15:length(indlist)]){ # loop through individuals
  temp.d <- dinoc[which(dinoc$tag_local_identifier==ind & dinoc$SolarElevationCorrected_ForATMRefrctionDeg > 0),]
  temp.n <- dinoc[which(dinoc$tag_local_identifier==ind & dinoc$SolarElevationCorrected_ForATMRefrctionDeg < -6),]
  for (cycno in unique(temp.d$CycNum[which(temp.d$CycNum %in% temp.n$CycNum)])){ # loop through cycle numbers
    temp.d.cn <- temp.d[which(temp.d$CycNum==cycno),]
    temp.n.cn <- temp.n[which(temp.n$CycNum==cycno),]
    cyc.dist <- data.frame()
    for (i in 1:nrow(temp.d.cn)){ #loop through diurnal point locations, measure distance to all nocturnal point locations
      # combine single diurnal point location with all nocturnal from same cycle
      temp.dist <- rbind(temp.d.cn[i,], temp.n.cn)
      # convert to sfc object in UTM 10N projection
      coords <- st_as_sf(temp.dist, coords = c("BestLong", "BestLat"),
                         crs = 4326, stringsAsFactors = FALSE) %>%
        st_transform(crs = 32610)
      pt.int <- as.numeric(st_intersects(coords, ss.hexgrid))
      if (!(is.na(pt.int[1]))){ # if the origin point's hexgrid is NA, skip
        temp.dist$PtInt <- pt.int
        # remove NAs in endpoints, if they occur
        if (sum(is.na(pt.int))>0){
          pt.int <- pt.int[!(is.na(pt.int))] 
          temp.dist <- temp.dist[which(!(is.na(temp.dist$PtInt))),]
        } 
        # remove duplicate distances if multiple end points fall in the same hexagon
        pt.int <- c(pt.int[1], pt.int[!(duplicated(pt.int))][-1])
        temp.dist.nd <- rbind(temp.dist[1,], temp.dist[!(duplicated(temp.dist$PtInt)),][-1,])
        temp.dist.nd$DistTrav  <- c(0, distances(gr.ss, pt.int[1], pt.int[-1]))
        temp.dist2 <- rbind(temp.dist.nd[1,], temp.dist.nd[which(temp.dist.nd$DistTrav==max(temp.dist.nd$DistTrav)),][1,])
        # save the starting point and point at greatest distance
        cyc.dist <- rbind(cyc.dist, data.frame(tag_local_identifier=temp.dist2$tag_local_identifier, latitude=temp.dist2$BestLat, 
                                             longitude=temp.dist2$BestLong, CycNum=temp.dist2$CycNum, DistTrav=temp.dist2$DistTrav, Region=temp.dist2$Region, argos_lc=temp.dist2$argos_lc, Pair=i))
      } # end NA test on origin pt.int
    }
    # get the pair number with the greatest distance traveled (take the first in case of ties)
    pairno <- cyc.dist$Pair[which(cyc.dist$DistTrav==max(cyc.dist$DistTrav))][1]
    # save the start and end points of the longest distance traveled for each individual in each duty cycle
    dist.trav <- rbind(dist.trav, cyc.dist[which(cyc.dist$Pair==pairno),])
  }
}
Sys.time() - start.time
# save the file with movement centroids for each individual in each duty cycle
dist.trav$MvmtUniq <- paste(dist.trav$tag_local_identifier, dist.trav$CycNum, sep=".")
write.csv(dist.trav, file="MaxDistanceTraveled_Individual_DutyCycle.csv")

missind <- indlist[which(!(indlist %in% unique(dist.trav$tag_local_identifier)))]
# individuals without D and N detections on the same cycle: 43906 43910 43912 43892 43909 43914

########################################################

### Create 100 random points within error buffers for each observed distance traveled endpoint
### load random points back in, merge with original points, then loop through 1000 times and measure distances traveled

rand.dist <- read.csv(file="MaxDist_RandomPts.csv")
obs.dist <- read.csv(file="MaxDist_ObservedPts.csv")

obs.dist$MvmtUniq <- paste(obs.dist$tag_local_,obs.dist$CycNum,obs.dist$Pair, sep="-")
obs.dist$StartEnd <- 0
obs.dist$StartEnd[which(obs.dist$DistTrav>0)] <- 1

# merge files to get remaining fields, remove extra records from rand.dist
rand.dist <- merge(rand.dist, obs.dist[,c("Field1", "CycNum", "Region","Pair","MvmtUniq", "DistTrav","StartEnd")], by="Field1")

pairlist <- unique(rand.dist$MvmtUniq)
pairlist <- pairlist[18:length(pairlist)]
pair.df <- data.frame(MvmtUniq=pairlist, Cluster=c(rep(1,23), rep(2,23),rep(3,23),rep(4,23),rep(5,23),rep(6,23),rep(7,23),rep(8,23)))


CalcRandDist <- function(x, pair.df, rand.dist, ss.hexgrid, gr.ss,...){

  pairlist <- pair.df$MvmtUniq[which(pair.df$Cluster==x)]
  dist.trav.r <- data.frame()
  for (pr in pairlist){ # loop through individuals
    temp.0.all <- rand.dist[which(rand.dist$MvmtUniq==pr & rand.dist$StartEnd==0),]
    temp.1.all <- rand.dist[which(rand.dist$MvmtUniq==pr & rand.dist$StartEnd==1),]
    for (i in 1:3){ # sample a random start and end point 500 times
      temp.pr <- rbind(temp.0.all[sample(c(1:100),1),], temp.1.all[sample(c(1:100),1),])
      coords <- st_as_sf(temp.pr, coords = c("PS_XCoord", "PS_YCoord"),
                crs = 32610, stringsAsFactors = FALSE)
      pt.int <- as.numeric(st_intersects(coords, ss.hexgrid))
      disttrav  <- distances(gr.ss, pt.int[1], pt.int[-1])
      dist.trav.r <- rbind(dist.trav.r, data.frame(MvmtUniq=pr, DistTrav=disttrav))
    }
  }
  # save the file with movement centroids for each individual in each duty cycle
  write.csv(dist.trav.r, file=paste0("DistanceTraveled_RandomPts_C",x,".csv"))
  
}


#########################################
## RUN FUNCTION TO ESTIMATE DISTANCE TRAVELED BETWEEN ENDPOINTS (USING RANDOM POINTS WITHIN ERROR BUFFERS)
starttime <- Sys.time()
sfStop()

# Start parallel processing

# set number of cores here. I'm starting with 4, to avoid memory overloads
numCPUs = 8
sfInit(parallel=TRUE, cpus=numCPUs) # leave as is

# load the packages you need (not their dependencies) into the parallel environment using sfLibrary()
sfLibrary(knitr)
sfLibrary(sf)
sfLibrary(sp)
sfLibrary(igraph)
sfLibrary(raster)
sfLibrary(rgeos)

sfExport(list=c('pair.df', 'rand.dist', 'ss.hexgrid', 'gr.ss')) # list here any lists or data files 
#   you're iterating through or passing to the function on the fly

sfClusterApplyLB(
  1:8, # this command tells snowfall which list to iterate through and how
  fun=CalcRandDist , # name of your function
  pair.df=pair.df,
  rand.dist=rand.dist,
  ss.hexgrid=ss.hexgrid,
  gr.ss=gr.ss
)


sfStop()

#### load distance estimates, calculate means per individual
dist.trav <- read.csv(file="DistanceTraveled_RandomPts_1_17.csv")
for (i in c(1:8)){
  temp <- read.csv(file=paste0("DistanceTraveled_RandomPts_C",i,".csv"))
  dist.trav <- rbind(dist.trav, temp)
}
dist.trav$X <- NULL

# merge in Region and individual id from obs.dist
obs.dist <- read.csv(file="MaxDist_ObservedPts.csv")
obs.dist$MvmtUniq <- paste(obs.dist$tag_local_,obs.dist$CycNum,obs.dist$Pair, sep="-")
dist.trav <- merge(dist.trav, obs.dist[which(obs.dist$DistTrav>0),c("MvmtUniq","tag_local_","Region")], by="MvmtUniq", all.x=T)

write.csv(dist.trav, "DistanceTraveled_RandomPts_all.csv")

### calculate mean distance traveled per movement
library(Rmisc)
dist.trav.summ <- summarySE(dist.trav, measurevar="DistTrav", groupvars="MvmtUniq")
summary(dist.trav.summ$DistTrav)

# calculate overall mean and 95% CI
mean.overall <- mean(dist.trav.summ$DistTrav)
dist.trav.summ$var <- dist.trav.summ$sd^2
sd.overall <- sqrt(sum(dist.trav.summ$var)/nrow(dist.trav.summ))
se.overall <- sd.overall/sqrt(dist.trav.summ$N[1])
l95.overall <- mean.overall - 1.96*se.overall
u95.overall <- mean.overall + 1.96*se.overall

#######################
## BUILD MODEL TO TEST FOR REGIONAL DIFFERENCES IN DISTANCE TRAVELED

library(MuMIn)
library(AICcmodavg)
library(ggplot2)
library(nlme)

dist.trav <- read.csv("DistanceTraveled_RandomPts_all.csv")

# determine best fit distribution
library(fitdistrplus)
plotdist(dist.trav$DistTrav, histo=T, demp=T)
plotdist(log(dist.trav$DistTrav+1), histo=T, demp=T)
descdist(dist.trav$DistTrav, discrete=FALSE, boot=500)
fit_ln <- fitdist((dist.trav$DistTrav+1), "lnorm")
fit_n <- fitdist((dist.trav$DistTrav+1), "norm")
qqcomp(list(fit_ln, fit_n))

# to log or not to log: LOG
mod <- lme(DistTrav ~ Region, random=~1|tag_local_identifier/MvmtUniq, weights=varIdent(form=~1|Region), data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
plot(mod)
qqnorm(residuals(mod, type="pearson"))

dist.trav$LnDistTrav <- log(dist.trav$DistTrav + 1)
mod.ln <- lme(LnDistTrav ~ Region, random=~1|tag_local_identifier/MvmtUniq, weights=varIdent(form=~1|Region), data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
plot(mod.ln)
qqnorm(residuals(mod.ln, type="pearson")) # better

# to model heteroscedascity or not
mod.ln.nw <- lme(LnDistTrav ~ Region, random=~1|tag_local_identifier/MvmtUniq, data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
AIC(mod.ln, mod.ln.nw) # keep weights

mod.ln.gls.nw <- gls(LnDistTrav ~ Region, data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
mod.ln.gls.nw.0 <- gls(LnDistTrav ~ 1, data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
mod.ln.nw <- lme(LnDistTrav ~ Region, random=~1|tag_local_identifier/MvmtUniq, data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
mod.ln.nw.0 <- lme(LnDistTrav ~ 1, random=~1|tag_local_identifier/MvmtUniq, data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
mod.ln.0 <- lme(LnDistTrav ~ 1, random=~1|tag_local_identifier/MvmtUniq, weights=varIdent(form=~1|Region), data=dist.trav, method="ML", control = lmeControl(opt = "optim"))
mod.ln.gls <- gls(LnDistTrav ~ Region, data=dist.trav, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod.ln.gls.0 <- gls(LnDistTrav ~ 1, data=dist.trav, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))

modtab <- data.frame(modname=c("mod.ln","mod.ln.gls","mod.ln.gls.0","mod.ln.nw","mod.ln.nw.0","mod.ln.0","mod.ln.gls.nw","mod.ln.gls.nw.0"))
modtab$AICc <- c(AICc(mod.ln), AICc(mod.ln.gls), AICc(mod.ln.gls.0), AICc(mod.ln.nw), AICc(mod.ln.nw.0), AICc(mod.ln.0), AICc(mod.ln.gls.nw), AICc(mod.ln.gls.nw.0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(mod.ln), AICc(mod.ln.gls), AICc(mod.ln.gls.0), AICc(mod.ln.nw), AICc(mod.ln.nw.0), AICc(mod.ln.0), AICc(mod.ln.gls.nw), AICc(mod.ln.gls.nw.0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]

write.csv(modtab, file="Modtab_ER_RegionalDiffMovementDists.csv")

# final model - fit with REML
library(lmtest)
mod.ln <- lme(LnDistTrav ~ Region, random=~1|tag_local_identifier/MvmtUniq, weights=varIdent(form=~1|Region), data=dist.trav, method="REML", control = lmeControl(opt = "optim"))
summary(mod.ln)
mod.int <- lme(LnDistTrav ~ 1, random=~1|tag_local_identifier/MvmtUniq, weights=varIdent(form=~1|Region), data=dist.trav, method="REML", control = lmeControl(opt = "optim"))
lrtest(mod.ln, mod.int)
# Regional differences are significant
#Df LogLik Df  Chisq Pr(>Chisq)   
#1  10 -46265                        
#2   7 -46273 -3 14.701   0.002091 **

# conduct post-hoc tests
# run post-hoc tests
library(emmeans)
mod.emm <- emmeans(mod.ln, ~Region, type="response")
mod.emm.contrasts <- contrast(mod.emm, "pairwise", simple = "each", combine = TRUE, adjust = "mvt")
mod.emm.contrasts
# distances traveled longer in South PS than Central PS
#contrast                                estimate    SE  df t.ratio p.value
#Central Puget Sound - North Puget Sound  -0.2343 0.152 170 -1.541  0.3944 
#Central Puget Sound - South Puget Sound  -0.6564 0.204 170 -3.223  0.0073 *
#Central Puget Sound - Strait of Georgia  -0.5612 0.366 170 -1.532  0.3998 
#North Puget Sound - South Puget Sound    -0.4221 0.204 170 -2.068  0.1544 
#North Puget Sound - Strait of Georgia    -0.3269 0.368 170 -0.888  0.7987 
#South Puget Sound - Strait of Georgia     0.0952 0.396 170  0.241  0.9947 

# plot mean distance traveled by region
library(AICcmodavg)
dist.trav.un <- dist.trav[!(duplicated(dist.trav$MvmtUniq)),c("MvmtUniq","tag_local_identifier","Region")]

library(Rmisc)
# use predict with random effects to calculate SEs w/ individual variation
dist.trav.un$Pred <- unique(exp(predict(mod.ln, type="response"))-1)
pred.summ <- summarySE(dist.trav.un, measurevar="Pred", groupvars="Region")
pred.summ
pred.summ$RegOrder <- c(3,2,4,1)
pred.summ <- pred.summ[order(pred.summ$RegOrder),]
pred.summ$RegOrder <- as.factor(pred.summ$RegOrder)

library(RColorBrewer)
colorlist <- brewer.pal(n=4, name="Dark2")
colorlist

### PLOT MEAN REGIONAL DISTANCE TRAVELED
# plot distances traveled
png(file="SurfScoter_DistanceTraveled_LongNames.png", width=5, height=5.3, units="in" , res=600)
p <- ggplot(pred.summ, aes(x=RegOrder, y=Pred, fill=RegOrder))+  
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Pred-ci, ymax=Pred+ci), width=0.2)+
  scale_x_discrete(breaks=c("1","2","3","4"), labels=c("Northern &\nCentral\nStrait of\nGeorgia","Southern\nStrait of\nGeorgia &\nNorthern\nPuget Sound", "Central\nPuget\nSound", "Southern\nPuget\nSound"))+
  scale_fill_manual(values=colorlist)+
  labs(x="", y="Distance Traveled (m, with 95% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 14), 
        axis.title.x=element_text(size = 14))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
print(p)
dev.off()



#############################################################################
## OBJ. 2 ANALYZE NOCTURNAL HABITAT RELATIONSHIPS
##
## Use GAMMs with individual as random effect. Repeat 1000 times, each time randomly selecting
##   one use point within each buffer
################################################################################

########
# estimate mean movement distance between nocturnal use locations 
#    movement distances used to determine individual home ranges, from which availability sample is taken

# read in observed use points, measure Euclidean distances traveled between locations 
nocuse <- read.csv("NocUseLocsHabitat.csv")
length(unique(nocuse$tag_local_identifier)) # 34 individuals

# remove 2 use points in western Strait of Juan de Fuca (outside study area)
temp <- nocuse[which(nocuse$tag_local_identifier==53988 & nocuse$location_long< -123),]
nocuse <- nocuse[which(!(nocuse$OBJECTID %in% temp$OBJECTID)),]

# load the hex grid (created with 100m sides in ArcMap)
ss.hexgrid <- st_read("SalishSeaHexGrid_100mClip.shp")
edgelist.ss <- as.matrix(as.data.frame(st_intersects(ss.hexgrid)))
gr.ss <- graph_from_edgelist(edgelist.ss)
E(gr.ss)$weight <- 173

indlist <- unique(nocuse$tag_local_identifier)

dist.trav <- data.frame()
start.time <- Sys.time()
for (ind in indlist){ # loop through individuals
  temp <- nocuse[which(nocuse$tag_local_identifier==ind),]
  coords <- st_as_sf(temp, coords = c("location_long", "location_lat"),
                     crs = 4326, stringsAsFactors = FALSE) %>%
    st_transform(crs = 32610)
  disttrav <- vector()
  for (i in 1:(nrow(temp)-1)){ # loop through observations, measure movement distances between successive observations
    pt.int <- as.numeric(st_intersects(coords[c(i, i+1),], ss.hexgrid))
    if (!(is.na(pt.int[1])) & !(is.na(pt.int[2]))){ # if the origin point's hexgrid is NA, skip
      disttrav<- c(disttrav, distances(gr.ss, pt.int[1], pt.int[-1]))
    } # end NA test on origin pt.int
  } # ind loop through individual's observations
  dist.trav <- rbind(dist.trav, data.frame(tag_local_identifier=ind, MeanMvmt=mean(disttrav, na.rm=TRUE)))
}
Sys.time() - start.time
write.csv(dist.trav, file="MeanMovementDistance_NocturnalUseLocations.csv")

# add movement distance to nocuse dataset
nocuse <- merge(nocuse, dist.trav, by="tag_local_identifier", all.x=T)

# get number of observations per individual
library(plyr)
nocuse$dum <- 1
nocuse.numind <- ddply(nocuse, c("tag_local_identifier"), summarise, NumObs=sum(dum))
write.csv(nocuse.numind, file="NocUse_NumberObservationsByIndividual.csv")
nocuse$dum <- NULL

#######################
## CREATE AND MERGE MCPs and BUFFERS AROUND POINTS BASED ON INDIVIDUAL'S MEAN MOVEMENT DISTANCE
library(adehabitatHR)

WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coords.sp <- nocuse[,c("tag_local_identifier", "location_long", "location_lat")]
coordinates(coords.sp) <- c("location_long", "location_lat")
proj4string(coords.sp) <- WGS84
coords.sp <- spTransform(coords.sp, crs(coords))
coords.mcp <- mcp(coords.sp, percent=100)
shapefile(coords.mcp, file="MinimumConvexPolygons_AllNocUse.shp")

coords.mcp.sf <- st_read("MinimumConvexPolygons_AllNocUse.shp")
for (ind in indlist){
  tempind <- nocuse[which(nocuse$tag_local_identifier==ind),]
  coords <- st_as_sf(tempind, coords = c("location_long", "location_lat"),
                     crs = 4326, stringsAsFactors = FALSE) %>%
    st_transform(crs = 32610)
  coords_buffer <- st_buffer(coords, mean(tempind$MeanMvmt, na.rm=T))
  #plot(st_geometry(coords_buffer))
  ind.merge <- st_union(coords_buffer, coords.mcp.sf[which(coords.mcp.sf$id==ind),])
  ind.diss <- st_union(ind.merge)
  #ind.diss.id <- st_sf(data.frame(tag_local_identifier=ind, geom=ind.diss))
  if (ind==indlist[1]){
    MCP_Buffers_Merged <- ind.diss
  } else {
    MCP_Buffers_Merged <- c(MCP_Buffers_Merged, ind.diss)
  }
}
MCP_Buffers_Merged.id <- st_sf(data.frame(tag_local_identifier=indlist, geom=MCP_Buffers_Merged))
st_write(MCP_Buffers_Merged.id, paste0("MCP_Buffers_Merged_NocUse.shp"))

# read in use and non use data with habitat variables
avail1 <- read.csv("Noc_Availability_VesselDensity_DistToCoast.csv")
avail2 <- read.csv("Noc_Availability_Depth_TidalCurrent.csv")
avail3 <- read.csv("Noc_Availability_Regions.csv")
randpts1 <- read.csv("RandomUsePts_Depth_ShorelineDist.csv")
randpts2 <- read.csv("RandomUsePts_Current_VesselDensity.csv")
randpts3 <- read.csv("RandomUsePts_Region.csv")

# for randpts layers: OBJECTID holds IDs for unique observations.
#   For each observation, an error buffer was created and 100 random points were generated within each buffer
# merge random use points data
randpts <- merge(randpts1, randpts2[,c("FID","TidalCurrent","VesselDensity")], by="FID")
randpts <- randpts[,c("OBJECTID", "PS_XCoord", "PS_YCoord", "Depth", "DistanceToShore", "TidalCurrent", "VesselDensity")]
# merge in region. remove 2 points in western Strait of Juan de Fuca
randpts <- merge(randpts, randpts3[!(duplicated(randpts3$OBJECTID)),c("OBJECTID", "Region")], by="OBJECTID", all.x=F)

# add in individual ID
nocobs <- read.csv(file="NocUseLocsHabitat.csv")
randpts <- merge(randpts, nocobs[,c("tag_local_identifier","OBJECTID")], by="OBJECTID", all.x=T)
randpts$UseDetermination <- 1
# replace Depth=-9999 with 0
randpts$Depth[which(randpts$Depth< -9000)] <- 0
write.csv(randpts, file="RandomNocturnalUsePoints_Merged.csv")

# merge nonuse pt datasets

# fix extra point in regions dataset (borders 2 regions)
avail3[duplicated(avail3$FID_NocUse),]
avail3[which(avail3$FID_NocUse==7390),]
avail3$Region[which(avail3$FID_NocUse==7390)] <- "Central Puget Sound"
write.csv(avail3, file="Noc_Availability_Regions.csv")

# replace all OBJECTID in nonuse with '0'
availpts <- merge(avail1, avail2[,c("FID","Depth","TidalCurrent")], by="FID")
availpts <- merge(availpts, avail3[,c("FID","Region")], by="FID")
availpts <- availpts[,c("PS_XCoord", "PS_YCoord", "tag_local_identifier", "Depth", "DistanceToShore", "TidalCurrent", "VesselDensity","Region")]
availpts$OBJECTID <- 0
availpts$UseDetermination <- 0

# duplicate 4 points to bring total to 10000
availpts <- rbind(availpts, availpts[sample(c(1:9996),4,replace=F),])
# replace Depth=-9999 with NA
availpts$Depth[which(availpts$Depth< -9000)] <- 0
write.csv(availpts, file="NocturnalAvailabilityPoints_Merged.csv")


###############################
## BEGIN RUNNING GAMMs
##############################

library(mgcViz)
library(ape)
library(pROC)
library(mgcv)
library(Hmisc)
library(raster)

# read in observed use, random use, and availability data
nocuse <- read.csv("NocUseLocsHabitat.csv")
temp <- nocuse[which(nocuse$tag_local_identifier==53988 & nocuse$location_long< -123),]
nocuse <- nocuse[which(!(nocuse$OBJECTID %in% temp$OBJECTID)),]
randpts <- read.csv(file="RandomNocturnalUsePoints_Merged.csv")
randpts$X <- NULL
availpts <- read.csv(file="NocturnalAvailabilityPoints_Merged.csv")
availpts$X <- NULL
availpts$X.1 <- NULL

# update regions in nocuse
nocuse$Region <- NULL
nocuse <- merge(nocuse, randpts[!(duplicated(randpts$OBJECTID)),c("OBJECTID","Region")], by="OBJECTID", all.x=T)
nocuse <- nocuse[,c("OBJECTID", "Use_Determination", "tag_local_identifier", "location_long", "location_lat", "RMS_Tidal_", "DEPTH", "RASTERVALU", "NEAR_DIST", "Region")]
names(nocuse) <- c("OBJECTID", "UseDetermination", "tag_local_identifier", "PS_XCoord", "PS_YCoord", "TidalCurrent", "Depth", "VesselDensity", "DistanceToShore", "Region")

# fix Depth in nocuse
nocuse$Depth[which(nocuse$Depth < -9000)] <- 0

rsfdat <- rbind(nocuse, availpts)

rsfdat$Region <- as.factor(rsfdat$Region)
library(AICcmodavg)
library(MuMIn)

# build the GAMM - model selection
#m <- gam(UseDetermination ~ s(TidalCurrent, k=-1), 
#                 data=rsfdat, family=binomial, method="GCV.Cp")
#m <- gam(UseDetermination ~ s(TidalCurrent, k=500), 
#         data=rsfdat, family=binomial, method="ML")
#m <- gam(UseDetermination ~ s(Depth, k=-1) + s(VesselDensity, k=-1) + s(DistanceToShore, k=-1) 
#         +s(TidalCurrent, k=-1) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
#         data=rsfdat, family=binomial, method="GCV.Cp")
#m <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
#         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
#         data=rsfdat, family=binomial, method="ML")
# check model fit
#gam.check(m)
#plot.gam(m)
#summary(m)

# compare models with various combinations of random effects
m2 <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
          +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
          data=rsfdat, family=binomial, method="ML")
m1.r <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
          +s(TidalCurrent, k=10) + s(Region, bs="re"), 
          data=rsfdat, family=binomial, method="ML")
m1.i <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
            +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re"), 
            data=rsfdat, family=binomial, method="ML")
m0 <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
          +s(TidalCurrent, k=10) , 
          data=rsfdat, family=binomial, method="ML")
modtab <- data.frame(modname=c("m2","m1.r", "m1.i", "m0"))
modtab$AICc <- c(AICc(m2), AICc(m1.r), AICc(m1.i), AICc(m0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(m2), AICc(m1.r), AICc(m1.i), AICc(m0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]

write.csv(modtab, file="Modtab_ER_NoctHabSeln_RandomEffects.csv")
summary(m2)[10]
summary(m1.r)[10]
summary(m1.i)[10]
summary(m0)[10]

# compare models with various fixed effect smooths vs linear predictors
m <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
          +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
          data=rsfdat, family=binomial, method="ML")
m.d <- gam(UseDetermination ~ Depth + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
         data=rsfdat, family=binomial, method="ML")
m.v <- gam(UseDetermination ~ s(Depth, k=30) + VesselDensity + s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
         data=rsfdat, family=binomial, method="ML")
m.s <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + DistanceToShore 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
         data=rsfdat, family=binomial, method="ML")
m.c <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
         + TidalCurrent + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
         data=rsfdat, family=binomial, method="ML")
modtab <- data.frame(modname=c("m","m.d", "m.v", "m.s", "m.c"))
modtab$AICc <- c(AICc(m), AICc(m.d), AICc(m.v), AICc(m.s), AICc(m.c))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(m), AICc(m.d), AICc(m.v), AICc(m.s), AICc(m.c)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_NoctHabSeln_FixedSmooths.csv")
summary(m.v)[10]
summary(m.c)[10]
summary(m.s)[10]
summary(m.d)[10]

## full subsets fixed effects selection
m <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.vsc <- gam(UseDetermination ~ s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.dsc <- gam(UseDetermination ~ s(Depth, k=30) + s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.dvc <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.dvs <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
        + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.dv <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.ds <- gam(UseDetermination ~ s(Depth, k=30) + s(DistanceToShore, k=10) 
         + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.dc <- gam(UseDetermination ~ s(Depth, k=30) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.vs <- gam(UseDetermination ~ s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
          + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.vc <- gam(UseDetermination ~  s(VesselDensity, k=10)  
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.sc <- gam(UseDetermination ~  s(DistanceToShore, k=10) 
         +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.d <- gam(UseDetermination ~ s(Depth, k=30)  + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.v <- gam(UseDetermination ~ s(VesselDensity, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.s <- gam(UseDetermination ~ s(DistanceToShore, k=10)  + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m.c <- gam(UseDetermination ~ s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
m0 <- gam(UseDetermination ~ s(tag_local_identifier, bs="re")  + s(Region, bs="re"), data=rsfdat, family=binomial, method="ML")
modtab <- data.frame(modname=c("m", "m.vsc", "m.dsc", "m.dvc", "m.dvs", "m.dv", "m.ds", "m.dc", "m.vs", "m.vc", "m.sc", "m.d", "m.v", "m.s", "m.c", "m0"))
modtab$AICc <- c(AICc(m), AICc(m.vsc), AICc(m.dsc), AICc(m.dvc), AICc(m.dvs), AICc(m.dv), AICc(m.ds), AICc(m.dc), AICc(m.vs), AICc(m.vc), AICc(m.sc), AICc(m.d), AICc(m.v), AICc(m.s), AICc(m.c), AICc(m0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(m), AICc(m.vsc), AICc(m.dsc), AICc(m.dvc), AICc(m.dvs), AICc(m.dv), AICc(m.ds), AICc(m.dc), AICc(m.vs), AICc(m.vc), AICc(m.sc), AICc(m.d), AICc(m.v), AICc(m.s), AICc(m.c), AICc(m0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_NoctHabSeln_FixedEffects.csv")

summary(m.vsc)[10]
summary(m.dsc)[10]
summary(m.dvc)[10]
summary(m.dvs)[10]
summary(m.dv)[10]
summary(m.ds)[10]
summary(m.dc)[10]
summary(m.vs)[10]
summary(m.vc)[10]
summary(m.sc)[10]
summary(m.d)[10]
summary(m.v)[10]
summary(m.s)[10]
summary(m.c)[10]
summary(m0)[10]

## check for collinearity among predictors
cormat <- cor(rsfdat[,c("Depth", "DistanceToShore","TidalCurrent", "VesselDensity")], use="complete.obs")
max(cormat[which(cormat<1)]) # 0.0452155
min(cormat[which(!(is.na(cormat)))]) # -0.06465363
write.csv(cormat, file="SurfScoter_CorMatrix.csv")

# check for multicollinearity
concurvity(m)
# values: 0-1, 0 = no problem, 1 = total lack of identifiability
#          para   s(Depth) s(VesselDensity) s(DistanceToShore) s(TidalCurrent) s(tag_local_identifier) s(Region)
#worst       1 0.16358545       0.11352952         0.15547885       0.2074182               0.9898745     0.1326337
#observed    1 0.04160140       0.01943402         0.02269246       0.1673254               0.9898745     0.1326337
#estimate    1 0.00642889       0.02446919         0.05524210       0.1042997               0.9898745     0.4434816

# loop through 500 model iterations, select one random use point for each observation and comparing with paired set of availability points
#gam.modfit <- data.frame()
start.time <- Sys.time()
for (i in 75:500){
  tempuse <- data.frame()
  for (obs in unique(randpts$OBJECTID)){
    tempobs <- randpts[which(randpts$OBJECTID==obs),]
    tempuse <- rbind(tempuse, tempobs[sample(c(1:100),1),])
  }
  rsfdat <- rbind(tempuse, availpts)
  # create 5 latitudinal folds & 2 longitudinal folds for spatial cross-validation  
  rsfdat$LatFold <- as.numeric(cut2(rsfdat$PS_YCoord, g=5))
  rsfdat$LongFold <- as.numeric(cut2(rsfdat$PS_XCoord, g=2))
  rsfdat$GridFold <- paste(rsfdat$LatFold, rsfdat$LongFold, sep=".")
  # build the GAMM
  m <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
           +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
          data=rsfdat, family=binomial, method="REML")
  save(m, file=paste0("gamms/gam_",i,".rda"))
  ### calculate Moran's i to evalute spatial correlation in residuals
  gam.resid<-residuals(m,type="pearson")
  gamresxy<-data.frame(gam.resid,rsfdat$PS_XCoord,rsfdat$PS_YCoord)
  colnames(gamresxy) <- c("Resid","Longitude","Latitude")
  # generate distance matrix, invert it, and replace diagonals with 0
  gamres.dists <- as.matrix(dist(cbind(gamresxy$Longitude, gamresxy$Latitude)))
  gamres.dists.inv <- 1/gamres.dists
  diag(gamres.dists.inv) <- 0
  gamres.dists.inv[is.infinite(gamres.dists.inv)] <- 0
  # calculate Moran's I of residuals. I = observed. I > 0.3 is evidence of problematic spatial autocorrelation (Lichstein et al., 2002)
  Moran <- Moran.I(gamresxy$Resid, gamres.dists.inv)$observed 
  
  # calculate self and cross-validated AUC
  #Pred <- plogis(predict.gam(m,rsfdat))
  auc.self <- auc(roc(rsfdat$UseDetermination ~ plogis(predict.gam(m,rsfdat))))
  auc.cv <- vector()
  for (gf in unique(rsfdat$GridFold)){
    preddat <- rsfdat[which(rsfdat$GridFold==gf),]
    cvdat <- rsfdat[which(rsfdat$GridFold!=gf),]
    m.cv <- gam(UseDetermination ~ s(Depth, k=30) + s(VesselDensity, k=10) + s(DistanceToShore, k=10) 
             +s(TidalCurrent, k=10) + s(tag_local_identifier, bs="re")  + s(Region, bs="re"), 
             data=cvdat, family=binomial, method="REML")
    if (sum(preddat$UseDetermination)>0){
      auc.cv <- c(auc.cv, auc(roc(preddat$UseDetermination ~ plogis(predict.gam(m,preddat)))))
    }
  }
  auc.cv.mean <- mean(auc.cv)
  auc.cv.se <- sd(auc.cv, na.rm=T)/sqrt(length(auc.cv))
  
  # use full model to predict from raster stack
  #rpred <- predict(rstack, m)
  #rpred.df <- as.data.frame(rpred)
  #write.csv(rpred.df, file=paste0("gamms/gamPreds_",i,".csv"))
  
  # save modelfit  stats
  gam.modfit <- rbind(gam.modfit, data.frame(Rep=i, AUC.self=auc.self, AUC.cv.mean=auc.cv.mean, AUC.cv.se=auc.cv.se, MoranI=Moran))

  print(paste0("iteration ",i))
}
Sys.time()-start.time
write.csv(gam.modfit, file="GAM_ModelFitStats.csv")


#########################
### predict from GAM

# load covariates
Depth <- raster("Depth_SalishSea.tif")
DistanceToShore <- raster("DistanceToShore_SalishSea.tif")
VesselDensity <- raster("VesselDensity_SalishSea2.tif")
TidalCurrent <- raster("TidalCurrent_SalishSea.tif")
Region <- raster("SalishSea_Regions.tif")
# fix rasters with incorrect extents or number of rows
dep.pts <- as(Depth, 'SpatialPointsDataFrame')
temp <- Region
temp <- resample(Region, Depth, "bilinear")
crs(temp) <- crs(Depth)
assign("Region", temp)

tag_local_identifier <- raster("SalishSea_zeroes.tif")
rstack <- stack(Depth, DistanceToShore, VesselDensity, TidalCurrent, Region, tag_local_identifier)
summary(rstack)

names(rstack) <- c("Depth", "DistanceToShore", "VesselDensity", "TidalCurrent", "Region", "tag_local_identifier")

# convert to df
dep.pts <- as(Depth, 'SpatialPointsDataFrame')
dts.pts <- as(DistanceToShore, 'SpatialPointsDataFrame')
vd.pts <- as(VesselDensity, 'SpatialPointsDataFrame')
tc.pts <- as(TidalCurrent, 'SpatialPointsDataFrame')
reg.pts <- as(Region, 'SpatialPointsDataFrame')
preddat <- as.data.frame(dep.pts)
preddat$xy <- paste(preddat$x, preddat$y, sep=".")
dts.pts$xy <- paste(dts.pts$x, dts.pts$y, sep=".")
dts.pts <- dts.pts[,c("xy","DistanceToShore_SalishSea")]
vd.pts$xy <- paste(vd.pts$x, vd.pts$y, sep=".")
vd.pts <- vd.pts[,c("xy","VesselDensity_SalishSea2")]
tc.pts$xy <- paste(tc.pts$x, tc.pts$y, sep=".")
tc.pts <- tc.pts[,c("xy","TidalCurrent_SalishSea")]
reg.pts$xy <- paste(reg.pts$x, reg.pts$y, sep=".")
reg.pts <- reg.pts[,c("xy","SalishSea_Regions")]
preddat <- merge(preddat, dts.pts, by="xy")
preddat <- merge(preddat, vd.pts, by="xy")
preddat <- merge(preddat, tc.pts, by="xy")
preddat <- merge(preddat, reg.pts, by="xy")
preddat$tag_local_identifier <- 0
names(preddat) <- c("xy", "Depth", "DistanceToShore", "VesselDensity", "TidalCurrent", "Region", "PS_XCoord", "PS_YCoord", "tag_local_identifier")

# make Region values factors
preddat$Region[which(preddat$Region<2)] <- "Strait of Georgia"
preddat$Region[which(preddat$Region>1.9 & preddat$Region<3)] <- "North Puget Sound"
preddat$Region[which(preddat$Region>2.9 & preddat$Region<4)] <- "Central Puget Sound"
preddat$Region[which(preddat$Region==4)] <- "South Puget Sound"
preddat$Region <- as.factor(preddat$Region)

# use full model to predict from raster stack
# predict from 200 randomly-selected models # could use 43906
preddat$tag_local_identifier <- 43906
modlist <- sample(c(1:500), 200, replace=F)
for (i in modlist){
  load(file=paste0("gamms/gam_",i,".rda"))
  rpred <- plogis(predict.gam(m, preddat)) # start 8:33
  prednames <- names(preddat) 
  preddat$rpred <- rpred
  names(preddat) <- c(prednames, paste0("Pred",i))
  
}
preddat$PredMean <- rowMeans(preddat[,c(10:ncol(preddat))], na.rm=T)
summary(preddat$PredMean)

library(FIACH)
preddat$PredSE <- rowsd(preddat[,c(10:ncol(preddat))])/sqrt(500)
write.csv(preddat[,c(1:9,210:ncol(preddat))], file=paste0("gamms/PredictedOccurrenceProbability.csv"))

# write predicted occurrence raster
predrast=rasterFromXYZ(preddat[,c("PS_XCoord","PS_YCoord","PredMean")]) 
crs(predrast)=crs(Depth)
res(predrast)=res(Depth)
extent(predrast) = extent(Depth)
writeRaster(predrast, paste("SurfScoterNocturnalHabitatRaster.tif",sep=""), format="GTiff", overwrite=TRUE)

###############
#### PRODUCE RESPONSE PLOTS
randpts <- read.csv(file="RandomNocturnalUsePoints_Merged.csv")
randpts$X <- NULL
randpts$X.1 <- NULL
availpts <- read.csv(file="NocturnalAvailabilityPoints_Merged.csv")
availpts$X <- NULL
availpts$X.1 <- NULL

alldat <- rbind(randpts, availpts)

# build newdat files for response curves
newdat.depth <- data.frame(Depth=seq(min(alldat$Depth), max(alldat$Depth), by=((max(alldat$Depth)-min(alldat$Depth))/100)),
                           DistanceToShore=rep(mean(alldat$DistanceToShore),101), TidalCurrent=rep(mean(alldat$TidalCurrent), 101),
                           VesselDensity=rep(mean(alldat$VesselDensity), 101), Region=rep("North Puget Sound", 101), 
                           tag_local_identifier=rep(43906, 101))
newdat.shore <- data.frame(Depth=rep(mean(alldat$Depth), 101),
                           DistanceToShore=seq(min(alldat$DistanceToShore), max(alldat$DistanceToShore), by=((max(alldat$DistanceToShore)-min(alldat$DistanceToShore))/100)), 
                           TidalCurrent=rep(mean(alldat$TidalCurrent), 101),
                           VesselDensity=rep(mean(alldat$VesselDensity), 101), Region=rep("North Puget Sound", 101), 
                           tag_local_identifier=rep(43906, 101))
newdat.current <- data.frame(Depth=rep(mean(alldat$Depth), 101), DistanceToShore=rep(mean(alldat$DistanceToShore),101), 
                             TidalCurrent=seq(min(alldat$TidalCurrent), max(alldat$TidalCurrent), by=((max(alldat$TidalCurrent)-min(alldat$TidalCurrent))/100)), 
                             VesselDensity=rep(mean(alldat$VesselDensity), 101), Region=rep("North Puget Sound", 101), 
                             tag_local_identifier=rep(43906, 101))
newdat.vessel <- data.frame(Depth=rep(mean(alldat$Depth), 101), DistanceToShore=rep(mean(alldat$DistanceToShore),101), 
                            TidalCurrent=rep(mean(alldat$TidalCurrent), 101),
                            VesselDensity=seq(min(alldat$VesselDensity), max(alldat$VesselDensity), by=((max(alldat$VesselDensity)-min(alldat$VesselDensity))/100)), 
                            Region=rep("North Puget Sound", 101), 
                            tag_local_identifier=rep(43906, 101))

newdat.depth.se <- newdat.depth
newdat.shore.se <- newdat.shore
newdat.current.se <- newdat.current
newdat.vessel.se <- newdat.vessel

### extract parameters and response curves
start.time <- Sys.time()
for (i in 1:500){
  load(file=paste0("gamms/gam_",i,".rda"))
  depth.names <- names(newdat.depth)
  depth.names.se <- names(newdat.depth.se)
  newdat.depth$Pred <- predict.gam(m, newdat.depth, type="response")
  newdat.depth.se$Pred.se <- predict.gam(m, newdat.depth.se, type="response", se.fit=T)$se.fit
  names(newdat.depth) <- c(depth.names, paste0("Pred",i))
  names(newdat.depth.se) <- c(depth.names.se, paste0("PredSE",i))
  shore.names <- names(newdat.shore)
  shore.names.se <- names(newdat.shore.se)
  newdat.shore$Pred <- predict.gam(m, newdat.shore, type="response")
  newdat.shore.se$Pred.se <- predict.gam(m, newdat.shore.se, type="response", se.fit=T)$se.fit
  names(newdat.shore) <- c(shore.names, paste0("Pred",i))
  names(newdat.shore.se) <- c(shore.names.se, paste0("PredSE",i))
  current.names <- names(newdat.current)
  current.names.se <- names(newdat.current.se)
  newdat.current$Pred <- predict.gam(m, newdat.current, type="response")
  newdat.current.se$Pred.se <- predict.gam(m, newdat.current.se, type="response", se.fit=T)$se.fit 
  names(newdat.current) <- c(current.names, paste0("Pred",i))
  names(newdat.current.se) <- c(current.names.se, paste0("PredSE",i))
  vessel.names <- names(newdat.vessel)
  vessel.names.se <- names(newdat.vessel.se)
  newdat.vessel$Pred <- predict.gam(m, newdat.vessel, type="response")
  newdat.vessel.se$Pred.se <- predict.gam(m, newdat.vessel.se, type="response", se.fit=T)$se.fit
  names(newdat.vessel) <- c(vessel.names, paste0("Pred",i))
  names(newdat.vessel.se) <- c(vessel.names.se, paste0("Pred",i))
}
paramdat <- data.frame(Iteration=seq(1,500,1), Param=rep(c("Depth","Vessel","Shore","Current","Indiv","Region"),500), EDF=edfdat, ChiSq=chisqdat, P=pdat, Rsq=rsqdat, DevExpl=devexpl)

# summarize mean and SE
library(FIACH)
colMeans(paramdat[which(paramdat$Param=="Depth"),c(3:ncol(paramdat))])
#  EDF        ChiSq            P          Rsq      DevExpl 
#1.466493e+01 2.674866e+02 1.525773e-42 2.080057e-01 2.399648e-01
colsd(paramdat[which(paramdat$Param=="Depth"),c(3:ncol(paramdat))])/sqrt(500)
#2.813072e-02 5.375416e-01 6.923242e-43 2.115873e-04 2.159944e-04

colMeans(paramdat[which(paramdat$Param=="Shore"),c(3:ncol(paramdat))])
#  EDF        ChiSq            P          Rsq      DevExpl 
#6.422300e+00 1.895738e+02 1.019224e-32 2.080057e-01 2.399648e-01
colsd(paramdat[which(paramdat$Param=="Shore"),c(3:ncol(paramdat))])/sqrt(500)
#6.158777e-03 5.165444e-01 7.416839e-33 2.115873e-04 2.159944e-04

colMeans(paramdat[which(paramdat$Param=="Current"),c(3:ncol(paramdat))])
#  EDF        ChiSq            P          Rsq      DevExpl 
#6.445994e+00 2.038186e+02 2.650782e-34 2.080057e-01 2.399648e-01
colsd(paramdat[which(paramdat$Param=="Current"),c(3:ncol(paramdat))])/sqrt(500)
#1.533025e-02 4.610530e-01 1.824724e-34 2.115873e-04 2.159944e-04

colMeans(paramdat[which(paramdat$Param=="Vessel"),c(3:ncol(paramdat))])
#  EDF        ChiSq            P          Rsq      DevExpl 
#4.308298e+00 7.633960e+01 6.298837e-13 2.080057e-01 2.399648e-01
colsd(paramdat[which(paramdat$Param=="Vessel"),c(3:ncol(paramdat))])/sqrt(500)
#3.990493e-02 2.951885e-01 1.718452e-13 2.115873e-04 2.159944e-04

write.csv(paramdat, file="GAM_Parameters.csv")

depth.resp <- data.frame(Depth=newdat.depth$Depth, Pred=rowMeans(newdat.depth[,c(7:ncol(newdat.depth))]), SE=rowMeans(newdat.depth.se[,c(7:ncol(newdat.depth.se))]))
shore.resp <- data.frame(DistanceToShore=newdat.shore$DistanceToShore, Pred=rowMeans(newdat.shore[,c(7:ncol(newdat.shore))]), SE=rowMeans(newdat.shore.se[,c(7:ncol(newdat.shore.se))]))
current.resp <- data.frame(TidalCurrent=newdat.current$TidalCurrent, Pred=rowMeans(newdat.current[,c(7:ncol(newdat.current))]), SE=rowMeans(newdat.current.se[,c(7:ncol(newdat.current.se))]))
vessel.resp <- data.frame(VesselDensity=newdat.vessel$VesselDensity, Pred=rowMeans(newdat.vessel[,c(7:ncol(newdat.vessel))]), SE=rowMeans(newdat.vessel.se[,c(7:ncol(newdat.vessel.se))]))

write.csv(depth.resp, file="ResponseData_Depth.csv")
write.csv(shore.resp, file="ResponseData_DistanceToShore.csv")
write.csv(current.resp, file="ResponseData_TidalCurrent.csv")
write.csv(vessel.resp, file="ResponseData_VesselDensity.csv")
Sys.time() - start.time


# plot responses
depth.resp <- read.csv(file="ResponseData_Depth.csv")
shore.resp <- read.csv(file="ResponseData_DistanceToShore.csv")
current.resp <- read.csv(file="ResponseData_TidalCurrent.csv")
vessel.resp <- read.csv(file="ResponseData_VesselDensity.csv")

## Depth
depth.resp$L95=depth.resp$Pred-(1.96*depth.resp$SE)
depth.resp$U95=depth.resp$Pred+(1.96*depth.resp$SE)
depth.resp$Depth[which(depth.resp$Depth>0)] <- 0
depth.resp$Depth <- abs(depth.resp$Depth)

p.depth <- ggplot(data=depth.resp) +
  geom_ribbon(aes(x=Depth, ymin=L95, ymax=U95), fill="#120be0", alpha=0.3) +
  geom_line(aes(x=Depth, y=Pred), colour="#120be0", lwd=1.5) +
  xlab("Depth (m below sea level)") + ylab("Predicted occurrence probability") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.title.x=element_text(size = 12))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p.depth)

## Distance to shore
shore.resp$L95=shore.resp$Pred-(1.96*shore.resp$SE)
shore.resp$U95=shore.resp$Pred+(1.96*shore.resp$SE)

p.shore <- ggplot(data=shore.resp) +
  geom_ribbon(aes(x=DistanceToShore, ymin=L95, ymax=U95), fill="#1a524d", alpha=0.3) +
  geom_line(aes(x=DistanceToShore, y=Pred), colour="#1a524d", lwd=1.5) +
  xlab("Distance to shore (m)") + ylab("Predicted occurrence probability") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.title.x=element_text(size = 12))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p.shore)

## Tidal current
current.resp$L95=current.resp$Pred-(1.96*current.resp$SE)
current.resp$U95=current.resp$Pred+(1.96*current.resp$SE)

p.current <- ggplot(data=current.resp) +
  geom_ribbon(aes(x=TidalCurrent, ymin=L95, ymax=U95), fill="#694c32", alpha=0.3) +
  geom_line(aes(x=TidalCurrent, y=Pred), colour="#694c32", lwd=1.5) +
  xlab("Tidal current") + ylab("Predicted occurrence probability") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.title.x=element_text(size = 12))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p.current)

## Vessel density
vessel.resp$L95=vessel.resp$Pred-(1.96*vessel.resp$SE)
vessel.resp$U95=vessel.resp$Pred+(1.96*vessel.resp$SE)

p.vessel <- ggplot(data=vessel.resp) +
  geom_ribbon(aes(x=VesselDensity, ymin=L95, ymax=U95), fill="#91110d", alpha=0.3) +
  geom_line(aes(x=VesselDensity, y=Pred), colour="#91110d", lwd=1.5) +
  xlab("Vessel density") + ylab("Predicted occurrence probability") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.title.x=element_text(size = 12))+
  theme(legend.position = "bottom", legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p.vessel)


# Arranging the plot
png(file="SurfScoter_ResponsePlotPanel.png", width=9, height=9, units="in" , res=600)
p.all <- ggarrange(p.depth, p.shore, p.current, p.vessel, 
                   ncol = 2, nrow = 2,  align = "hv", 
                   widths = c(2, 2), heights = c(2, 2),
                   labels = c("1","2","3","4"),
                   common.legend = FALSE)
print(p.all)
dev.off()


###############################################################################################
## OBJ: ASSESS DIFFERENCES IN DEPTH, DISTANCE FROM SHORE, VESSEL DENSITY, AND 
##      TIDAL CURRENTS BETWEEN DIURNAL AND NOCTURNAL LOCATIONS
###############################################################################################

# read the diurnal/nocturnal locations
dinoc <- read.csv("DiNocLocations_ToPair.csv")

# separate into diurnal and nocturnal files
# exclude dawn/dusk observations (solar elevation between -6:0)
dinoc$Diurnal <- "Dusk"
dinoc$Diurnal[which(dinoc$SolarElevationCorrected_ForATMRefrctionDeg > 0)] <- "Day"
dinoc$Diurnal[which(dinoc$SolarElevationCorrected_ForATMRefrctionDeg < -6)] <- "Night"

library(nlme)
library(MuMIn)
library(AICcmodavg)
library(ggplot2)

# Depth
summary(dinoc$Depth)
dinoc2 <- dinoc[which(dinoc$Depth>-9999 & dinoc$Diurnal != "Dusk"),]
hist(dinoc2$Depth)
dinoc2$Depth0 <- dinoc2$Depth
dinoc2$Depth0[which(dinoc2$Depth0>0)] <- 0
dinoc2$Depth0 <- dinoc2$Depth0 + 0.1
dinoc2$LnDepth <- log(abs(dinoc2$Depth0))
hist(dinoc2$LnDepth)

# fit and select models
mod.r.w <- lme(LnDepth ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.w.0 <- lme(LnDepth ~ 1, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r <- lme(LnDepth ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.0 <- lme(LnDepth ~ 1, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w <- lme(LnDepth ~ Diurnal, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w.0 <- lme(LnDepth ~ 1, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri <- lme(LnDepth ~ Diurnal, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.0 <- lme(LnDepth ~ 1, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w <- lme(LnDepth ~ Diurnal, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w.0 <- lme(LnDepth ~ 1, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr <- lme(LnDepth ~ Diurnal, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.0 <- lme(LnDepth ~ 1, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.w <- gls(LnDepth ~ Diurnal, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod.w.0 <- gls(LnDepth ~ 1, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod <- gls(LnDepth ~ Diurnal, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.0 <- gls(LnDepth ~ 1, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
modtab <- data.frame(modname=c("mod.r.w","mod.r.w.0","mod.r","mod.r.0","mod.ri.w","mod.ri.w.0","mod.ri","mod.ri.0","mod.rr.w","mod.rr.w.0","mod.rr","mod.rr.0","mod.w","mod.w.0","mod","mod.0"))
modtab$AICc <- c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_DiurnalNocturnalDiff_Depth.csv")

# summarize best model
summary(mod.r.w)
#                  Value  Std.Error   DF   t-value p-value
#(Intercept)  1.943170 0.27118257 4121  7.165542       0
#DiurnalNight 1.320852 0.04287767 4121 30.805120       0

newdata <- rbind(dinoc2[which(dinoc2$Diurnal=="Day"),][1,], dinoc2[which(dinoc2$Diurnal=="Night"),][1,])
dinocsumm <- data.frame(Var=rep("Depth",2), TimeOfDay=c("Day","Night"))
dinocsumm$Fit <- exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit)-0.1
dinocsumm$L95 <- exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit - 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit))-0.1
dinocsumm$U95 <- exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit + 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit))-0.1



#### Distance to shore
summary(dinoc$NEAR_DIST_SHORELINE)
dinoc2 <- dinoc[which(dinoc$Diurnal != "Dusk"),]
hist(dinoc2$NEAR_DIST_SHORELINE)
dinoc2$LnDist <- log(dinoc2$NEAR_DIST_SHORELINE)
hist(dinoc2$LnDist)

# fit and select models
mod.r.w <- lme(LnDist ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.w.0 <- lme(LnDist ~ 1, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r <- lme(LnDist ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.0 <- lme(LnDist ~ 1, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w <- lme(LnDist ~ Diurnal, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w.0 <- lme(LnDist ~ 1, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri <- lme(LnDist ~ Diurnal, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.0 <- lme(LnDist ~ 1, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w <- lme(LnDist ~ Diurnal, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w.0 <- lme(LnDist ~ 1, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr <- lme(LnDist ~ Diurnal, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.0 <- lme(LnDist ~ 1, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.w <- gls(LnDist ~ Diurnal, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod.w.0 <- gls(LnDist ~ 1, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod <- gls(LnDist ~ Diurnal, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.0 <- gls(LnDist ~ 1, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
modtab <- data.frame(modname=c("mod.r.w","mod.r.w.0","mod.r","mod.r.0","mod.ri.w","mod.ri.w.0","mod.ri","mod.ri.0","mod.rr.w","mod.rr.w.0","mod.rr","mod.rr.0","mod.w","mod.w.0","mod","mod.0"))
modtab$AICc <- c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_DiurnalNocturnalDiff_DistToShore.csv")

# summarize best model
summary(mod.r.w)
#                  Value Std.Error   DF  t-value p-value
#(Intercept)  6.047105 0.2761993 4163 21.89399       0
#DiurnalNight 0.971750 0.0286077 4163 33.96812       0

dinocsumm <- rbind(dinocsumm, data.frame(Var=rep("DistToShore",2), TimeOfDay=c("Day","Night"), 
                                         Fit=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit), 
                                         L95=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit - 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit)),
                                         U95=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit + 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit))))


#### Tidal Current
summary(dinoc$RMS_TIDAL_SPEED)
dinoc2 <- dinoc[which(dinoc$Diurnal != "Dusk" & dinoc$RMS_TIDAL_SPEED > -9999),]
hist(dinoc2$RMS_TIDAL_SPEED)
dinoc2$LnCurr <- log(dinoc2$RMS_TIDAL_SPEED + 0.001)
hist(dinoc2$LnCurr)

# fit and select models
mod.r.w <- lme(LnCurr ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.w.0 <- lme(LnCurr ~ 1, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r <- lme(LnCurr ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.0 <- lme(LnCurr ~ 1, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w <- lme(LnCurr ~ Diurnal, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w.0 <- lme(LnCurr ~ 1, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri <- lme(LnCurr ~ Diurnal, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.0 <- lme(LnCurr ~ 1, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w <- lme(LnCurr ~ Diurnal, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w.0 <- lme(LnCurr ~ 1, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr <- lme(LnCurr ~ Diurnal, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.0 <- lme(LnCurr ~ 1, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.w <- gls(LnCurr ~ Diurnal, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod.w.0 <- gls(LnCurr ~ 1, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod <- gls(LnCurr ~ Diurnal, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.0 <- gls(LnCurr ~ 1, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
modtab <- data.frame(modname=c("mod.r.w.0","mod.r","mod.r.0","mod.ri.w","mod.ri.w.0","mod.ri","mod.ri.0","mod.rr.w.0","mod.rr","mod.rr.0","mod","mod.0"))
modtab$AICc <- c(AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod), AICc(mod.0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod), AICc(mod.0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_DiurnalNocturnalDiff_TidalCurrents.csv")

# summarize best model
summary(mod.ri.w)
#                 Value  Std.Error   DF   t-value p-value
#(Intercept)  -6.884246 0.00950371 4076 -724.3746       0
#DiurnalNight  0.000000 0.00000000 4076    0.0000       1

dinocsumm <- rbind(dinocsumm, data.frame(Var=rep("TidalCurrent",2), TimeOfDay=c("Day","Night"), 
                                          Fit=exp(predictSE.lme(mod.ri.w, newdata, se.fit=T)$fit), 
                                          L95=exp(predictSE.lme(mod.ri.w, newdata, se.fit=T)$fit - 1.96*(predictSE.lme(mod.ri.w, newdata, se.fit=T)$se.fit)),
                                          U95=exp(predictSE.lme(mod.ri.w, newdata, se.fit=T)$fit + 1.96*(predictSE.lme(mod.ri.w, newdata, se.fit=T)$se.fit))))
                                         

#### Vessel density
summary(dinoc$RASTERVALU)
dinoc2 <- dinoc[which(dinoc$Diurnal != "Dusk"),]
hist(dinoc2$RASTERVALU)
dinoc2$LnDens <- log(dinoc2$RASTERVALU)
hist(dinoc2$LnDens)

# fit and select models
mod.r.w <- lme(LnDens ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.w.0 <- lme(LnDens ~ 1, random=list(~1|Region, ~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r <- lme(LnDens ~ Diurnal, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.r.0 <- lme(LnDens ~ 1, random=list(~1|Region, ~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w <- lme(LnDens ~ Diurnal, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.w.0 <- lme(LnDens ~ 1, random=list(~1|tag_local_identifier), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri <- lme(LnDens ~ Diurnal, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.ri.0 <- lme(LnDens ~ 1, random=list(~1|tag_local_identifier), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w <- lme(LnDens ~ Diurnal, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.w.0 <- lme(LnDens ~ 1, random=list(~1|Region), weights=varIdent(form=~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr <- lme(LnDens ~ Diurnal, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.rr.0 <- lme(LnDens ~ 1, random=list(~1|Region), data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.w <- gls(LnDens ~ Diurnal, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod.w.0 <- gls(LnDens ~ 1, data=dinoc2, weights=varIdent(form=~1|Region), method="ML", control = lmeControl(opt = "optim"))
mod <- gls(LnDens ~ Diurnal, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
mod.0 <- gls(LnDens ~ 1, data=dinoc2, method="ML", control = lmeControl(opt = "optim"))
modtab <- data.frame(modname=c("mod.r.w","mod.r.w.0","mod.r","mod.r.0","mod.ri.w","mod.ri.w.0","mod.ri","mod.ri.0","mod.rr.w","mod.rr.w.0","mod.rr","mod.rr.0","mod.w","mod.w.0","mod","mod.0"))
modtab$AICc <- c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0))
modtab$DAICc <- modtab$AICc-min(modtab$AICc)
modtab$Weight <- Weights(c(AICc(mod.r.w),AICc(mod.r.w.0), AICc(mod.r), AICc(mod.r.0), AICc(mod.ri.w),AICc(mod.ri.w.0), AICc(mod.ri), AICc(mod.ri.0), AICc(mod.rr.w),AICc(mod.rr.w.0), AICc(mod.rr), AICc(mod.rr.0), AICc(mod.w), AICc(mod.w.0), AICc(mod), AICc(mod.0)))
modtab$ER <- 1/(exp(-1/2 * modtab$DAICc))
modtab[order(modtab$DAICc),]
write.csv(modtab, file="Modtab_ER_DiurnalNocturnalDiff_VesselDensity.csv")

# summarize best model
summary(mod.r.w)
#                   Value  Std.Error   DF    t-value p-value
#(Intercept)  -2.6237773 0.20173286 4163 -13.006197  0.0000
#DiurnalNight -0.0188919 0.02090849 4163  -0.903552  0.3663

dinocsumm <- rbind(dinocsumm, data.frame(Var=rep("VesselDensity",2), TimeOfDay=c("Day","Night"), 
                                         Fit=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit), 
                                         L95=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit - 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit)),
                                         U95=exp(predictSE.lme(mod.r.w, newdata, se.fit=T)$fit + 1.96*(predictSE.lme(mod.r.w, newdata, se.fit=T)$se.fit))))



# save file
write.csv(dinocsumm, "DiurnalNocturnalCovariateSummary.csv")

## produce plots

# combine into a single plot
library(ggpubr)

# define color scheme - night/day
#colorlist <- c("#f0d848", "#303048")



colorlist <- c("#b8b6f6", "#120be0")
# depth
p.depth <- ggplot(dinocsumm[which(dinocsumm$Var=="Depth"),], aes(x=TimeOfDay, y=Fit, fill=TimeOfDay))+  
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=L95, ymax=U95), width=0.2)+
  scale_x_discrete(breaks=c("Day","Night"), labels=c("Day","Night"))+
  guides(fill=guide_legend(title="Time of Day"))+
  scale_fill_manual(values=colorlist)+
  labs(x="", y="Depth (m, with 95% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
print(p.depth)

colorlist <- c("#bbcccb", "#1a524d")
# distance to shore
p.dist <- ggplot(dinocsumm[which(dinocsumm$Var=="DistToShore"),], aes(x=TimeOfDay, y=Fit, fill=TimeOfDay))+  
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=L95, ymax=U95), width=0.2)+
  scale_x_discrete(breaks=c("Day","Night"), labels=c("Day","Night"))+
  guides(fill=guide_legend(title="Time of Day"))+
  scale_fill_manual(values=colorlist)+
  labs(x="", y="Distance to shore (m, with 95% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 14), 
        axis.title.x=element_text(size = 15))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
print(p.dist)

colorlist <- c("#d2cac2", "#694c32")
# tidal current
p.curr <- ggplot(dinocsumm[which(dinocsumm$Var=="TidalCurrent"),], aes(x=TimeOfDay, y=Fit, fill=TimeOfDay))+  
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=L95, ymax=U95), width=0.2)+
  scale_x_discrete(breaks=c("Day","Night"), labels=c("Day","Night"))+
  guides(fill=guide_legend(title="Time of Day"))+
  scale_fill_manual(values=colorlist)+
  labs(x="", y="Tidal Current (with 95% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
print(p.curr)

colorlist <- c("#deb8b7", "#91110d")
# vessel density
p.vess <- ggplot(dinocsumm[which(dinocsumm$Var=="VesselDensity"),], aes(x=TimeOfDay, y=Fit, fill=TimeOfDay))+  
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=L95, ymax=U95), width=0.2)+
  scale_x_discrete(breaks=c("Day","Night"), labels=c("Day","Night"))+
  guides(fill=guide_legend(title="Time of Day"))+
  scale_fill_manual(values=colorlist)+
  labs(x="", y="Vessel density (with 95% CI)") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y=element_text(size = 15), 
        axis.title.x=element_text(size = 15))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
print(p.vess)


# Arranging the plot
png(file="SurfScoter_DiurnalNocturnalPanel.png", width=7, height=9, units="in" , res=600)
p.all <- ggarrange(p.depth, p.dist, p.curr, p.vess, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 2), heights = c(2, 2),
          #labels = c("1","2","3","4"),
          common.legend = FALSE)
print(p.all)
dev.off()
