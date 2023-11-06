### Title: SingleSpOccupancy
### Author: Courtney L. Larson
### Description: Single species occupancy analysis for Larson et al. 2023. Reptile responses to outdoor recreation in urban habitat fragments.

#################### packages
#install.packages(reshape2, dplyr, unmarked, AICcmodavg)
library(reshape2)
library(dplyr)
library(unmarked)
library(AICcmodavg)

#################### read in data

setwd("C:/Courtney/PhD_work/HerpAnalysis/UrbanEcol submission/DataCode_forRepository")
herpcounts <- read.csv("All-sp_All-season_Dets.csv",header=T,sep=",") # this is a spreadsheet of all detections, with visit number, date, visit ID, species, point
pts <- read.csv("Sampling-points-list.csv",header=TRUE) # this is a list of all the sampling points
siteCovs <- read.csv("SiteCovs.csv",header=TRUE) #site-level covariates
obsCovs <- read.csv("ObsCovs.csv",header=TRUE) #observation-level covariates

#################### functions

#Goodness of fit test for the top model; Place the function you create (named fitstats here) in the bootstrapping code below
fitstats <- function(mod) {
  observed <- getY(mod@data)
  expected <- fitted(mod)
  resids <- residuals(mod)
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)}

# a helpful function for working with strings from https://stackoverflow.com/questions/13863599/insert-a-character-at-a-specific-location-in-a-string

substr1 <- function(x,y) {
  z <- sapply(strsplit(as.character(x),''),function(w) paste(na.omit(w[y]),collapse=''))
  dim(z) <- dim(x)
  return(z) }

`substr1<-` <- function(x,y,value) {
  names(y) <- c(value,rep('',length(y)-length(value)))
  z <- sapply(strsplit(as.character(x),''),function(w) {
    v <- seq(w)
    names(v) <- w
    paste(names(sort(c(y,v[setdiff(v,y)]))),collapse='') })
  dim(z) <- dim(x)
  return(z) }

#################### clean and prepare reptile detection data for occupancy modeling

# cleaning: drop problematic visits and sample points
herpcounts <- filter(herpcounts, Visit != "remove") # drop visits to remove (flagged as such in Access)
herpcounts <- filter(herpcounts, PointGrp != "BMT_P01" & PointGrp != "BMT_P02") # drop BMTP01 and BMTP02
herpcounts <- filter(herpcounts,!Visit %in% c("s5v3","s5v4"))

herpcounts$Occ <- rep(1, nrow(herpcounts))

# make separate dfs for 3 focal species
utst_counts <- filter(herpcounts,SpeciesCode == "UTST")
ashy_counts <- filter(herpcounts,SpeciesCode == "ASHY")
scoc_counts <- filter(herpcounts,SpeciesCode == "SCOC")

# create a dataframe with multiple repeats of each site name, for inserting data from each visit. This will allow inserting the appropriate zeroes for non-detections
n=11 #this is the number of visits
pts <- filter(pts,PointGrp != "BMT_P01" & PointGrp != "BMT_P02")
pts.seas.empty <- pts[rep(row.names(pts),n),]
pts.seas.empty <- pts.seas.empty[order(pts.seas.empty$PointGrp),] #sort the rows by PointGrp

# make a visit list and join to the site name df
visits.list.exp <- rep(c("s1v1","s1v2","s1v3","s2v1","s2v2","s2v3","s3v1","s3v2","s3v3","s5v1","s5v2"),nrow(pts))
pts.seas.empty$Visit <- as.factor(as.character(visits.list.exp))
pts.seas.empty$PointGrpVis <- paste(pts.seas.empty$PointGrp,pts.seas.empty$Visit,sep="-")

# join each species count data to PointGrp*Season df
utstcounts2 <- left_join(pts.seas.empty,utst_counts,by=c("PointGrp","Visit"))
ashycounts2 <- left_join(pts.seas.empty,ashy_counts,by=c("PointGrp","Visit"))
scoccounts2 <- left_join(pts.seas.empty,scoc_counts,by=c("PointGrp","Visit"))

utstcounts2$Season <- substr1(utstcounts2$Visit,c(1:2))
utstcounts2$VisNo <- substr1(utstcounts2$Visit,c(3:4))
ashycounts2$Season <- substr1(ashycounts2$Visit,c(1:2))
ashycounts2$VisNo <- substr1(ashycounts2$Visit,c(3:4))
scoccounts2$Season <- substr1(scoccounts2$Visit,c(1:2))
scoccounts2$VisNo <- substr1(scoccounts2$Visit,c(3:4))

# collapse multiple detections on a given visit
utstcounts3 <- utstcounts2 %>%
  group_by(PointGrpVis,PointGrp,Visit,Season,VisNo) %>%
  summarize(Occ2=sum(Occ))

ashycounts3 <- ashycounts2 %>%
  group_by(PointGrpVis,PointGrp,Visit,Season,VisNo) %>%
  summarize(Occ2=sum(Occ))

scoccounts3 <- scoccounts2 %>%
  group_by(PointGrpVis,PointGrp,Visit,Season,VisNo) %>%
  summarize(Occ2=sum(Occ))

# change NAs to 0s
b = which(utstcounts3$Occ2 > 0)
utstcounts3$Occ2[-b] = 0

b = which(ashycounts3$Occ2 > 0)
ashycounts3$Occ2[-b] = 0

b = which(scoccounts3$Occ2 > 0)
scoccounts3$Occ2[-b] = 0

# put in NAs for sites that were not visited for particular surveys
utstcounts3 <- as.data.frame(utstcounts3)
row.names(utstcounts3) <- make.names(utstcounts3$PointGrpVis,unique=TRUE)
utstcounts3["LPQ_P08.s3v1","Occ2"] = NA
utstcounts3["LPQ_P01.s3v2","Occ2"] = NA
utstcounts3["LPQ_P10.s3v3","Occ2"] = NA

ashycounts3 <- as.data.frame(ashycounts3)
row.names(ashycounts3) <- make.names(ashycounts3$PointGrpVis,unique=TRUE)
ashycounts3["LPQ_P08.s3v1","Occ2"] = NA
ashycounts3["LPQ_P01.s3v2","Occ2"] = NA
ashycounts3["LPQ_P10.s3v3","Occ2"] = NA

scoccounts3 <- as.data.frame(scoccounts3)
row.names(scoccounts3) <- make.names(scoccounts3$PointGrpVis,unique=TRUE)
scoccounts3["LPQ_P08.s3v1","Occ2"] = NA
scoccounts3["LPQ_P01.s3v2","Occ2"] = NA
scoccounts3["LPQ_P10.s3v3","Occ2"] = NA

# make a version in which counts are collapsed to 1 for occupancy
utstcounts4 <- utstcounts3
b = which(utstcounts4$Occ2 > 0)
utstcounts4$Occ2[b] = 1

ashycounts4 <- ashycounts3
b = which(ashycounts4$Occ2 > 0)
ashycounts4$Occ2[b] = 1

scoccounts4 <- scoccounts3
b = which(scoccounts4$Occ2 > 0)
scoccounts4$Occ2[b] = 1

# get detection histories into format for unmarked
utstDH <- reshape2::dcast(utstcounts4,PointGrp*Season ~ VisNo)
rownames(utstDH) <- paste(utstDH$PointGrp,utstDH$Season,sep="-")
utstM <- as.matrix(utstDH[,3:5])

ashyDH <- reshape2::dcast(ashycounts4,PointGrp*Season ~ VisNo)
rownames(ashyDH) <- paste(ashyDH$PointGrp,ashyDH$Season,sep="-")
ashyM <- as.matrix(ashyDH[,3:5])

scocDH <- reshape2::dcast(scoccounts4,PointGrp*Season ~ VisNo)
rownames(scocDH) <- paste(scocDH$PointGrp,scocDH$Season,sep="-")
scocM <- as.matrix(scocDH[,3:5])

#################### prepare covariate data

# scale the numeric covariates (site-level covariates)
ind <- sapply(siteCovs, is.numeric)
siteCovs[ind] <- lapply(siteCovs[ind], scale)


# deal with missing observations from the observation-level covariates for model selection purposes (see https://cran.r-project.org/web/packages/unmarked/vignettes/colext.pdf)
obsCovsNA <- matrix()
for (i in 1:12){
  l=i*3-2
  cols <- obsCovs[,(l+2):(l+4)]
  cols[is.na(utstM) != is.na(cols)] <- NA
  for (j in 1:12){
    k=j*3-2
    cols[is.na(obsCovs[,(k+2):(k+4)]) != is.na(cols)] <- NA
  }
  obsCovsNA <- cbind(obsCovsNA,cols)
}

#################### set up unmarked data frames for occupancy modeling

# UTST 
utst_umf <- unmarkedFrameOccu(y=utstM, 
                              siteCovs=siteCovs,
                              obsCovs=list(temp=obsCovsNA[,2:4],
                                           sky=obsCovsNA[,5:7],
                                           time=obsCovsNA[,8:10],
                                           jul=obsCovsNA[,11:13],
                                           obs=obsCovsNA[,14:16],
                                           obsn=obsCovsNA[,17:19],
                                           damp=obsCovsNA[,20:22],
                                           wind=obsCovsNA[,23:25],
                                           vtlen=obsCovsNA[,26:28],
                                           trwid=obsCovsNA[,32:34],
                                           seas=obsCovsNA[,35:37],
                                           hum=obsCovsNA[,29:31]))

utst_umf@obsCovs$temp <- scale(utst_umf@obsCovs$temp)
utst_umf@obsCovs$time <- scale(utst_umf@obsCovs$time)
utst_umf@obsCovs$jul <- scale(utst_umf@obsCovs$jul)
utst_umf@obsCovs$obsn <- as.numeric(utst_umf@obsCovs$obsn)
utst_umf@obsCovs$obsn <- scale(utst_umf@obsCovs$obsn)
utst_umf@obsCovs$wind <- as.numeric(utst_umf@obsCovs$wind)
utst_umf@obsCovs$wind <- scale(utst_umf@obsCovs$wind)
utst_umf@obsCovs$hum <- scale(utst_umf@obsCovs$hum)

ashy_umf <- unmarkedFrameOccu(y=ashyM, 
                              siteCovs=siteCovs,
                              obsCovs=list(temp=obsCovsNA[,2:4],
                                           sky=obsCovsNA[,5:7],
                                           time=obsCovsNA[,8:10],
                                           jul=obsCovsNA[,11:13],
                                           obs=obsCovsNA[,14:16],
                                           obsn=obsCovsNA[,17:19],
                                           damp=obsCovsNA[,20:22],
                                           wind=obsCovsNA[,23:25],
                                           vtlen=obsCovsNA[,26:28],
                                           trwid=obsCovsNA[,32:34],
                                           seas=obsCovsNA[,35:37],
                                           hum=obsCovsNA[,29:31]))

ashy_umf@obsCovs$temp <- scale(ashy_umf@obsCovs$temp)
ashy_umf@obsCovs$time <- scale(ashy_umf@obsCovs$time)
ashy_umf@obsCovs$jul <- scale(ashy_umf@obsCovs$jul)
ashy_umf@obsCovs$obsn <- as.numeric(ashy_umf@obsCovs$obsn)
ashy_umf@obsCovs$obsn <- scale(ashy_umf@obsCovs$obsn)
ashy_umf@obsCovs$wind <- as.numeric(ashy_umf@obsCovs$wind)
ashy_umf@obsCovs$wind <- scale(ashy_umf@obsCovs$wind)
ashy_umf@obsCovs$hum <- scale(ashy_umf@obsCovs$hum)

scoc_umf <- unmarkedFrameOccu(y=scocM, 
                              siteCovs=siteCovs,
                              obsCovs=list(temp=obsCovsNA[,2:4],
                                           sky=obsCovsNA[,5:7],
                                           time=obsCovsNA[,8:10],
                                           jul=obsCovsNA[,11:13],
                                           obs=obsCovsNA[,14:16],
                                           obsn=obsCovsNA[,17:19],
                                           damp=obsCovsNA[,20:22],
                                           wind=obsCovsNA[,23:25],
                                           vtlen=obsCovsNA[,26:28],
                                           trwid=obsCovsNA[,32:34],
                                           seas=obsCovsNA[,35:37],
                                           hum=obsCovsNA[,29:31]))

scoc_umf@obsCovs$temp <- scale(scoc_umf@obsCovs$temp)
scoc_umf@obsCovs$time <- scale(scoc_umf@obsCovs$time)
scoc_umf@obsCovs$jul <- scale(scoc_umf@obsCovs$jul)
scoc_umf@obsCovs$obsn <- as.numeric(scoc_umf@obsCovs$obsn)
scoc_umf@obsCovs$obsn <- scale(scoc_umf@obsCovs$obsn)
scoc_umf@obsCovs$wind <- as.numeric(scoc_umf@obsCovs$wind)
scoc_umf@obsCovs$wind <- scale(scoc_umf@obsCovs$wind)
scoc_umf@obsCovs$hum <- scale(scoc_umf@obsCovs$hum)

#################### run occupancy models for ASHY

## test null model to get ballpark for p and psi estimates
ashy.null <- occu(~1 ~1, data=ashy_umf)
summary(ashy.null)

## gof test
#ashy.global <- occu(~temp+I(temp^2)+jul+I(jul^2)+obsn ~mds1+mds2+fire+ndvi+elev10+chili10+hik.ln, data=ashy_umf)
#ashy.obs.boot <- mb.gof.test(ashy.global,nsim=1000) 
#ashy.obs.boot

## models w/varying detection
ashyp1 <- occu(~sky ~mds1+mds2+fire, data=ashy_umf)
ashyp2 <- occu(~wind ~mds1+mds2+fire, data=ashy_umf)
ashyp3 <- occu(~temp+I(temp^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp4 <- occu(~jul+I(jul^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp5 <- occu(~time+I(time^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp6 <- occu(~obs ~mds1+mds2+fire, data=ashy_umf)
ashyp7 <- occu(~trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp8 <- occu(~obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp9 <- occu(~hum ~mds1+mds2+fire, data=ashy_umf)

ashyp10 <- occu(~sky+jul+I(jul^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp11 <- occu(~sky+time+I(time^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp12 <- occu(~sky+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp13 <- occu(~sky+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp14 <- occu(~sky+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp15 <- occu(~sky+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp16 <- occu(~wind+jul+I(jul^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp17 <- occu(~wind+time+I(time^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp18 <- occu(~wind+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp19 <- occu(~wind+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp20 <- occu(~wind+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp21 <- occu(~wind+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp22 <- occu(~temp+I(temp^2)+jul+I(jul^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp23 <- occu(~temp+I(temp^2)+time+I(time^2) ~mds1+mds2+fire, data=ashy_umf)
ashyp24 <- occu(~temp+I(temp^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp25 <- occu(~temp+I(temp^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp26 <- occu(~temp+I(temp^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp27 <- occu(~temp+I(temp^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp28 <- occu(~jul+I(jul^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp29 <- occu(~jul+I(jul^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp30 <- occu(~jul+I(jul^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp31 <- occu(~jul+I(jul^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp32 <- occu(~time+I(time^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp33 <- occu(~time+I(time^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp34 <- occu(~time+I(time^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp35 <- occu(~time+I(time^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp36 <- occu(~obs+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp37 <- occu(~trwid+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp38 <- occu(~obsn+hum ~mds1+mds2+fire, data=ashy_umf)

ashyp39 <- occu(~sky+jul+I(jul^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp40 <- occu(~sky+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp41 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp42 <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp43 <- occu(~sky+time+I(time^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp44 <- occu(~sky+time+I(time^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp45 <- occu(~sky+time+I(time^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp46 <- occu(~sky+time+I(time^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp47 <- occu(~wind+jul+I(jul^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp48 <- occu(~wind+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp49 <- occu(~wind+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp50 <- occu(~wind+jul+I(jul^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp51 <- occu(~wind+time+I(time^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp52 <- occu(~wind+time+I(time^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp53 <- occu(~wind+time+I(time^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp54 <- occu(~wind+time+I(time^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp55 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp56 <- occu(~temp+I(temp^2)+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp57 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp58 <- occu(~temp+I(temp^2)+jul+I(jul^2)+hum ~mds1+mds2+fire, data=ashy_umf)
ashyp59 <- occu(~temp+I(temp^2)+time+I(time^2)+obs ~mds1+mds2+fire, data=ashy_umf)
ashyp60 <- occu(~temp+I(temp^2)+time+I(time^2)+trwid ~mds1+mds2+fire, data=ashy_umf)
ashyp61 <- occu(~temp+I(temp^2)+time+I(time^2)+obsn ~mds1+mds2+fire, data=ashy_umf)
ashyp62 <- occu(~temp+I(temp^2)+time+I(time^2)+hum ~mds1+mds2+fire, data=ashy_umf)

ashy.p.models <- list(ashyp1,ashyp2,ashyp3,ashyp4,ashyp5,ashyp6,ashyp7,ashyp8,ashyp9,ashyp10,ashyp11,ashyp12,ashyp13,ashyp14,ashyp15,ashyp16,ashyp17,ashyp18,ashyp19,ashyp20,ashyp21,ashyp22,ashyp23,ashyp24,ashyp25,ashyp26,ashyp27,ashyp28,ashyp29,ashyp30,ashyp31,ashyp32,ashyp33,ashyp34,ashyp35,ashyp36,ashyp37,ashyp38,ashyp39,ashyp40,ashyp41,ashyp42,ashyp43,ashyp44,ashyp45,ashyp46,ashyp47,ashyp48,ashyp49,ashyp50,ashyp51,ashyp52,ashyp53,ashyp54,ashyp55,ashyp56,ashyp57,ashyp58,ashyp59,ashyp60,ashyp61,ashyp62)
ashy.p.aic <- aictab(ashy.p.models,nobs=92,c.hat=1.74)

## models w/varying psi
ashyt1 <- occu(~sky+jul+I(jul^2)+obsn ~mds1, data=ashy_umf)
ashyt2 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2, data=ashy_umf)
ashyt3 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+fire, data=ashy_umf)
ashyt4 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+ndvi, data=ashy_umf)
ashyt5 <- occu(~sky+jul+I(jul^2)+obsn ~mds2, data=ashy_umf)
ashyt6 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+fire, data=ashy_umf)
ashyt7 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+ndvi, data=ashy_umf)
ashyt8 <- occu(~sky+jul+I(jul^2)+obsn ~fire, data=ashy_umf)
ashyt9 <- occu(~sky+jul+I(jul^2)+obsn ~fire+ndvi, data=ashy_umf)
ashyt10 <- occu(~sky+jul+I(jul^2)+obsn ~ndvi, data=ashy_umf)

ashyt11 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+elev10, data=ashy_umf)
ashyt12 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+elev10, data=ashy_umf)
ashyt13 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+fire+elev10, data=ashy_umf)
ashyt14 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+ndvi+elev10, data=ashy_umf)
ashyt15 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+elev10, data=ashy_umf)
ashyt16 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+fire+elev10, data=ashy_umf)
ashyt17 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+ndvi+elev10, data=ashy_umf)
ashyt18 <- occu(~sky+jul+I(jul^2)+obsn ~fire+elev10, data=ashy_umf)
ashyt19 <- occu(~sky+jul+I(jul^2)+obsn ~fire+ndvi+elev10, data=ashy_umf)
ashyt20 <- occu(~sky+jul+I(jul^2)+obsn ~ndvi+elev10, data=ashy_umf)

ashyt21 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+chili10, data=ashy_umf)
ashyt22 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+chili10, data=ashy_umf)
ashyt23 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+fire+chili10, data=ashy_umf)
ashyt24 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+ndvi+chili10, data=ashy_umf)
ashyt25 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+chili10, data=ashy_umf)
ashyt26 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+fire+chili10, data=ashy_umf)
ashyt27 <- occu(~sky+jul+I(jul^2)+obsn ~mds2+ndvi+chili10, data=ashy_umf)
ashyt28 <- occu(~sky+jul+I(jul^2)+obsn ~fire+chili10, data=ashy_umf)
ashyt29 <- occu(~sky+jul+I(jul^2)+obsn ~fire+ndvi+chili10, data=ashy_umf)
ashyt30 <- occu(~sky+jul+I(jul^2)+obsn ~ndvi+chili10, data=ashy_umf)

ashyt1b <- occu(~sky+jul+I(jul^2) ~mds1, data=ashy_umf)
ashyt2b <- occu(~sky+jul+I(jul^2) ~mds1+mds2, data=ashy_umf)
ashyt3b <- occu(~sky+jul+I(jul^2) ~mds1+fire, data=ashy_umf)
ashyt4b <- occu(~sky+jul+I(jul^2) ~mds1+ndvi, data=ashy_umf)
ashyt5b <- occu(~sky+jul+I(jul^2) ~mds2, data=ashy_umf)
ashyt6b <- occu(~sky+jul+I(jul^2) ~mds2+fire, data=ashy_umf)
ashyt7b <- occu(~sky+jul+I(jul^2) ~mds2+ndvi, data=ashy_umf)
ashyt8b <- occu(~sky+jul+I(jul^2) ~fire, data=ashy_umf)
ashyt9b <- occu(~sky+jul+I(jul^2) ~fire+ndvi, data=ashy_umf)
ashyt10b<- occu(~sky+jul+I(jul^2) ~ndvi, data=ashy_umf)

ashyt11b <- occu(~sky+jul+I(jul^2) ~mds1+elev10, data=ashy_umf)
ashyt12b <- occu(~sky+jul+I(jul^2) ~mds1+mds2+elev10, data=ashy_umf)
ashyt13b <- occu(~sky+jul+I(jul^2) ~mds1+fire+elev10, data=ashy_umf)
ashyt14b <- occu(~sky+jul+I(jul^2) ~mds1+ndvi+elev10, data=ashy_umf)
ashyt15b <- occu(~sky+jul+I(jul^2) ~mds2+elev10, data=ashy_umf)
ashyt16b <- occu(~sky+jul+I(jul^2) ~mds2+fire+elev10, data=ashy_umf)
ashyt17b <- occu(~sky+jul+I(jul^2) ~mds2+ndvi+elev10, data=ashy_umf)
ashyt18b <- occu(~sky+jul+I(jul^2) ~fire+elev10, data=ashy_umf)
ashyt19b <- occu(~sky+jul+I(jul^2) ~fire+ndvi+elev10, data=ashy_umf)
ashyt20b <- occu(~sky+jul+I(jul^2) ~ndvi+elev10, data=ashy_umf)

ashyt21b <- occu(~sky+jul+I(jul^2) ~mds1+chili10, data=ashy_umf)
ashyt22b <- occu(~sky+jul+I(jul^2) ~mds1+mds2+chili10, data=ashy_umf)
ashyt23b <- occu(~sky+jul+I(jul^2) ~mds1+fire+chili10, data=ashy_umf)
ashyt24b <- occu(~sky+jul+I(jul^2) ~mds1+ndvi+chili10, data=ashy_umf)
ashyt25b <- occu(~sky+jul+I(jul^2) ~mds2+chili10, data=ashy_umf)
ashyt26b <- occu(~sky+jul+I(jul^2) ~mds2+fire+chili10, data=ashy_umf)
ashyt27b <- occu(~sky+jul+I(jul^2) ~mds2+ndvi+chili10, data=ashy_umf)
ashyt28b <- occu(~sky+jul+I(jul^2) ~fire+chili10, data=ashy_umf)
ashyt29b <- occu(~sky+jul+I(jul^2) ~fire+ndvi+chili10, data=ashy_umf)
ashyt30b <- occu(~sky+jul+I(jul^2) ~ndvi+chili10, data=ashy_umf)

ashy.t.models <- list("ashyt1"=ashyt1,"ashyt2"=ashyt2,"ashyt3"=ashyt3,"ashyt4"=ashyt4,"ashyt5"=ashyt5,"ashyt6"=ashyt6,"ashyt7"=ashyt7,"ashyt8"=ashyt8,"ashyt9"=ashyt9,"ashyt10"=ashyt10,"ashyt11"=ashyt11,"ashyt12"=ashyt12,"ashyt13"=ashyt13,"ashyt14"=ashyt14,"ashyt15"=ashyt15,"ashyt16"=ashyt16,"ashyt17"=ashyt17,"ashyt18"=ashyt18,"ashyt19"=ashyt19,"ashyt20"=ashyt20,"ashyt21"=ashyt21,"ashyt22"=ashyt22,"ashyt23"=ashyt23,"ashyt24"=ashyt24,"ashyt25"=ashyt25,"ashyt26"=ashyt26,"ashyt27"=ashyt27,"ashyt28"=ashyt28,"ashyt29"=ashyt29,"ashyt30"=ashyt30,"ashyt1b"=ashyt1b,"ashyt2b"=ashyt2b,"ashyt3b"=ashyt3b,"ashyt4b"=ashyt4b,"ashyt5b"=ashyt5b,"ashyt6b"=ashyt6b,"ashyt7b"=ashyt7b,"ashyt8b"=ashyt8b,"ashyt9b"=ashyt9b,"ashyt10b"=ashyt10b,"ashyt11b"=ashyt11b,"ashyt12b"=ashyt12b,"ashyt13b"=ashyt13b,"ashyt14b"=ashyt14b,"ashyt15b"=ashyt15b,"ashyt16b"=ashyt16b,"ashyt17b"=ashyt17b,"ashyt18b"=ashyt18b,"ashyt19b"=ashyt19b,"ashyt20b"=ashyt20b,"ashyt21b"=ashyt21b,"ashyt22b"=ashyt22b,"ashyt23b"=ashyt23b,"ashyt24b"=ashyt24b,"ashyt25b"=ashyt25b,"ashyt26b"=ashyt26b,"ashyt27b"=ashyt27b,"ashyt28b"=ashyt28b,"ashyt29b"=ashyt29b,"ashyt30b"=ashyt30b)
ashy.t.aic <- aictab(ashy.t.models,nobs=92,c.hat=1.74)

ashytop1 <- ashyt5b
ashytop2 <- ashyt2b
ashytop3 <- ashyt5
ashytop4 <- ashyt1b
ashytop5 <- ashyt2
ashytop6 <- ashyt8b
ashytop7 <- ashyt12b
ashytop8 <- ashyt6b

## adding human vars to psi 
ashytop1.cyc <- occu(~sky+jul+I(jul^2) ~mds2+cyc.ln, data=ashy_umf)
ashytop1.hik <- occu(~sky+jul+I(jul^2) ~mds2+hik.ln, data=ashy_umf)
ashytop1.hum <- occu(~sky+jul+I(jul^2) ~mds2+allhum.ln, data=ashy_umf)

ashytop2.cyc <- occu(~sky+jul+I(jul^2) ~mds1+mds2+cyc.ln, data=ashy_umf)
ashytop2.hik <- occu(~sky+jul+I(jul^2) ~mds1+mds2+hik.ln, data=ashy_umf)
ashytop2.hum <- occu(~sky+jul+I(jul^2) ~mds1+mds2+allhum.ln, data=ashy_umf)

ashytop3.cyc <- occu(~sky+jul+I(jul^2)+obsn ~mds2+cyc.ln, data=ashy_umf)
ashytop3.hik <- occu(~sky+jul+I(jul^2)+obsn ~mds2+hik.ln, data=ashy_umf)
ashytop3.hum <- occu(~sky+jul+I(jul^2)+obsn ~mds2+allhum.ln, data=ashy_umf)

ashytop4.cyc <- occu(~sky+jul+I(jul^2) ~mds1+cyc.ln, data=ashy_umf)
ashytop4.hik <- occu(~sky+jul+I(jul^2) ~mds1+hik.ln, data=ashy_umf)
ashytop4.hum <- occu(~sky+jul+I(jul^2) ~mds1+allhum.ln, data=ashy_umf)

ashytop5.cyc <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+cyc.ln, data=ashy_umf)
ashytop5.hik <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+hik.ln, data=ashy_umf)
ashytop5.hum <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+allhum.ln, data=ashy_umf)

ashytop6.cyc <- occu(~sky+jul+I(jul^2) ~fire+cyc.ln, data=ashy_umf)
ashytop6.hik <- occu(~sky+jul+I(jul^2) ~fire+hik.ln, data=ashy_umf)
ashytop6.hum <- occu(~sky+jul+I(jul^2) ~fire+allhum.ln, data=ashy_umf)

ashytop7.cyc <- occu(~sky+jul+I(jul^2) ~mds1+mds2+elev10+cyc.ln, data=ashy_umf)
ashytop7.hik <- occu(~sky+jul+I(jul^2) ~mds1+mds2+elev10+hik.ln, data=ashy_umf)
ashytop7.hum <-occu(~sky+jul+I(jul^2) ~mds1+mds2+elev10+allhum.ln, data=ashy_umf) 

ashytop8.cyc <- occu(~sky+jul+I(jul^2) ~mds2+fire+cyc.ln, data=ashy_umf)
ashytop8.hik <- occu(~sky+jul+I(jul^2) ~mds2+fire+hik.ln, data=ashy_umf)
ashytop8.hum <- occu(~sky+jul+I(jul^2) ~mds2+fire+allhum.ln, data=ashy_umf)

## adding human vars to both p and psi
ashytop1.cycp <- occu(~sky+jul+I(jul^2)+hum ~mds2+cyc.ln, data=ashy_umf)
ashytop1.hikp <- occu(~sky+jul+I(jul^2)+hum ~mds2+hik.ln, data=ashy_umf)
ashytop1.hump <- occu(~sky+jul+I(jul^2)+hum ~mds2+allhum.ln, data=ashy_umf)
ashytop2.cycp <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+cyc.ln, data=ashy_umf)
ashytop2.hikp <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+hik.ln, data=ashy_umf)
ashytop2.hump <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+allhum.ln, data=ashy_umf)
ashytop3.cycp <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds2+cyc.ln, data=ashy_umf)
ashytop3.hikp <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds2+hik.ln, data=ashy_umf)
ashytop3.hump <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds2+allhum.ln, data=ashy_umf)
ashytop4.cycp <- occu(~sky+jul+I(jul^2)+hum ~mds1+cyc.ln, data=ashy_umf)
ashytop4.hikp <- occu(~sky+jul+I(jul^2)+hum ~mds1+hik.ln, data=ashy_umf)
ashytop4.hump <- occu(~sky+jul+I(jul^2)+hum ~mds1+allhum.ln, data=ashy_umf)
ashytop5.cycp <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds1+mds2+cyc.ln, data=ashy_umf)
ashytop5.hikp <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds1+mds2+hik.ln, data=ashy_umf)
ashytop5.hump <- occu(~sky+jul+I(jul^2)+obsn+hum ~mds1+mds2+allhum.ln, data=ashy_umf)
ashytop6.cycp <- occu(~sky+jul+I(jul^2)+hum ~fire+cyc.ln, data=ashy_umf)
ashytop6.hikp <- occu(~sky+jul+I(jul^2)+hum ~fire+hik.ln, data=ashy_umf)
ashytop6.hump <- occu(~sky+jul+I(jul^2)+hum ~fire+allhum.ln, data=ashy_umf)
ashytop7.cycp <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+elev10+cyc.ln, data=ashy_umf)
ashytop7.hikp <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+elev10+hik.ln, data=ashy_umf)
ashytop7.hump <- occu(~sky+jul+I(jul^2)+hum ~mds1+mds2+elev10+allhum.ln, data=ashy_umf) 
ashytop8.cycp <- occu(~sky+jul+I(jul^2)+hum ~mds2+fire+cyc.ln, data=ashy_umf)
ashytop8.hikp <- occu(~sky+jul+I(jul^2)+hum ~mds2+fire+hik.ln, data=ashy_umf)
ashytop8.hump <- occu(~sky+jul+I(jul^2)+hum ~mds2+fire+allhum.ln, data=ashy_umf)

ashy.v.models <- list("ashytop6.cyc"=ashytop6.cyc,
                      "ashytop6.hik"=ashytop6.hik,
                      "ashytop6.hum"=ashytop6.hum,
                      "ashytop6"=ashytop6,
                      "ashytop6.cycp"=ashytop6.cycp,
                      "ashytop6.hikp"=ashytop6.hikp,
                      "ashytop6.hump"=ashytop6.hump,
                      "ashytop4.cyc"=ashytop4.cyc,
                      "ashytop4.hik"=ashytop4.hik,
                      "ashytop4.hum"=ashytop4.hum,
                      "ashytop4"=ashytop4,
                      "ashytop4.cycp"=ashytop4.cycp,
                      "ashytop4.hikp"=ashytop4.hikp,
                      "ashytop4.hump"=ashytop4.hump,
                      "ashytop1.cyc"=ashytop1.cyc,
                      "ashytop1.hik"=ashytop1.hik,
                      "ashytop1.hum"=ashytop1.hum,
                      "ashytop1"=ashytop1,
                      "ashytop1.cycp"=ashytop1.cycp,
                      "ashytop1.hikp"=ashytop1.hikp,
                      "ashytop1.hump"=ashytop1.hump,
                      "ashytop3.cyc"=ashytop3.cyc,
                      "ashytop3.hik"=ashytop3.hik,
                      "ashytop3.hum"=ashytop3.hum,
                      "ashytop3"=ashytop3,
                      "ashytop3.cycp"=ashytop3.cycp,
                      "ashytop3.hikp"=ashytop3.hikp,
                      "ashytop3.hump"=ashytop3.hump,
                      "ashytop8.cyc"=ashytop8.cyc,
                      "ashytop8.hik"=ashytop8.hik,
                      "ashytop8.hum"=ashytop8.hum,
                      "ashytop8"=ashytop8,
                      "ashytop8.cycp"=ashytop8.cycp,
                      "ashytop8.hikp"=ashytop8.hikp,
                      "ashytop8.hump"=ashytop8.hump,
                      "ashytop2.cyc"=ashytop2.cyc,
                      "ashytop2.hik"=ashytop2.hik,
                      "ashytop2.hum"=ashytop2.hum,
                      "ashytop2"=ashytop2,
                      "ashytop2.cycp"=ashytop2.cycp,
                      "ashytop2.hikp"=ashytop2.hikp,
                      "ashytop2.hump"=ashytop2.hump,
                      "ashytop5.cyc"=ashytop5.cyc,
                      "ashytop5.hik"=ashytop5.hik,
                      "ashytop5.hum"=ashytop5.hum,
                      "ashytop5"=ashytop5,
                      "ashytop5.cycp"=ashytop5.cycp,
                      "ashytop5.hikp"=ashytop5.hikp,
                      "ashytop5.hump"=ashytop5.hump,
                      "ashytop7.cyc"=ashytop7.cyc,
                      "ashytop7.hik"=ashytop7.hik,
                      "ashytop7.hum"=ashytop7.hum,
                      "ashytop7"=ashytop7,
                      "ashytop7.cycp"=ashytop7.cycp,
                      "ashytop7.hikp"=ashytop7.hikp,
                      "ashytop7.hump"=ashytop7.hump,
                      "(.)p(.)"=ashy.null)

ashy.v.aic <- aictab(ashy.v.models,nobs=92,c.hat=1.74)

ashytophik <- ashytop1.hik
ashytopcyc <- ashytop1.cyc
ashytophum <- ashytop1.hum

ashy.topmds1 <- ashytop2
ashy.topmds2 <- ashytop1
ashy.topfire <- ashytop6
ashy.topelev <- ashytop7

confint(ashytophik,type="state")
confint(ashytopcyc,type="state")
confint(ashytophum,type="state")


#################### run occupancy models for UTST

## test null model to get ballpark for p and psi estimates
utst.null <- occu(~1 ~1, data=utst_umf)

#gof test
#utst.global <- occu(~sky+obs ~mds1+mds2+fire+ndvi+elev10+chili10+hik.ln, data=utst_umf)
#utst.obs.boot <- mb.gof.test(utst.global,nsim=1000) 
#utst.obs.boot

utstp1 <- occu(~sky ~mds1+mds2+fire, data=utst_umf)
utstp2 <- occu(~wind ~mds1+mds2+fire, data=utst_umf)
utstp3 <- occu(~temp+I(temp^2) ~mds1+mds2+fire, data=utst_umf)
utstp4 <- occu(~jul+I(jul^2) ~mds1+mds2+fire, data=utst_umf)
utstp5 <- occu(~time+I(time^2) ~mds1+mds2+fire, data=utst_umf)
utstp6 <- occu(~obs ~mds1+mds2+fire, data=utst_umf)
utstp7 <- occu(~trwid ~mds1+mds2+fire, data=utst_umf)
utstp8 <- occu(~obsn ~mds1+mds2+fire, data=utst_umf)

utstp9 <- occu(~sky+jul+I(jul^2) ~mds1+mds2+fire, data=utst_umf)
utstp10 <- occu(~sky+time+I(time^2) ~mds1+mds2+fire, data=utst_umf)
utstp11 <- occu(~sky+obs ~mds1+mds2+fire, data=utst_umf)
utstp12 <- occu(~sky+trwid ~mds1+mds2+fire, data=utst_umf)
utstp13 <- occu(~sky+obsn ~mds1+mds2+fire, data=utst_umf)

utstp14 <- occu(~wind+jul+I(jul^2) ~mds1+mds2+fire, data=utst_umf)
utstp15 <- occu(~wind+time+I(time^2) ~mds1+mds2+fire, data=utst_umf)
utstp16 <- occu(~wind+obs ~mds1+mds2+fire, data=utst_umf)
utstp17 <- occu(~wind+trwid ~mds1+mds2+fire, data=utst_umf)
utstp18 <- occu(~wind+obsn ~mds1+mds2+fire, data=utst_umf)

utstp19 <- occu(~temp+I(temp^2)+jul+I(jul^2) ~mds1+mds2+fire, data=utst_umf)
utstp20 <- occu(~temp+I(temp^2)+time+I(time^2) ~mds1+mds2+fire, data=utst_umf)
utstp21 <- occu(~temp+I(temp^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp22 <- occu(~temp+I(temp^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp23 <- occu(~temp+I(temp^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp24 <- occu(~jul+I(jul^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp25 <- occu(~jul+I(jul^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp26 <- occu(~jul+I(jul^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp27 <- occu(~time+I(time^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp28 <- occu(~time+I(time^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp29 <- occu(~time+I(time^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp30 <- occu(~sky+jul+I(jul^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp31 <- occu(~sky+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp32 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp33 <- occu(~sky+time+I(time^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp34 <- occu(~sky+time+I(time^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp35 <- occu(~sky+time+I(time^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp36 <- occu(~wind+jul+I(jul^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp37 <- occu(~wind+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp38 <- occu(~wind+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp39 <- occu(~wind+time+I(time^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp40 <- occu(~wind+time+I(time^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp41 <- occu(~wind+time+I(time^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp42 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp43 <- occu(~temp+I(temp^2)+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp44 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utstp45 <- occu(~temp+I(temp^2)+time+I(time^2)+obs ~mds1+mds2+fire, data=utst_umf)
utstp46 <- occu(~temp+I(temp^2)+time+I(time^2)+trwid ~mds1+mds2+fire, data=utst_umf)
utstp47 <- occu(~temp+I(temp^2)+time+I(time^2)+obsn ~mds1+mds2+fire, data=utst_umf)

utst.p.models <- list(utstp1,utstp2,utstp3,utstp4,utstp5,utstp6,utstp7,utstp8,utstp9,utstp10,utstp11,utstp12,utstp13,utstp14,utstp15,utstp16,utstp17,utstp18,utstp19,utstp20,utstp21,utstp22,utstp23,utstp24,utstp25,utstp26,utstp27,utstp28,utstp29,utstp30,utstp31,utstp32,utstp33,utstp34,utstp35,utstp36,utstp37,utstp38,utstp39,utstp40,utstp41,utstp42,utstp43,utstp44,utstp45,utstp46,utstp47)

utst.p.aic <- aictab(utst.p.models,nobs=92,c.hat=1.74)

utstt1 <- occu(~sky+obs ~mds1, data=utst_umf)
utstt2 <- occu(~sky+obs ~mds1+mds2, data=utst_umf)
utstt3 <- occu(~sky+obs ~mds1+fire, data=utst_umf)
utstt4 <- occu(~sky+obs ~mds1+ndvi, data=utst_umf)
utstt5 <- occu(~sky+obs ~mds2, data=utst_umf)
utstt6 <- occu(~sky+obs ~mds2+fire, data=utst_umf)
utstt7 <- occu(~sky+obs ~mds2+ndvi, data=utst_umf)
utstt8 <- occu(~sky+obs ~fire, data=utst_umf)
utstt9 <- occu(~sky+obs ~fire+ndvi, data=utst_umf)
utstt10 <- occu(~sky+obs ~ndvi, data=utst_umf)

utstt11 <- occu(~sky+obs ~mds1+elev10, data=utst_umf)
utstt12 <- occu(~sky+obs ~mds1+mds2+elev10, data=utst_umf)
utstt13 <- occu(~sky+obs ~mds1+fire+elev10, data=utst_umf)
utstt14 <- occu(~sky+obs ~mds1+ndvi+elev10, data=utst_umf)
utstt15 <- occu(~sky+obs ~mds2+elev10, data=utst_umf)
utstt16 <- occu(~sky+obs ~mds2+fire+elev10, data=utst_umf)
utstt17 <- occu(~sky+obs ~mds2+ndvi+elev10, data=utst_umf)
utstt18 <- occu(~sky+obs ~fire+elev10, data=utst_umf)
utstt19 <- occu(~sky+obs ~fire+ndvi+elev10, data=utst_umf)
utstt20 <- occu(~sky+obs ~ndvi+elev10, data=utst_umf)

utstt21 <- occu(~sky+obs ~mds1+chili10, data=utst_umf)
utstt22 <- occu(~sky+obs ~mds1+mds2+chili10, data=utst_umf)
utstt23 <- occu(~sky+obs ~mds1+fire+chili10, data=utst_umf)
utstt24 <- occu(~sky+obs ~mds1+ndvi+chili10, data=utst_umf)
utstt25 <- occu(~sky+obs ~mds2+chili10, data=utst_umf)
utstt26 <- occu(~sky+obs ~mds2+fire+chili10, data=utst_umf)
utstt27 <- occu(~sky+obs ~mds2+ndvi+chili10, data=utst_umf)
utstt28 <- occu(~sky+obs ~fire+chili10, data=utst_umf)
utstt29 <- occu(~sky+obs ~fire+ndvi+chili10, data=utst_umf)
utstt30 <- occu(~sky+obs ~ndvi+chili10, data=utst_umf)

utst.t.models <- list(utstt1,utstt2,utstt3,utstt4,utstt5,utstt6,utstt7,utstt8,utstt9,utstt10,utstt11,utstt12,utstt13,utstt14,utstt15,utstt16,utstt17,utstt18,utstt19,utstt20,utstt21,utstt22,utstt23,utstt24,utstt25,utstt26,utstt27,utstt28,utstt29,utstt30)
utst.t.aic <- aictab(utst.t.models,nobs=92,c.hat=1.74)

utsttop1 <- utstt26
utsttop2 <- utstt28
utsttop3 <- utstt16
utsttop4 <- utstt6
utsttop5 <- utstt8
utsttop6 <- utstt18

## adding human vars to psi
utsttop1.cyc <- occu(~sky+obs ~mds2+fire+chili10+cyc.ln, data=utst_umf)
utsttop1.hik <- occu(~sky+obs ~mds2+fire+chili10+hik.ln, data=utst_umf)
utsttop1.hum <- occu(~sky+obs ~mds2+fire+chili10+allhum.ln, data=utst_umf)

utsttop2.cyc <- occu(~sky+obs ~fire+chili10+cyc.ln, data=utst_umf)
utsttop2.hik <- occu(~sky+obs ~fire+chili10+hik.ln, data=utst_umf)
utsttop2.hum <- occu(~sky+obs ~fire+chili10+allhum.ln, data=utst_umf)

utsttop3.cyc <- occu(~sky+obs ~mds2+fire+elev10+cyc.ln, data=utst_umf)
utsttop3.hik <- occu(~sky+obs ~mds2+fire+elev10+hik.ln, data=utst_umf)
utsttop3.hum <- occu(~sky+obs ~mds2+fire+elev10+allhum.ln, data=utst_umf)

utsttop4.cyc <- occu(~sky+obs ~mds2+fire+cyc.ln, data=utst_umf)
utsttop4.hik <- occu(~sky+obs ~mds2+fire+hik.ln, data=utst_umf)
utsttop4.hum <- occu(~sky+obs ~mds2+fire+allhum.ln, data=utst_umf)

utsttop5.cyc <- occu(~sky+obs ~fire+cyc.ln, data=utst_umf)
utsttop5.hik <- occu(~sky+obs ~fire+hik.ln, data=utst_umf)
utsttop5.hum <- occu(~sky+obs ~fire+allhum.ln, data=utst_umf)

utsttop6.cyc <- occu(~sky+obs ~fire+elev10+cyc.ln, data=utst_umf)
utsttop6.hik <- occu(~sky+obs ~fire+elev10+hik.ln, data=utst_umf)
utsttop6.hum <- occu(~sky+obs ~fire+elev10+allhum.ln, data=utst_umf)

## adding human vars to psi and p

utsttop1.cycp<- occu(~sky+obs+hum ~mds2+fire+chili10+cyc.ln, data=utst_umf)
utsttop1.hikp<- occu(~sky+obs+hum ~mds2+fire+chili10+hik.ln, data=utst_umf)
utsttop1.hump<- occu(~sky+obs+hum ~mds2+fire+chili10+allhum.ln, data=utst_umf)
utsttop2.cycp<- occu(~sky+obs+hum ~fire+chili10+cyc.ln, data=utst_umf)
utsttop2.hikp<- occu(~sky+obs+hum ~fire+chili10+hik.ln, data=utst_umf)
utsttop2.hump<- occu(~sky+obs+hum ~fire+chili10+allhum.ln, data=utst_umf)
utsttop3.cycp<- occu(~sky+obs+hum ~mds2+fire+elev10+cyc.ln, data=utst_umf)
utsttop3.hikp<- occu(~sky+obs+hum ~mds2+fire+elev10+hik.ln, data=utst_umf)
utsttop3.hump<- occu(~sky+obs+hum ~mds2+fire+elev10+allhum.ln, data=utst_umf)
utsttop4.cycp<- occu(~sky+obs+hum ~mds2+fire+cyc.ln, data=utst_umf)
utsttop4.hikp<- occu(~sky+obs+hum ~mds2+fire+hik.ln, data=utst_umf)
utsttop4.hump<- occu(~sky+obs+hum ~mds2+fire+allhum.ln, data=utst_umf)
utsttop5.cycp<- occu(~sky+obs+hum ~fire+cyc.ln, data=utst_umf)
utsttop5.hikp<- occu(~sky+obs+hum ~fire+hik.ln, data=utst_umf)
utsttop5.hump<- occu(~sky+obs+hum ~fire+allhum.ln, data=utst_umf)
utsttop6.cycp<- occu(~sky+obs+hum ~fire+elev10+cyc.ln, data=utst_umf)
utsttop6.hikp<- occu(~sky+obs+hum ~fire+elev10+hik.ln, data=utst_umf)
utsttop6.hump<- occu(~sky+obs+hum ~fire+elev10+allhum.ln, data=utst_umf)

utst.v.models <- list("utsttop2.cyc"=utsttop2.cyc,
                      "utsttop2.hik"=utsttop2.hik,
                      "utsttop2.hum"=utsttop2.hum,
                      "utsttop2"=utsttop2,
                      "utsttop2.cycp"=utsttop2.cycp,
                      "utsttop2.hikp"=utsttop2.hikp,
                      "utsttop2.hump"=utsttop2.hump,
                      "utsttop1.cyc"=utsttop1.cyc,
                      "utsttop1.hik"=utsttop1.hik,
                      "utsttop1.hum"=utsttop1.hum,
                      "utsttop1"=utsttop1,
                      "utsttop1.cycp"=utsttop1.cycp,
                      "utsttop1.hikp"=utsttop1.hikp,
                      "utsttop1.hump"=utsttop1.hump,
                      "utsttop3.cyc"=utsttop3.cyc,
                      "utsttop3.hik"=utsttop3.hik,
                      "utsttop3.hum"=utsttop3.hum,
                      "utsttop3"=utsttop3,
                      "utsttop3.cycp"=utsttop3.cycp,
                      "utsttop3.hikp"=utsttop3.hikp,
                      "utsttop3.hump"=utsttop3.hump,
                      "utsttop4.cyc"=utsttop4.cyc,
                      "utsttop4.hik"=utsttop4.hik,
                      "utsttop4.hum"=utsttop4.hum,
                      "utsttop4"=utsttop4,
                      "utsttop4.cycp"=utsttop4.cycp,
                      "utsttop4.hikp"=utsttop4.hikp,
                      "utsttop4.hump"=utsttop4.hump,
                      "utsttop5.cyc"=utsttop5.cyc,
                      "utsttop5.hik"=utsttop5.hik,
                      "utsttop5.hum"=utsttop5.hum,
                      "utsttop5"=utsttop5,
                      "utsttop5.cycp"=utsttop5.cycp,
                      "utsttop5.hikp"=utsttop5.hikp,
                      "utsttop5.hump"=utsttop5.hump,
                      "utsttop6.cyc"=utsttop6.cyc,
                      "utsttop6.hik"=utsttop6.hik,
                      "utsttop6.hum"=utsttop6.hum,
                      "utsttop6"=utsttop6,
                      "utsttop6.cycp"=utsttop6.cycp,
                      "utsttop6.hikp"=utsttop6.hikp,
                      "utsttop6.hump"=utsttop6.hump,
                      "utst.null"=utst.null)
utst.v.aic <- aictab(utst.v.models,nobs=92,c.hat=1.74)

utsttophik <- utsttop5.hik
utsttopcyc <- utsttop5.cycp
utsttophum <- utsttop5.hum

utst.topfire <- utsttop5.hum
utst.topchili <- utsttop2.hik
utst.topmds2 <- utsttop4.hum
utst.topelev <- utsttop6.hik

confint(utsttophik, type="state")
confint(utsttopcyc, type="state")
confint(utsttophum, type="state")

#################### run occupancy models for SCOC

## test null model to get ballpark for p and psi estimates
scoc.null <- occu(~1 ~1, data=scoc_umf)
#summary(scoc.null)

## gof test
#scoc.global <- occu(~temp+I(temp^2)+jul+I(jul^2)+obsn ~mds1+mds2+fire+ndvi+elev10+chili10+hik.ln, data=scoc_umf)
#scoc.obs.boot <- mb.gof.test(scoc.global,nsim=1000) 
#scoc.obs.boot

## models w/varying detection
scocp1 <- occu(~sky ~mds1+mds2+fire, data=scoc_umf)
scocp2 <- occu(~wind ~mds1+mds2+fire, data=scoc_umf)
scocp3 <- occu(~temp+I(temp^2) ~mds1+mds2+fire, data=scoc_umf)
scocp4 <- occu(~jul+I(jul^2) ~mds1+mds2+fire, data=scoc_umf)
scocp5 <- occu(~time+I(time^2) ~mds1+mds2+fire, data=scoc_umf)
scocp6 <- occu(~obs ~mds1+mds2+fire, data=scoc_umf)
scocp7 <- occu(~trwid ~mds1+mds2+fire, data=scoc_umf)
scocp8 <- occu(~obsn ~mds1+mds2+fire, data=scoc_umf)

scocp9 <- occu(~sky+jul+I(jul^2) ~mds1+mds2+fire, data=scoc_umf)
scocp10 <- occu(~sky+time+I(time^2) ~mds1+mds2+fire, data=scoc_umf)
scocp11 <- occu(~sky+obs ~mds1+mds2+fire, data=scoc_umf)
scocp12 <- occu(~sky+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp13 <- occu(~sky+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp14 <- occu(~wind+jul+I(jul^2) ~mds1+mds2+fire, data=scoc_umf)
scocp15 <- occu(~wind+time+I(time^2) ~mds1+mds2+fire, data=scoc_umf)
scocp16 <- occu(~wind+obs ~mds1+mds2+fire, data=scoc_umf)
scocp17 <- occu(~wind+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp18 <- occu(~wind+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp19 <- occu(~temp+I(temp^2)+jul+I(jul^2) ~mds1+mds2+fire, data=scoc_umf)
scocp20 <- occu(~temp+I(temp^2)+time+I(time^2) ~mds1+mds2+fire, data=scoc_umf)
scocp21 <- occu(~temp+I(temp^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp22 <- occu(~temp+I(temp^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp23 <- occu(~temp+I(temp^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp24 <- occu(~jul+I(jul^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp25 <- occu(~jul+I(jul^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp26 <- occu(~jul+I(jul^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp27 <- occu(~time+I(time^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp28 <- occu(~time+I(time^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp29 <- occu(~time+I(time^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp30 <- occu(~sky+jul+I(jul^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp31 <- occu(~sky+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp32 <- occu(~sky+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp33 <- occu(~sky+time+I(time^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp34 <- occu(~sky+time+I(time^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp35 <- occu(~sky+time+I(time^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp36 <- occu(~wind+jul+I(jul^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp37 <- occu(~wind+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp38 <- occu(~wind+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp39 <- occu(~wind+time+I(time^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp40 <- occu(~wind+time+I(time^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp41 <- occu(~wind+time+I(time^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp42 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp43 <- occu(~temp+I(temp^2)+jul+I(jul^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp44 <- occu(~temp+I(temp^2)+jul+I(jul^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scocp45 <- occu(~temp+I(temp^2)+time+I(time^2)+obs ~mds1+mds2+fire, data=scoc_umf)
scocp46 <- occu(~temp+I(temp^2)+time+I(time^2)+trwid ~mds1+mds2+fire, data=scoc_umf)
scocp47 <- occu(~temp+I(temp^2)+time+I(time^2)+obsn ~mds1+mds2+fire, data=scoc_umf)

scoc.p.models <- list(scocp1,scocp2,scocp3,scocp4,scocp5,scocp6,scocp7,scocp8,scocp9,scocp10,scocp11,scocp12,scocp13,scocp14,scocp15,scocp16,scocp17,scocp18,scocp19,scocp20,scocp21,scocp22,scocp23,scocp24,scocp25,scocp26,scocp27,scocp28,scocp29,scocp30,scocp31,scocp32,scocp33,scocp34,scocp35,scocp36,scocp37,scocp38,scocp39,scocp40,scocp41,scocp42,scocp43,scocp44,scocp45,scocp46,scocp47)

scoc.p.aic <- aictab(scoc.p.models,nobs=92,c.hat=1.92)

# varying psi
scoct1 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1, data=scoc_umf)
scoct2 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+mds2, data=scoc_umf)
scoct3 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+fire, data=scoc_umf)
scoct4 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+ndvi, data=scoc_umf)
scoct5 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2, data=scoc_umf)
scoct6 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+fire, data=scoc_umf)
scoct7 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+ndvi, data=scoc_umf)
scoct8 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire, data=scoc_umf)
scoct9 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire+ndvi, data=scoc_umf)
scoct10 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ndvi, data=scoc_umf)

scoct11 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+elev10, data=scoc_umf)
scoct12 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+mds2+elev10, data=scoc_umf)
scoct13 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+fire+elev10, data=scoc_umf)
scoct14 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+ndvi+elev10, data=scoc_umf)
scoct15 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+elev10, data=scoc_umf)
scoct16 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+fire+elev10, data=scoc_umf)
scoct17 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+ndvi+elev10, data=scoc_umf)
scoct18 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire+elev10, data=scoc_umf)
scoct19 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire+ndvi+elev10, data=scoc_umf)
scoct20 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ndvi+elev10, data=scoc_umf)

scoct21 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+chili10, data=scoc_umf)
scoct22 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+mds2+chili10, data=scoc_umf)
scoct23 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+fire+chili10, data=scoc_umf)
scoct24 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds1+ndvi+chili10, data=scoc_umf)
scoct25 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+chili10, data=scoc_umf)
scoct26 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+fire+chili10, data=scoc_umf)
scoct27 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~mds2+ndvi+chili10, data=scoc_umf)
scoct28 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire+chili10, data=scoc_umf)
scoct29 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~fire+ndvi+chili10, data=scoc_umf)
scoct30 <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ndvi+chili10, data=scoc_umf)

scoct1b <- occu(~wind + jul + I(jul^2) ~mds1, data=scoc_umf)
scoct2b <- occu(~wind + jul + I(jul^2) ~mds1+mds2, data=scoc_umf)
scoct3b <- occu(~wind + jul + I(jul^2) ~mds1+fire, data=scoc_umf)
scoct4b <- occu(~wind + jul + I(jul^2) ~mds1+ndvi, data=scoc_umf)
scoct5b <- occu(~wind + jul + I(jul^2) ~mds2, data=scoc_umf)
scoct6b <- occu(~wind + jul + I(jul^2) ~mds2+fire, data=scoc_umf)
scoct7b <- occu(~wind + jul + I(jul^2) ~mds2+ndvi, data=scoc_umf)
scoct8b <- occu(~wind + jul + I(jul^2) ~fire, data=scoc_umf)
scoct9b <- occu(~wind + jul + I(jul^2) ~fire+ndvi, data=scoc_umf)
scoct10b <- occu(~wind + jul + I(jul^2) ~ndvi, data=scoc_umf)

scoct11b <- occu(~wind + jul + I(jul^2) ~mds1+elev10, data=scoc_umf)
scoct12b <- occu(~wind + jul + I(jul^2) ~mds1+mds2+elev10, data=scoc_umf)
scoct13b <- occu(~wind + jul + I(jul^2) ~mds1+fire+elev10, data=scoc_umf)
scoct14b <- occu(~wind + jul + I(jul^2) ~mds1+ndvi+elev10, data=scoc_umf)
scoct15b <- occu(~wind + jul + I(jul^2) ~mds2+elev10, data=scoc_umf)
scoct16b <- occu(~wind + jul + I(jul^2) ~mds2+fire+elev10, data=scoc_umf)
scoct17b <- occu(~wind + jul + I(jul^2) ~mds2+ndvi+elev10, data=scoc_umf)
scoct18b <- occu(~wind + jul + I(jul^2) ~fire+elev10, data=scoc_umf)
scoct19b <- occu(~wind + jul + I(jul^2) ~fire+ndvi+elev10, data=scoc_umf)
scoct20b <- occu(~wind + jul + I(jul^2) ~ndvi+elev10, data=scoc_umf)

scoct21b <- occu(~wind + jul + I(jul^2) ~mds1+chili10, data=scoc_umf)
scoct22b <- occu(~wind + jul + I(jul^2) ~mds1+mds2+chili10, data=scoc_umf)
scoct23b <- occu(~wind + jul + I(jul^2) ~mds1+fire+chili10, data=scoc_umf)
scoct24b <- occu(~wind + jul + I(jul^2) ~mds1+ndvi+chili10, data=scoc_umf)
scoct25b <- occu(~wind + jul + I(jul^2) ~mds2+chili10, data=scoc_umf)
scoct26b <- occu(~wind + jul + I(jul^2) ~mds2+fire+chili10, data=scoc_umf)
scoct27b <- occu(~wind + jul + I(jul^2) ~mds2+ndvi+chili10, data=scoc_umf)
scoct28b <- occu(~wind + jul + I(jul^2) ~fire+chili10, data=scoc_umf)
scoct29b <- occu(~wind + jul + I(jul^2) ~fire+ndvi+chili10, data=scoc_umf)
scoct30b <- occu(~wind + jul + I(jul^2) ~ndvi+chili10, data=scoc_umf)

scoc.t.models <- list(scoct1,scoct2,scoct3,scoct4,scoct5,scoct6,scoct7,scoct8,scoct9,scoct10,scoct11,scoct12,scoct13,scoct14,scoct15,scoct16,scoct17,scoct18,scoct19,scoct20,scoct21,scoct22,scoct23,scoct24,scoct25,scoct26,scoct27,scoct28,scoct29,scoct30,scoct1b,scoct2b,scoct3b,scoct4b,scoct5b,scoct6b,scoct7b,scoct8b,scoct9b,scoct10b,scoct11b,scoct12b,scoct13b,scoct14b,scoct15b,scoct16b,scoct17b,scoct18b,scoct19b,scoct20b,scoct21b,scoct22b,scoct23b,scoct24b,scoct25b,scoct26b,scoct27b,scoct28b,scoct29b,scoct30b)
scoc.t.aic <- aictab(scoc.t.models,nobs=92,c.hat=1.92)

scoctop1 <- scoct10
scoctop2 <- scoct20
scoctop3 <- scoct20b
scoctop4 <- scoct10b
scoctop5 <- scoct4

## adding human vars to top model(s)
# to psi
scoctop1.cyc <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi+cyc.ln, data=scoc_umf)
scoctop1.hik <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi+hik.ln, data=scoc_umf)
scoctop1.hum <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi+allhum.ln, data=scoc_umf)

scoctop2.cyc <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi + elev10+cyc.ln, data=scoc_umf)
scoctop2.hik <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi + elev10+hik.ln, data=scoc_umf)
scoctop2.hum <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ ndvi + elev10+allhum.ln, data=scoc_umf)

scoctop3.cyc <- occu(~wind + jul + I(jul^2) ~ ndvi + elev10+cyc.ln, data=scoc_umf)
scoctop3.hik <- occu(~wind + jul + I(jul^2) ~ ndvi + elev10+hik.ln, data=scoc_umf)
scoctop3.hum <- occu(~wind + jul + I(jul^2) ~ ndvi + elev10+allhum.ln, data=scoc_umf)

scoctop4.cyc <- occu(~wind + jul + I(jul^2) ~ ndvi +cyc.ln, data=scoc_umf)
scoctop4.hik <- occu(~wind + jul + I(jul^2) ~ ndvi +hik.ln, data=scoc_umf)
scoctop4.hum <- occu(~wind + jul + I(jul^2) ~ ndvi +allhum.ln, data=scoc_umf)

scoctop5.cyc <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ mds1 + ndvi +cyc.ln, data=scoc_umf)
scoctop5.hik <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ mds1 + ndvi +hik.ln, data=scoc_umf)
scoctop5.hum <- occu(~temp + I(temp^2) + jul + I(jul^2) ~ mds1 + ndvi +allhum.ln, data=scoc_umf)

# to both
scoctop1.cycp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi+cyc.ln, data=scoc_umf)
scoctop1.hikp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi+hik.ln, data=scoc_umf)
scoctop1.hump <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi+allhum.ln, data=scoc_umf)

scoctop2.cycp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi + elev10+cyc.ln, data=scoc_umf)
scoctop2.hikp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi + elev10+hik.ln, data=scoc_umf)
scoctop2.hump <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ ndvi + elev10+allhum.ln, data=scoc_umf)

scoctop3.cycp <- occu(~wind + jul + I(jul^2)+hum ~ ndvi + elev10+cyc.ln, data=scoc_umf)
scoctop3.hikp <- occu(~wind + jul + I(jul^2)+hum ~ ndvi + elev10+hik.ln, data=scoc_umf)
scoctop3.hump <- occu(~wind + jul + I(jul^2)+hum ~ ndvi + elev10+allhum.ln, data=scoc_umf)

scoctop4.cycp <- occu(~wind + jul + I(jul^2)+hum ~ ndvi +cyc.ln, data=scoc_umf)
scoctop4.hikp <- occu(~wind + jul + I(jul^2)+hum ~ ndvi +hik.ln, data=scoc_umf)
scoctop4.hump <- occu(~wind + jul + I(jul^2)+hum ~ ndvi +allhum.ln, data=scoc_umf)

scoctop5.cycp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ mds1 + ndvi +cyc.ln, data=scoc_umf)
scoctop5.hikp <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ mds1 + ndvi +hik.ln, data=scoc_umf)
scoctop5.hump <- occu(~temp + I(temp^2) + jul + I(jul^2)+hum ~ mds1 + ndvi +allhum.ln, data=scoc_umf)

scoc.v.models <- list("scoctop1.cyc"=scoctop1.cyc,
                      "scoctop1.hik"=scoctop1.hik,
                      "scoctop1.hum"=scoctop1.hum,
                      "scoctop1"=scoctop1,
                      "scoctop2.cyc"=scoctop2.cyc,
                      "scoctop2.hik"=scoctop2.hik,
                      "scoctop2.hum"=scoctop2.hum,
                      "scoctop2"=scoctop2,
                      "scoctop3.cyc"=scoctop3.cyc,
                      "scoctop3.hik"=scoctop3.hik,
                      "scoctop3.hum"=scoctop3.hum,
                      "scoctop3"=scoctop3,
                      "scoctop4.cyc"=scoctop4.cyc,
                      "scoctop4.hik"=scoctop4.hik,
                      "scoctop4.hum"=scoctop4.hum,
                      "scoctop4"=scoctop4,
                      "scoctop5.cyc"=scoctop5.cyc,
                      "scoctop5.hik"=scoctop5.hik,
                      "scoctop5.hum"=scoctop5.hum,
                      "scoctop5"=scoctop5,
                      "scoctop1.cycp"=scoctop1.cycp,
                      "scoctop1.hikp"=scoctop1.hikp,
                      "scoctop1.hump"=scoctop1.hump,
                      "scoctop2.cycp"=scoctop2.cycp,
                      "scoctop2.hikp"=scoctop2.hikp,
                      "scoctop2.hump"=scoctop2.hump,
                      "scoctop3.cycp"=scoctop3.cycp,
                      "scoctop3.hikp"=scoctop3.hikp,
                      "scoctop3.hump"=scoctop3.hump,
                      "scoctop4.cycp"=scoctop4.cycp,
                      "scoctop4.hikp"=scoctop4.hikp,
                      "scoctop4.hump"=scoctop4.hump,
                      "scoctop5.cycp"=scoctop5.cycp,
                      "scoctop5.hikp"=scoctop5.hikp,
                      "scoctop5.hump"=scoctop5.hump,
                      "scoc.null"=scoc.null)
scoc.v.aic <- aictab(scoc.v.models,nobs=92,c.hat=1.92)

scoctophum <- scoctop1.hum
scoctophik <- scoctop1.hikp
scoctopcyc <- scoctop2.cycp

scoc.topndvi <- scoctop1.hum
scoc.topmds1 <- scoctop5.cycp
scoc.topelev <- scoctop2.cycp

confint(scoctophum, type="state")
confint(scoctophik, type="state")
confint(scoctopcyc, type="state")
