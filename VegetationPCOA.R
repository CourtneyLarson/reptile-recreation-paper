### Title: Vegetation PCOA
### Author: Courtney L. Larson
### Description: Vegetation analysis for Larson et al. 2023. Reptile responses to outdoor recreation in urban habitat fragments. Urban Ecosystems.

# packages
#install.packages("dplyr","vegan")
library(dplyr)
library(vegan)

# read in data
veg <- read.csv("C:/Courtney/PhD_work/HerpAnalysis/UrbanEcol submission/DataCode_forRepository/VegSurveys_compiled.csv",header=TRUE)

# sum up the number of times each species detected at each point (summing across hits)
veg.sum <- veg[,c(1,3:66)] %>%
  group_by(Survey) %>%
  summarize_all(funs(sum),na.rm=TRUE) %>%
  mutate(meas = 40-missing)

# use the sums to calculate % cover of each species at each point
veg.prop <- veg.sum %>%
  mutate_if(is.numeric,funs(./meas))
veg.prop <- within(veg.prop, rm(hit,missing))
rowSums(veg.prop[,2:63]) #this is the total proportion cover value, all over 1 because of taking multiple hits at each measurement location
veg.chars <- dplyr::select(veg.chars, -meas)
veg.chars <- as.data.frame(veg.chars)

# drop species with very few observations -  < 5 sites
veg.prev <- ifelse(veg.chars[,2:65] > 0, 1, 0)
veg.prev <- rbind(veg.prev,colSums(veg.prev))
veg.morethan5 <- as.data.frame(veg.prev[,which(veg.prev[93,] > 5)])
veg.chars3 <- dplyr::select(veg.chars, one_of(names(veg.morethan5)))

## run PCOA
pcoa <- capscale(veg.chars3~1,distance="bray")
summary(pcoa)

screeplot(pcoa)

# look at results in order by MDS axes to help interpret
top3 <- as.data.frame(pcoa$CA$v[,1:3])
top3$sp <- row.names(top3)
arrange(top3,MDS1)
arrange(top3,MDS2)
arrange(top3,MDS3)

pcoaSite <- as.data.frame(pcoa$CA$u[,1:2])
pcoaSite$PointGrp <- veg.chars$Survey

# write out results
write.csv(pcoaSite, "MDS.csv")


