### Title: Riptile richness and trait analysis
### Author: Courtney L. Larson
### Description: Reptile species richness and community trait analyses for Larson et al. 2023. Reptile responses to outdoor recreation in urban habitat fragments.

#################### packages
#install.packages("dplyr", "vegan")
library(dplyr)
library(vegan)


#################### read in data
setwd("C:/Courtney/PhD_work/HerpAnalysis/UrbanEcol submission/DataCode_forRepository")
herp <- read.csv("SiteXSpecies.csv",header=TRUE)
siteCovs <- read.csv("SiteCovs.csv",header=TRUE)

#################### data setup for sp richness & comm comp

siteCovs2 <- siteCovs %>%
  group_by(PointGrp) %>%
  summarise_if(is.numeric,mean,na.rm=TRUE)

#################### total species richness models

## estimate total sp richness
rich.tot <- specpool(herp[,3:24])

## estimate richness by site
siteS <- estimateR(herp[,3:24])
Chao.sp <- data.frame("ChaoRich"=siteS[2,],"ChaoSE"=siteS[3,],"PointGrp"=herp$Point,"Sp.obs"=siteS[1,])

Chao.sp <- full_join(Chao.sp,siteCovs2,by="PointGrp")
Chao.sp$ChaoRich <- ifelse(is.na(Chao.sp$ChaoRich),0,Chao.sp$ChaoRich)
Chao.sp$ChaoSE <- ifelse(is.na(Chao.sp$ChaoSE),0,Chao.sp$ChaoSE)
Chao.sp$Sp.obs <- ifelse(is.na(Chao.sp$Sp.obs),0,Chao.sp$Sp.obs)

# scale numeric covariates
Chao.sp[,5:26] <- lapply(Chao.sp[,5:26], scale)
Chao.sp[,5:26] <- sapply(Chao.sp[,5:26], as.numeric)

## linear regressions using Chao sp as response
# fit models with up to 2 habitat vars and 1 (or no) topo var
r1 <- lm(ChaoRich~mds1+mds2, data=Chao.sp)
r2 <- lm(ChaoRich~mds1+fire, data=Chao.sp)
r3 <- lm(ChaoRich~mds1+ndvi, data=Chao.sp)
r4 <- lm(ChaoRich~mds2+fire, data=Chao.sp)
r5 <- lm(ChaoRich~mds2+ndvi, data=Chao.sp)
r6 <- lm(ChaoRich~fire+ndvi, data=Chao.sp)

r7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.sp)
r8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.sp)
r9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.sp)
r10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.sp)
r11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.sp)
r12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.sp)

r13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.sp)
r14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.sp)
r15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.sp)
r16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.sp)
r17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.sp)
r18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.sp)

r19 <- lm(ChaoRich~mds1, data=Chao.sp)
r20 <- lm(ChaoRich~mds2, data=Chao.sp)
r21 <- lm(ChaoRich~fire, data=Chao.sp)
r22 <- lm(ChaoRich~ndvi, data=Chao.sp)
r23 <- lm(ChaoRich~elev10, data=Chao.sp)
r24 <- lm(ChaoRich~chili10, data=Chao.sp)

env.list <- list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24)
env.aic <- aictab(env.list)

rich1 <- r15
rich2 <- r6
rich3 <- r18
rich4 <- r22

## adding human vars to best models
rich3.cyc <- lm(ChaoRich~fire+ndvi+elev10+cyc.ln, data=Chao.sp)
rich3.hik <- lm(ChaoRich~fire+ndvi+elev10+hik.ln, data=Chao.sp)
rich3.hum <- lm(ChaoRich~fire+ndvi+elev10+allhum.ln, data=Chao.sp)

rich2.cyc <- lm(ChaoRich~fire+ndvi+cyc.ln, data=Chao.sp)
rich2.hik <- lm(ChaoRich~fire+ndvi+hik.ln, data=Chao.sp)
rich2.hum <- lm(ChaoRich~fire+ndvi+allhum.ln, data=Chao.sp)

rich1.cyc <- lm(ChaoRich~mds1+ndvi+elev10+cyc.ln, data=Chao.sp)
rich1.hik <- lm(ChaoRich~mds1+ndvi+elev10+hik.ln, data=Chao.sp)
rich1.hum <- lm(ChaoRich~mds1+ndvi+elev10+allhum.ln, data=Chao.sp)

rich4.cyc <- lm(ChaoRich~ndvi+cyc.ln, data=Chao.sp)
rich4.hik <- lm(ChaoRich~ndvi+hik.ln, data=Chao.sp)
rich4.hum <- lm(ChaoRich~ndvi+allhum.ln, data=Chao.sp)

richlist <- list("fire+NDVI+elevation+**cyclist**"=rich3.cyc,
                 "fire+NDVI+elevation+**pedestrian**"=rich3.hik,
                 "fire+NDVI+elevation+**total human**"=rich3.hum,
                 "fire+NDVI+elevation"=rich3,
                 "fire+NDVI+**cyclist**"=rich2.cyc,
                 "fire+NDVI+**pedestrian**"=rich2.hik,
                 "fire+NDVI+**total human**"=rich2.hum,
                 "fire+NDVI"=rich2,
                 "nMDS1+NDVI+elevation+**cyclist**"=rich1.cyc,
                 "nMDS1+NDVI+elevation+**pedestrian**"=rich1.hik,
                 "nMDS1+NDVI+elevation+**total human**"=rich1.hum,
                 "nMDS1+NDVI+elevation"=rich1,
                 "NDVI+**cyclist**"=rich4.cyc,
                 "NDVI+**pedestrian**"=rich4.hik,
                 "NDVI+**total human**"=rich4.hum,
                 "NDVI$"=rich4)
rich.aic <- aictab(richlist)

rich.tophum <- rich2.hum
rich.tophik <- rich2.hik
rich.topcyc <- rich2.cyc

rich.topfire <- rich2.hum
rich.topndvi <- rich2.hum
rich.topmds1 <- rich1
rich.topelev <- rich1

# get 95% confidence interval for best human model
confint(rich.tophum)
confint(rich.tophik)
confint(rich.topcyc)

#################### lizard richness

liz <- dplyr::select(herp,Point,DailyMeanHuman,ANST,ASHY,ASTI,COVA,ELMU,PHBL,PLGI,PLSK,SCOC,SCOR,UTST)

rich.liz <- specpool(liz[,3:13])
siteliz <- estimateR(liz[,3:13])
Chao.liz <- data.frame("ChaoRich"=siteliz[2,],"ChaoSE"=siteliz[3,],"PointGrp"=liz$Point,"Sp.obs"=siteliz[1,])

Chao.liz <- full_join(Chao.liz,siteCovs2,by="PointGrp")
Chao.liz$ChaoRich <- ifelse(is.na(Chao.liz$ChaoRich),0,Chao.liz$ChaoRich)
Chao.liz$ChaoSE <- ifelse(is.na(Chao.liz$ChaoSE),0,Chao.liz$ChaoSE)
Chao.liz$Sp.obs <- ifelse(is.na(Chao.liz$Sp.obs),0,Chao.liz$Sp.obs)

# scale numeric covariates
Chao.liz[,5:26] <- lapply(Chao.liz[,5:26], scale)
Chao.liz[,5:26] <- sapply(Chao.liz[,5:26], as.numeric)

# fit linear models
liz1 <- lm(ChaoRich~mds1+mds2, data=Chao.liz)
liz2 <- lm(ChaoRich~mds1+fire, data=Chao.liz)
liz3 <- lm(ChaoRich~mds1+ndvi, data=Chao.liz)
liz4 <- lm(ChaoRich~mds2+fire, data=Chao.liz)
liz5 <- lm(ChaoRich~mds2+ndvi, data=Chao.liz)
liz6 <- lm(ChaoRich~fire+ndvi, data=Chao.liz)

liz7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.liz)
liz8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.liz)
liz9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.liz)
liz10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.liz)
liz11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.liz)
liz12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.liz)

liz13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.liz)
liz14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.liz)
liz15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.liz)
liz16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.liz)
liz17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.liz)
liz18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.liz)

liz19 <- lm(ChaoRich~mds1, data=Chao.liz)
liz20 <- lm(ChaoRich~mds2, data=Chao.liz)
liz21 <- lm(ChaoRich~fire, data=Chao.liz)
liz22 <- lm(ChaoRich~ndvi, data=Chao.liz)
liz23 <- lm(ChaoRich~elev10, data=Chao.liz)
liz24 <- lm(ChaoRich~chili10, data=Chao.liz)

liz.list1 <- list(liz1,liz2,liz3,liz4,liz5,liz6,liz7,liz8,liz9,liz10,liz11,liz12,liz13,liz14,liz15,liz16,liz17,liz18,liz19,liz20,liz21,liz22,liz23,liz24)
liz.aic <- aictab(liz.list1)

liztop1 <- liz15
liztop2 <- liz6
liztop3 <- liz17
liztop4 <- liz18
liztop5 <- liz22

liztop2.hum <- lm(ChaoRich~fire+ndvi+allhum.ln, data=Chao.liz)
liztop2.cyc <- lm(ChaoRich~fire+ndvi+cyc.ln, data=Chao.liz)
liztop2.hik <- lm(ChaoRich~fire+ndvi+hik.ln, data=Chao.liz)

liztop1.hum <- lm(ChaoRich~ndvi+allhum.ln, data=Chao.liz)
liztop1.cyc <- lm(ChaoRich~ndvi+cyc.ln, data=Chao.liz)
liztop1.hik <- lm(ChaoRich~ndvi+hik.ln, data=Chao.liz)

liztop4.hum <- lm(ChaoRich~fire+ndvi+elev10+allhum.ln, data=Chao.liz)
liztop4.cyc <- lm(ChaoRich~fire+ndvi+elev10+cyc.ln, data=Chao.liz)
liztop4.hik <- lm(ChaoRich~fire+ndvi+elev10+hik.ln, data=Chao.liz)

liztop5.hum <- lm(ChaoRich~mds1+ndvi+elev10+allhum.ln, data=Chao.liz)
liztop5.cyc <- lm(ChaoRich~mds1+ndvi+elev10+cyc.ln, data=Chao.liz)
liztop5.hik <- lm(ChaoRich~mds1+ndvi+elev10+hik.ln, data=Chao.liz)

liztop3.hum <- lm(ChaoRich~mds2+ndvi+elev10+allhum.ln, data=Chao.liz)
liztop3.cyc <- lm(ChaoRich~mds2+ndvi+elev10+cyc.ln, data=Chao.liz)
liztop3.hik <- lm(ChaoRich~mds2+ndvi+elev10+hik.ln, data=Chao.liz)

liz.list2 <- list("fire + NDVI + **total human**"=liztop2.hum,
                  "fire + NDVI + **cyclist**"=liztop2.cyc,
                  "fire + NDVI + **pedestrian**"=liztop2.hik,
                  "fire + NDVI"=liztop2,
                  "NDVI + **total human**"=liztop1.hum,
                  "NDVI + **cyclist**"=liztop1.cyc,
                  "NDVI + **pedestrian**"=liztop1.hik,
                  "NDVI**"=liztop1,
                  "fire + NDVI + elevation + **total human**"=liztop4.hum,
                  "fire + NDVI + elevation + **cyclist**"=liztop4.cyc,
                  "fire + NDVI + elevation + **pedestrian**"=liztop4.hik,
                  "fire + NDVI + elevation"=liztop4,
                  "nMDS1 + NDVI + elevation + **total human**"=liztop5.hum,
                  "nMDS1 + NDVI + elevation + **cyclist**"=liztop5.cyc,
                  "nMDS1 + NDVI + elevation + **pedestrian**"=liztop5.hik,
                  "nMDS1 + NDVI + elevation"=liztop5,
                  "nMDS2 + NDVI + elevation + **total human**"=liztop3.hum,
                  "nMDS2 + NDVI + elevation + **cyclist**"=liztop3.cyc,
                  "nMDS2 + NDVI + elevation + **pedestrian**"=liztop3.hik,
                  "nMDS2 + NDVI + elevation"=liztop3)
liz.aic.hum <- aictab(liz.list2)

liz.tophum <- liztop1.hum
liz.tophik <- liztop1.hik
liz.topcyc <- liztop2.cyc

liz.topfire <- liztop2.hum
liz.topndvi <- liztop1.hum
liz.topmds1 <- liztop5.hum
liz.topelev <- liztop4.hum

# get 95% conf int for best human model
confint(liztop1.hum)

#################### snake richness

snek <- dplyr::select(herp,Point,DailyMeanHuman,COFL,COLA,CROR,CRRU,DIPU,HYOC,LACA,PICA,RHLE,SAHE,THHA)

rich.snek <- specpool(snek[,3:13])
sitesnek <- estimateR(snek[,3:13])
Chao.snek <- data.frame("ChaoRich"=sitesnek[2,],"ChaoSE"=sitesnek[3,],"PointGrp"=snek$Point,"Sp.obs"=sitesnek[1,])

Chao.snek <- full_join(Chao.snek,siteCovs2,by="PointGrp")
Chao.snek$ChaoRich <- ifelse(is.na(Chao.snek$ChaoRich),0,Chao.snek$ChaoRich)
Chao.snek$ChaoSE <- ifelse(is.na(Chao.snek$ChaoSE),0,Chao.snek$ChaoSE)
Chao.snek$Sp.obs <- ifelse(is.na(Chao.snek$Sp.obs),0,Chao.snek$Sp.obs)

# scale numeric covariates
Chao.snek[,5:26] <- lapply(Chao.snek[,5:26], scale)
Chao.snek[,5:26] <- sapply(Chao.snek[,5:26], as.numeric)

# fit linear models
snek1 <- lm(ChaoRich~mds1+mds2, data=Chao.snek)
snek2 <- lm(ChaoRich~mds1+fire, data=Chao.snek)
snek3 <- lm(ChaoRich~mds1+ndvi, data=Chao.snek)
snek4 <- lm(ChaoRich~mds2+fire, data=Chao.snek)
snek5 <- lm(ChaoRich~mds2+ndvi, data=Chao.snek)
snek6 <- lm(ChaoRich~fire+ndvi, data=Chao.snek)

snek7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.snek)
snek8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.snek)
snek9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.snek)
snek10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.snek)
snek11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.snek)
snek12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.snek)

snek13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.snek)
snek14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.snek)
snek15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.snek)
snek16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.snek)
snek17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.snek)
snek18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.snek)

snek19 <- lm(ChaoRich~mds1, data=Chao.snek)
snek20 <- lm(ChaoRich~mds2, data=Chao.snek)
snek21 <- lm(ChaoRich~fire, data=Chao.snek)
snek22 <- lm(ChaoRich~ndvi, data=Chao.snek)
snek23 <- lm(ChaoRich~elev10, data=Chao.snek)
snek24 <- lm(ChaoRich~chili10, data=Chao.snek)

snek.list1 <- list(snek1,snek2,snek3,snek4,snek5,snek6,snek7,snek8,snek9,snek10,snek11,snek12,snek13,snek14,snek15,snek16,snek17,snek18,snek19,snek20,snek21,snek22,snek23,snek24)
snek.aic <- aictab(snek.list1)

snektop1 <- snek22
snektop2 <- snek19
snektop3 <- snek3
snektop4 <- snek5

snektop1.hum <- lm(ChaoRich~fire+allhum.ln, data=Chao.snek)
snektop4.hum <- lm(ChaoRich~mds2+ndvi+allhum.ln, data=Chao.snek)
snektop2.hum <- lm(ChaoRich~mds1+allhum.ln, data=Chao.snek)
snektop3.hum <- lm(ChaoRich~mds1+ndvi+allhum.ln, data=Chao.snek)

snek.list2 <- list(snektop1.hum,snektop2.hum,snektop3.hum,snektop4.hum,
                   snektop1,snektop2,snektop3,snektop4)
snek.aic.hum <- aictab(snek.list2)

snek.tophum <- snektop2.hum

# get 95% conf int for best human model
confint(snektop2.hum)

#################### richness of primarily transect-detected species
vt <- dplyr::select(herp,Point,DailyMeanHuman,ASHY,ASTI,COFL,COLA,CROR,CRRU,LACA,PHBL,SAHE,SCOC,SCOR,THHA,UTST)

rich.vt <- specpool(vt[,3:15])
sitevt <- estimateR(vt[,3:15])
Chao.vt <- data.frame("ChaoRich"=sitevt[2,],"ChaoSE"=sitevt[3,],"PointGrp"=vt$Point,"Sp.obs"=sitevt[1,])

Chao.vt <- full_join(Chao.vt,siteCovs2,by="PointGrp")
Chao.vt$ChaoRich <- ifelse(is.na(Chao.vt$ChaoRich),0,Chao.vt$ChaoRich)
Chao.vt$ChaoSE <- ifelse(is.na(Chao.vt$ChaoSE),0,Chao.vt$ChaoSE)
Chao.vt$Sp.obs <- ifelse(is.na(Chao.vt$Sp.obs),0,Chao.vt$Sp.obs)

# scale numeric covariates
Chao.vt[,5:26] <- lapply(Chao.vt[,5:26], scale)
Chao.vt[,5:26] <- sapply(Chao.vt[,5:26], as.numeric)

# fit linear models
vt1 <- lm(ChaoRich~mds1+mds2, data=Chao.vt)
vt2 <- lm(ChaoRich~mds1+fire, data=Chao.vt)
vt3 <- lm(ChaoRich~mds1+ndvi, data=Chao.vt)
vt4 <- lm(ChaoRich~mds2+fire, data=Chao.vt)
vt5 <- lm(ChaoRich~mds2+ndvi, data=Chao.vt)
vt6 <- lm(ChaoRich~fire+ndvi, data=Chao.vt)

vt7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.vt)
vt8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.vt)
vt9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.vt)
vt10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.vt)
vt11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.vt)
vt12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.vt)

vt13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.vt)
vt14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.vt)
vt15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.vt)
vt16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.vt)
vt17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.vt)
vt18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.vt)

vt19 <- lm(ChaoRich~mds1, data=Chao.vt)
vt20 <- lm(ChaoRich~mds2, data=Chao.vt)
vt21 <- lm(ChaoRich~fire, data=Chao.vt)
vt22 <- lm(ChaoRich~ndvi, data=Chao.vt)
vt23 <- lm(ChaoRich~elev10, data=Chao.vt)
vt24 <- lm(ChaoRich~chili10, data=Chao.vt)

vt.list1 <- list(vt1,vt2,vt3,vt4,vt5,vt6,vt7,vt8,vt9,vt10,vt11,vt12,vt13,vt14,vt15,vt16,vt17,vt18,vt19,vt20,vt21,vt22,vt23,vt24)
vt.aic <- aictab(vt.list1)

vttop1 <- vt23
vttop2 <- vt21

vttop1.hum <- lm(ChaoRich~elev10+allhum.ln, data=Chao.vt)
vttop2.hum <- lm(ChaoRich~fire+allhum.ln, data=Chao.vt)

vt.list2 <- list(vttop1.hum,vttop2.hum,
                 vttop1,vttop2)
vt.aic.hum <- aictab(vt.list2)

vt.tophum <- vttop1.hum

# get 95% conf int for best human model
confint(vttop1.hum)

#################### richness of primarily board-detected species
bd <- dplyr::select(herp,Point,DailyMeanHuman,ANST,COVA,DIPU,ELMU,HYOC,PICA,PLGI,PLSK,RHLE)

rich.bd <- specpool(bd[,3:11])
sitebd <- estimateR(bd[,3:11])
Chao.bd <- data.frame("ChaoRich"=sitebd[2,],"ChaoSE"=sitebd[3,],"PointGrp"=bd$Point,"Sp.obs"=sitebd[1,])

Chao.bd <- full_join(Chao.bd,siteCovs2,by="PointGrp")
Chao.bd$ChaoRich <- ifelse(is.na(Chao.bd$ChaoRich),0,Chao.bd$ChaoRich)
Chao.bd$ChaoSE <- ifelse(is.na(Chao.bd$ChaoSE),0,Chao.bd$ChaoSE)
Chao.bd$Sp.obs <- ifelse(is.na(Chao.bd$Sp.obs),0,Chao.bd$Sp.obs)

# scale numeric covariates
Chao.bd[,5:26] <- lapply(Chao.bd[,5:26], scale)
Chao.bd[,5:26] <- sapply(Chao.bd[,5:26], as.numeric)

# fit linear models
bd1 <- lm(ChaoRich~mds1+mds2, data=Chao.bd)
bd2 <- lm(ChaoRich~mds1+fire, data=Chao.bd)
bd3 <- lm(ChaoRich~mds1+ndvi, data=Chao.bd)
bd4 <- lm(ChaoRich~mds2+fire, data=Chao.bd)
bd5 <- lm(ChaoRich~mds2+ndvi, data=Chao.bd)
bd6 <- lm(ChaoRich~fire+ndvi, data=Chao.bd)

bd7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.bd)
bd8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.bd)
bd9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.bd)
bd10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.bd)
bd11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.bd)
bd12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.bd)

bd13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.bd)
bd14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.bd)
bd15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.bd)
bd16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.bd)
bd17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.bd)
bd18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.bd)

bd19 <- lm(ChaoRich~mds1, data=Chao.bd)
bd20 <- lm(ChaoRich~mds2, data=Chao.bd)
bd21 <- lm(ChaoRich~fire, data=Chao.bd)
bd22 <- lm(ChaoRich~ndvi, data=Chao.bd)
bd23 <- lm(ChaoRich~elev10, data=Chao.bd)
bd24 <- lm(ChaoRich~chili10, data=Chao.bd)

bd.list1 <- list(bd1,bd2,bd3,bd4,bd5,bd6,bd7,bd8,bd9,bd10,bd11,bd12,bd13,bd14,bd15,bd16,bd17,bd18,bd19,bd20,bd21,bd22,bd23,bd24)
bd.aic <- aictab(bd.list1)

bdtop1 <- bd22
bdtop2 <- bd3
bdtop3 <- bd15
bdtop4 <- bd6
bdtop5 <- bd5

bdtop1 <- bd3
bdtop2 <- bd15
bdtop3 <- bd22
bdtop4 <- bd6

bdtop3.hum <- lm(ChaoRich~ndvi+allhum.ln, data=Chao.bd)
bdtop1.hum <- lm(ChaoRich~mds1+ndvi+allhum.ln, data=Chao.bd)
bdtop2.hum <- lm(ChaoRich~mds1+ndvi+elev10+allhum.ln, data=Chao.bd)
bdtop4.hum <- lm(ChaoRich~fire+ndvi+allhum.ln, data=Chao.bd)

bd.list2 <- list(bdtop1.hum,bdtop2.hum,bdtop3.hum,bdtop4.hum,
                 bdtop1,bdtop2,bdtop3,bdtop4)
bd.aic.hum <- aictab(bd.list2)

bd.tophum <- bdtop1.hum

# get 95% conf int for best human model
confint(bdtop1.hum)

#################### specialist species richness
spec <- dplyr::select(herp,Point,DailyMeanHuman,COVA,ASHY,PHBL,PLGI,SCOR,ANST,THHA,SAHE,COFL,CRRU)

rich.spec <- specpool(spec[,3:12])
sitespec <- estimateR(spec[,3:12])
Chao.spec <- data.frame("ChaoRich"=sitespec[2,],"ChaoSE"=sitespec[3,],"PointGrp"=spec$Point,"Sp.obs"=sitespec[1,])

Chao.spec <- full_join(Chao.spec,siteCovs2,by="PointGrp")
Chao.spec$ChaoRich <- ifelse(is.na(Chao.spec$ChaoRich),0,Chao.spec$ChaoRich)
Chao.spec$ChaoSE <- ifelse(is.na(Chao.spec$ChaoSE),0,Chao.spec$ChaoSE)
Chao.spec$Sp.obs <- ifelse(is.na(Chao.spec$Sp.obs),0,Chao.spec$Sp.obs)

# scale numeric covariates
Chao.spec[,5:26] <- lapply(Chao.spec[,5:26], scale)
Chao.spec[,5:26] <- sapply(Chao.spec[,5:26], as.numeric)

# fit linear models
spec1 <- lm(ChaoRich~mds1+mds2, data=Chao.spec)
spec2 <- lm(ChaoRich~mds1+fire, data=Chao.spec)
spec3 <- lm(ChaoRich~mds1+ndvi, data=Chao.spec)
spec4 <- lm(ChaoRich~mds2+fire, data=Chao.spec)
spec5 <- lm(ChaoRich~mds2+ndvi, data=Chao.spec)
spec6 <- lm(ChaoRich~fire+ndvi, data=Chao.spec)

spec7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.spec)
spec8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.spec)
spec9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.spec)
spec10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.spec)
spec11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.spec)
spec12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.spec)

spec13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.spec)
spec14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.spec)
spec15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.spec)
spec16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.spec)
spec17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.spec)
spec18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.spec)

spec19 <- lm(ChaoRich~mds1, data=Chao.spec)
spec20 <- lm(ChaoRich~mds2, data=Chao.spec)
spec21 <- lm(ChaoRich~fire, data=Chao.spec)
spec22 <- lm(ChaoRich~ndvi, data=Chao.spec)
spec23 <- lm(ChaoRich~elev10, data=Chao.spec)
spec24 <- lm(ChaoRich~chili10, data=Chao.spec)

spec.list1 <- list(spec1,spec2,spec3,spec4,spec5,spec6,spec7,spec8,spec9,spec10,spec11,spec12,spec13,spec14,spec15,spec16,spec17,spec18,spec19,spec20,spec21,spec22,spec23,spec24)
spec.aic <- aictab(spec.list1)

spectop1 <- spec20
spectop2 <- spec5

spectop1.hum <- lm(ChaoRich~mds2+allhum.ln, data=Chao.spec)
spectop2.hum <- lm(ChaoRich~mds2+ndvi+allhum.ln, data=Chao.spec)

spec.list2 <- list(spectop1.hum,spectop2.hum,
                   spectop1,spectop2)
spec.aic.hum <- aictab(spec.list2)

spec.tophum <- spectop1.hum

# get 95% conf int for best human model
confint(spectop1.hum)

########### generalist species richness

gen <- dplyr::select(herp,Point,DailyMeanHuman,UTST,PLSK,SCOC,ASTI,ELMU,HYOC,DIPU,RHLE,LACA,CROR,COLA,PICA)

rich.gen <- specpool(gen[,3:14])
sitegen <- estimateR(gen[,3:14])
Chao.gen <- data.frame("ChaoRich"=sitegen[2,],"ChaoSE"=sitegen[3,],"PointGrp"=gen$Point,"Sp.obs"=sitegen[1,])

Chao.gen <- full_join(Chao.gen,siteCovs2,by="PointGrp")
Chao.gen$ChaoRich <- ifelse(is.na(Chao.gen$ChaoRich),0,Chao.gen$ChaoRich)
Chao.gen$ChaoSE <- ifelse(is.na(Chao.gen$ChaoSE),0,Chao.gen$ChaoSE)
Chao.gen$Sp.obs <- ifelse(is.na(Chao.gen$Sp.obs),0,Chao.gen$Sp.obs)

# scale numeric covariates
Chao.gen[,5:26] <- lapply(Chao.gen[,5:26], scale)
Chao.gen[,5:26] <- sapply(Chao.gen[,5:26], as.numeric)

# fit linear models
gen1 <- lm(ChaoRich~mds1+mds2, data=Chao.gen)
gen2 <- lm(ChaoRich~mds1+fire, data=Chao.gen)
gen3 <- lm(ChaoRich~mds1+ndvi, data=Chao.gen)
gen4 <- lm(ChaoRich~mds2+fire, data=Chao.gen)
gen5 <- lm(ChaoRich~mds2+ndvi, data=Chao.gen)
gen6 <- lm(ChaoRich~fire+ndvi, data=Chao.gen)

gen7 <- lm(ChaoRich~mds1+mds2+chili10, data=Chao.gen)
gen8 <- lm(ChaoRich~mds1+fire+chili10, data=Chao.gen)
gen9 <- lm(ChaoRich~mds1+ndvi+chili10, data=Chao.gen)
gen10 <- lm(ChaoRich~mds2+fire+chili10, data=Chao.gen)
gen11 <- lm(ChaoRich~mds2+ndvi+chili10, data=Chao.gen)
gen12 <- lm(ChaoRich~fire+ndvi+chili10, data=Chao.gen)

gen13 <- lm(ChaoRich~mds1+mds2+elev10, data=Chao.gen)
gen14 <- lm(ChaoRich~mds1+fire+elev10, data=Chao.gen)
gen15 <- lm(ChaoRich~mds1+ndvi+elev10, data=Chao.gen)
gen16 <- lm(ChaoRich~mds2+fire+elev10, data=Chao.gen)
gen17 <- lm(ChaoRich~mds2+ndvi+elev10, data=Chao.gen)
gen18 <- lm(ChaoRich~fire+ndvi+elev10, data=Chao.gen)

gen19 <- lm(ChaoRich~mds1, data=Chao.gen)
gen20 <- lm(ChaoRich~mds2, data=Chao.gen)
gen21 <- lm(ChaoRich~fire, data=Chao.gen)
gen22 <- lm(ChaoRich~ndvi, data=Chao.gen)
gen23 <- lm(ChaoRich~elev10, data=Chao.gen)
gen24 <- lm(ChaoRich~chili10, data=Chao.gen)

gen.list1 <- list(gen1,gen2,gen3,gen4,gen5,gen6,gen7,gen8,gen9,gen10,gen11,gen12,gen13,gen14,gen15,gen16,gen17,gen18,gen19,gen20,gen21,gen22,gen23,gen24)
gen.aic <- aictab(gen.list1)

gentop1 <- gen15
gentop2 <- gen13

gentop1.hum <- lm(ChaoRich~mds1+ndvi+elev10+allhum.ln, data=Chao.gen)
gentop2.hum <- lm(ChaoRich~mds1+mds2+elev10+allhum.ln, data=Chao.gen)


gen.list2 <- list(gentop1.hum,gentop2.hum,
                  gentop1,gentop2)
gen.aic.hum <- aictab(gen.list2)

gen.tophum <- gentop1.hum
summary(gen.tophum)

# get 95% conf int for best human model
confint(gentop1.hum)

########### organize data for body size and body temp analyses

herp3 <- reshape2::melt(herp,id.vars=c("Point","DailyMeanHuman"),value.name="Detections")
herp3 <- filter(herp3, Detections>0)
herp3 <- filter(herp3,variable != "TotalDetections" & variable != "Richness") # this is now a dataframe with a record for each species that was observed at a given site and the total detections of that species

herp4 <- left_join(herp3,herplist[,c(-6)],by=c("variable"="SpCode"))

mean.svl <- herp4 %>%
  filter(SpGroup == "Lizard") %>%
  group_by(Point) %>%
  mutate(SVL.detweighted = SVL_Amburgey*Detections) %>%
  summarise(svlsp.mean = mean(SVL_Amburgey), svlsp.med = median(SVL_Amburgey), svldet.mean = sum(SVL.detweighted)/sum(Detections))

mean.temp <- herp4 %>%
  filter(SpGroup == "Lizard") %>%
  group_by(Point) %>%
  mutate(temp.detweighted = MeanBodyTemp*Detections) %>%
  summarise(temp.spmean = mean(MeanBodyTemp), temp.spmed = median(MeanBodyTemp), tempdet.mean = sum(temp.detweighted)/sum(Detections))

span.temp <- herp4 %>%
  filter(SpGroup == "Lizard") %>%
  group_by(Point) %>%
  mutate(tempspan.detweighted = spanBodyTemp*Detections) %>%
  summarise(tempspan.spmean = mean(spanBodyTemp), tempspan.spmed = median(spanBodyTemp), tempspandet.mean = sum(tempspan.detweighted)/sum(Detections))

herptraits <- left_join(mean.svl,mean.temp,by="Point")
herptraits <- left_join(herptraits,span.temp,by="Point")
herptraits <- left_join(herptraits,siteCovs2,by=c("Point"="PointGrp")) ## VERY IMPORTANT: since there were no herp detections at MTRP08 and BMTP04, they are not included in this dataframe. This is appropriate for the body size and body temp models since n/as will be dropped anyway, but may not be for other things.

#standardize covariates
herptraits.unsc <- herptraits
herptraits[,11:32] <- scale(herptraits[,11:32])

########### body size linear regressions 

## linear regressions using mean det body size as response
# fit models with up to 2 habitat vars and 1 (or no) topo var
sizedet1 <- lm(svldet.mean~mds1+mds2, data=herptraits)
sizedet2 <- lm(svldet.mean~mds1+fire, data=herptraits)
sizedet3 <- lm(svldet.mean~mds1+ndvi, data=herptraits)
sizedet4 <- lm(svldet.mean~mds2+fire, data=herptraits)
sizedet5 <- lm(svldet.mean~mds2+ndvi, data=herptraits)
sizedet6 <- lm(svldet.mean~fire+ndvi, data=herptraits)

sizedet7 <- lm(svldet.mean~mds1+mds2+chili10, data=herptraits)
sizedet8 <- lm(svldet.mean~mds1+fire+chili10, data=herptraits)
sizedet9 <- lm(svldet.mean~mds1+ndvi+chili10, data=herptraits)
sizedet10 <- lm(svldet.mean~mds2+fire+chili10, data=herptraits)
sizedet11 <- lm(svldet.mean~mds2+ndvi+chili10, data=herptraits)
sizedet12 <- lm(svldet.mean~fire+ndvi+chili10, data=herptraits)

sizedet13 <- lm(svldet.mean~mds1+mds2+elev10, data=herptraits)
sizedet14 <- lm(svldet.mean~mds1+fire+elev10, data=herptraits)
sizedet15 <- lm(svldet.mean~mds1+ndvi+elev10, data=herptraits)
sizedet16 <- lm(svldet.mean~mds2+fire+elev10, data=herptraits)
sizedet17 <- lm(svldet.mean~mds2+ndvi+elev10, data=herptraits)
sizedet18 <- lm(svldet.mean~fire+ndvi+elev10, data=herptraits)

sizedet19 <- lm(svldet.mean~mds1, data=herptraits)
sizedet20 <- lm(svldet.mean~mds2, data=herptraits)
sizedet21 <- lm(svldet.mean~fire, data=herptraits)
sizedet22 <- lm(svldet.mean~ndvi, data=herptraits)
sizedet23 <- lm(svldet.mean~elev10, data=herptraits)
sizedet24 <- lm(svldet.mean~chili10, data=herptraits)

sizedet.list <- list(sizedet1,sizedet2,sizedet3,sizedet4,sizedet5,sizedet6,sizedet7,sizedet8,sizedet9,sizedet10,sizedet11,sizedet12,sizedet13,sizedet14,sizedet15,sizedet16,sizedet17,sizedet18,sizedet19,sizedet20,sizedet21,sizedet22,sizedet23,sizedet24)
sizedet.aic <- aictab(sizedet.list)

sizedettop1 <- sizedet12
sizedettop2 <- sizedet9
sizedettop3 <- sizedet17
sizedettop4 <- sizedet6
sizedettop5 <- sizedet3

sizedettop1.hum <- lm(svldet.mean~fire+ndvi+chili10+allhum.ln, data=herptraits)
sizedettop1.cyc <- lm(svldet.mean~fire+ndvi+chili10+cyc.ln, data=herptraits)
sizedettop1.hik <- lm(svldet.mean~fire+ndvi+chili10+hik.ln, data=herptraits)

sizedettop3.hum <- lm(svldet.mean~mds2+ndvi+elev10+allhum.ln, data=herptraits)
sizedettop3.cyc <- lm(svldet.mean~mds2+ndvi+elev10+cyc.ln, data=herptraits)
sizedettop3.hik <- lm(svldet.mean~mds2+ndvi+elev10+hik.ln, data=herptraits)

sizedettop4.hum <- lm(svldet.mean~fire+ndvi+allhum.ln, data=herptraits)
sizedettop4.cyc <- lm(svldet.mean~fire+ndvi+cyc.ln, data=herptraits)
sizedettop4.hik <- lm(svldet.mean~fire+ndvi+hik.ln, data=herptraits)

sizedettop2.hum <- lm(svldet.mean~mds1+ndvi+chili10+allhum.ln, data=herptraits)
sizedettop2.cyc <- lm(svldet.mean~mds1+ndvi+chili10+cyc.ln, data=herptraits)
sizedettop2.hik <- lm(svldet.mean~mds1+ndvi+chili10+hik.ln, data=herptraits)

sizedettop5.hum <- lm(svldet.mean~mds1+ndvi+allhum.ln, data=herptraits)
sizedettop5.cyc <- lm(svldet.mean~mds1+ndvi+cyc.ln, data=herptraits)
sizedettop5.hik <- lm(svldet.mean~mds1+ndvi+hik.ln, data=herptraits)

sizedet.list2 <- list("sizedettop1.hum"=sizedettop1.hum,
                      "sizedettop1.cyc"=sizedettop1.cyc,
                      "sizedettop1.hik"=sizedettop1.hik,
                      "sizedettop1"=sizedettop1,
                      "sizedettop2.hum"=sizedettop2.hum,
                      "sizedettop2.cyc"=sizedettop2.cyc,
                      "sizedettop2.hik"=sizedettop2.hik,
                      "sizedettop2"=sizedettop2,
                      "sizedettop3.hum"=sizedettop3.hum,
                      "sizedettop3.cyc"=sizedettop3.cyc,
                      "sizedettop3.hik"=sizedettop3.hik,
                      "sizedettop3"=sizedettop3,
                      "sizedettop4.hum"=sizedettop4.hum,
                      "sizedettop4.cyc"=sizedettop4.cyc,
                      "sizedettop4.hik"=sizedettop4.hik,
                      "sizedettop4"=sizedettop4,
                      "sizedettop5.hum"=sizedettop5.hum,
                      "sizedettop5.cyc"=sizedettop5.cyc,
                      "sizedettop5.hik"=sizedettop5.hik,
                      "sizedettop5"=sizedettop5)

sizedet.aic.hum <- aictab(sizedet.list2)

sizedet.tophum <- sizedettop4.hum

summary(sizedettop4.hum) #for R2 of best model

# get 95% conf int for best human model
confint(sizedet.tophum)

########### linear regressions using mean sp body size as response

# fit models with up to 2 habitat vars and 1 (or no) topo var
sizesp1 <- lm(svlsp.mean~mds1+mds2, data=herptraits)
sizesp2 <- lm(svlsp.mean~mds1+fire, data=herptraits)
sizesp3 <- lm(svlsp.mean~mds1+ndvi, data=herptraits)
sizesp4 <- lm(svlsp.mean~mds2+fire, data=herptraits)
sizesp5 <- lm(svlsp.mean~mds2+ndvi, data=herptraits)
sizesp6 <- lm(svlsp.mean~fire+ndvi, data=herptraits)

sizesp7 <- lm(svlsp.mean~mds1+mds2+chili10, data=herptraits)
sizesp8 <- lm(svlsp.mean~mds1+fire+chili10, data=herptraits)
sizesp9 <- lm(svlsp.mean~mds1+ndvi+chili10, data=herptraits)
sizesp10 <- lm(svlsp.mean~mds2+fire+chili10, data=herptraits)
sizesp11 <- lm(svlsp.mean~mds2+ndvi+chili10, data=herptraits)
sizesp12 <- lm(svlsp.mean~fire+ndvi+chili10, data=herptraits)

sizesp13 <- lm(svlsp.mean~mds1+mds2+elev10, data=herptraits)
sizesp14 <- lm(svlsp.mean~mds1+fire+elev10, data=herptraits)
sizesp15 <- lm(svlsp.mean~mds1+ndvi+elev10, data=herptraits)
sizesp16 <- lm(svlsp.mean~mds2+fire+elev10, data=herptraits)
sizesp17 <- lm(svlsp.mean~mds2+ndvi+elev10, data=herptraits)
sizesp18 <- lm(svlsp.mean~fire+ndvi+elev10, data=herptraits)

sizesp19 <- lm(svlsp.mean~mds1, data=herptraits)
sizesp20 <- lm(svlsp.mean~mds2, data=herptraits)
sizesp21 <- lm(svlsp.mean~fire, data=herptraits)
sizesp22 <- lm(svlsp.mean~ndvi, data=herptraits)
sizesp23 <- lm(svlsp.mean~elev10, data=herptraits)
sizesp24 <- lm(svlsp.mean~chili10, data=herptraits)

sizesp.list <- list(sizesp1,sizesp2,sizesp3,sizesp4,sizesp5,sizesp6,sizesp7,sizesp8,sizesp9,sizesp10,sizesp11,sizesp12,sizesp13,sizesp14,sizesp15,sizesp16,sizesp17,sizesp18,sizesp19,sizesp20,sizesp21,sizesp22,sizesp23,sizesp24)
sizesp.aic <- aictab(sizesp.list)

sizesptop2 <- sizesp3
sizesptop1 <- sizesp15
sizesptop3 <- sizesp22

sizesptop3.hum <- lm(svlsp.mean~ndvi+allhum.ln, data=herptraits)
sizesptop3.cyc <- lm(svlsp.mean~ndvi+cyc.ln, data=herptraits)
sizesptop3.hik <- lm(svlsp.mean~ndvi+hik.ln, data=herptraits)

sizesptop2.hum <- lm(svlsp.mean~mds1+ndvi+allhum.ln, data=herptraits)
sizesptop2.cyc <- lm(svlsp.mean~mds1+ndvi+cyc.ln, data=herptraits)
sizesptop2.hik <- lm(svlsp.mean~mds1+ndvi+hik.ln, data=herptraits)

sizesptop1.hum <- lm(svlsp.mean~mds1+ndvi+elev10+allhum.ln, data=herptraits)
sizesptop1.cyc <- lm(svlsp.mean~mds1+ndvi+elev10+cyc.ln, data=herptraits)
sizesptop1.hik <- lm(svlsp.mean~mds1+ndvi+elev10+hik.ln, data=herptraits)

sizesp.list2 <- list("NDVI + **total human**"=sizesptop3.hum,
                     "NDVI + **cyclist**"=sizesptop3.cyc,
                     "NDVI + **pedestrian**"= sizesptop3.hik,
                     "NDVI"=sizesptop3,
                     "nMDS1 + NDVI + **total human**"=sizesptop2.hum,
                     "nMDS1 + NDVI + **cyclist**"=sizesptop2.cyc,
                     "nMDS1 + NDVI + **pedestrian**"=sizesptop2.hik,
                     "nMDS1 + NDVI"=sizesptop2,
                     "nMDS1 + NDVI + elevation + **total human**"=sizesptop1.hum,
                     "nMDS1 + NDVI + elevation + **cyclist**"=sizesptop1.cyc,
                     "nMDS1 + NDVI + elevation + **pedestrian**"=sizesptop1.hik,
                     "nMDS1 + NDVI + elevation"=sizesptop1)

sizesp.aic.hum <- aictab(sizesp.list2)
sizesp.tophum <- sizesptop1.hum

summary(sizesptop1.hik) #for R2 of best model

# get 95% conf int for best human model
confint(sizesp.tophum)

########### linear regressions using span det body temp as response

tempdet1 <- lm(tempspandet.mean~mds1+mds2, data=herptraits)
tempdet2 <- lm(tempspandet.mean~mds1+fire, data=herptraits)
tempdet3 <- lm(tempspandet.mean~mds1+ndvi, data=herptraits)
tempdet4 <- lm(tempspandet.mean~mds2+fire, data=herptraits)
tempdet5 <- lm(tempspandet.mean~mds2+ndvi, data=herptraits)
tempdet6 <- lm(tempspandet.mean~fire+ndvi, data=herptraits)
tempdet7 <- lm(tempspandet.mean~mds1+mds2+chili10, data=herptraits)
tempdet8 <- lm(tempspandet.mean~mds1+fire+chili10, data=herptraits)
tempdet9 <- lm(tempspandet.mean~mds1+ndvi+chili10, data=herptraits)
tempdet10 <- lm(tempspandet.mean~mds2+fire+chili10, data=herptraits)
tempdet11 <- lm(tempspandet.mean~mds2+ndvi+chili10, data=herptraits)
tempdet12 <- lm(tempspandet.mean~fire+ndvi+chili10, data=herptraits)
tempdet13 <- lm(tempspandet.mean~mds1+mds2+elev10, data=herptraits)
tempdet14 <- lm(tempspandet.mean~mds1+fire+elev10, data=herptraits)
tempdet15 <- lm(tempspandet.mean~mds1+ndvi+elev10, data=herptraits)
tempdet16 <- lm(tempspandet.mean~mds2+fire+elev10, data=herptraits)
tempdet17 <- lm(tempspandet.mean~mds2+ndvi+elev10, data=herptraits)
tempdet18 <- lm(tempspandet.mean~fire+ndvi+elev10, data=herptraits)
tempdet19 <- lm(tempspandet.mean~mds1, data=herptraits)
tempdet20 <- lm(tempspandet.mean~mds2, data=herptraits)
tempdet21 <- lm(tempspandet.mean~fire, data=herptraits)
tempdet22 <- lm(tempspandet.mean~ndvi, data=herptraits)
tempdet23 <- lm(tempspandet.mean~elev10, data=herptraits)
tempdet24 <- lm(tempspandet.mean~chili10, data=herptraits)

tempdet.list <- list(tempdet1,tempdet2,tempdet3,tempdet4,tempdet5,tempdet6,tempdet7,tempdet8,tempdet9,tempdet10,tempdet11,tempdet12,tempdet13,tempdet14,tempdet15,tempdet16,tempdet17,tempdet18,tempdet19,tempdet20,tempdet21,tempdet22,tempdet23,tempdet24)
tempdet.aic <- aictab(tempdet.list)

tempdettop1 <- tempdet1
tempdettop2 <- tempdet13
tempdettop3 <- tempdet7

tempdettop1.hum <- lm(tempspandet.mean~mds1+mds2+allhum.ln, data=herptraits)
tempdettop1.cyc <- lm(tempspandet.mean~mds1+mds2+cyc.ln, data=herptraits)
tempdettop1.hik <- lm(tempspandet.mean~mds1+mds2+hik.ln, data=herptraits)

tempdettop2.hum <- lm(tempspandet.mean~mds1+mds2+elev10+allhum.ln, data=herptraits)
tempdettop2.cyc <- lm(tempspandet.mean~mds1+mds2+elev10+cyc.ln, data=herptraits)
tempdettop2.hik <- lm(tempspandet.mean~mds1+mds2+elev10+hik.ln, data=herptraits)

tempdettop3.hum <- lm(tempspandet.mean~mds1+mds2+chili10+allhum.ln, data=herptraits)
tempdettop3.cyc <- lm(tempspandet.mean~mds1+mds2+chili10+cyc.ln, data=herptraits)
tempdettop3.hik <- lm(tempspandet.mean~mds1+mds2+chili10+hik.ln, data=herptraits)

tempdet.list2 <- list("tempdettop1.hum"=tempdettop1.hum,
                      "tempdettop1.cyc"=tempdettop1.cyc,
                      "tempdettop1.hik"=tempdettop1.hik,
                      "tempdettop1"=tempdettop1,
                      "tempdettop2.hum"=tempdettop2.hum,
                      "tempdettop2.cyc"=tempdettop2.cyc,
                      "tempdettop2.hik"=tempdettop2.hik,
                      "tempdettop2"=tempdettop2,
                      "tempdettop3.hum"=tempdettop3.hum,
                      "tempdettop3.cyc"=tempdettop3.cyc,
                      "tempdettop3.hik"=tempdettop3.hik,
                      "tempdettop3"=tempdettop3)

tempdet.aic.hum <- aictab(tempdet.list2)
tempdet.tophum <- tempdettop1.hum

summary(tempdettop1.hum) #for R2 of best model

confint(tempdet.tophum)

########### linear regressions using range in species body temp as response

tempsp1 <- lm(tempspan.spmean~mds1+mds2, data=herptraits)
tempsp2 <- lm(tempspan.spmean~mds1+fire, data=herptraits)
tempsp3 <- lm(tempspan.spmean~mds1+ndvi, data=herptraits)
tempsp4 <- lm(tempspan.spmean~mds2+fire, data=herptraits)
tempsp5 <- lm(tempspan.spmean~mds2+ndvi, data=herptraits)
tempsp6 <- lm(tempspan.spmean~fire+ndvi, data=herptraits)
tempsp7 <- lm(tempspan.spmean~mds1+mds2+chili10, data=herptraits)
tempsp8 <- lm(tempspan.spmean~mds1+fire+chili10, data=herptraits)
tempsp9 <- lm(tempspan.spmean~mds1+ndvi+chili10, data=herptraits)
tempsp10 <- lm(tempspan.spmean~mds2+fire+chili10, data=herptraits)
tempsp11 <- lm(tempspan.spmean~mds2+ndvi+chili10, data=herptraits)
tempsp12 <- lm(tempspan.spmean~fire+ndvi+chili10, data=herptraits)
tempsp13 <- lm(tempspan.spmean~mds1+mds2+elev10, data=herptraits)
tempsp14 <- lm(tempspan.spmean~mds1+fire+elev10, data=herptraits)
tempsp15 <- lm(tempspan.spmean~mds1+ndvi+elev10, data=herptraits)
tempsp16 <- lm(tempspan.spmean~mds2+fire+elev10, data=herptraits)
tempsp17 <- lm(tempspan.spmean~mds2+ndvi+elev10, data=herptraits)
tempsp18 <- lm(tempspan.spmean~fire+ndvi+elev10, data=herptraits)
tempsp19 <- lm(tempspan.spmean~mds1, data=herptraits)
tempsp20 <- lm(tempspan.spmean~mds2, data=herptraits)
tempsp21 <- lm(tempspan.spmean~fire, data=herptraits)
tempsp22 <- lm(tempspan.spmean~ndvi, data=herptraits)
tempsp23 <- lm(tempspan.spmean~elev10, data=herptraits)
tempsp24 <- lm(tempspan.spmean~chili10, data=herptraits)

tempsp.list <- list(tempsp1,tempsp2,tempsp3,tempsp4,tempsp5,tempsp6,tempsp7,tempsp8,tempsp9,tempsp10,tempsp11,tempsp12,tempsp13,tempsp14,tempsp15,tempsp16,tempsp17,tempsp18,tempsp19,tempsp20,tempsp21,tempsp22,tempsp23,tempsp24)
tempsp.aic <- aictab(tempsp.list)

tempsptop1 <- tempsp19
tempsptop3 <- tempsp1

tempsptop1.hum <- lm(tempspan.spmean~mds1+allhum.ln, data=herptraits)
tempsptop1.cyc <- lm(tempspan.spmean~mds1+cyc.ln, data=herptraits)
tempsptop1.hik <- lm(tempspan.spmean~mds1+hik.ln, data=herptraits)

tempsptop3.hum <- lm(tempspan.spmean~mds1+mds2+allhum.ln, data=herptraits)
tempsptop3.cyc <- lm(tempspan.spmean~mds1+mds2+cyc.ln, data=herptraits)
tempsptop3.hik <- lm(tempspan.spmean~mds1+mds2+hik.ln, data=herptraits)

tempsp.list2 <- list("tempsptop1.hum"=tempsptop1.hum,
                     "tempsptop1.cyc"=tempsptop1.cyc,
                     "tempsptop1.hik"=tempsptop1.hik,
                     "tempsptop1"=tempsptop1,
                     "tempsptop3.hum"=tempsptop3.hum,
                     "tempsptop3.cyc"=tempsptop3.cyc,
                     "tempsptop3.hik"=tempsptop3.hik,
                     "tempsptop3"=tempsptop3)

tempsp.aic.hum <- aictab(tempsp.list2)

tempsp.tophum <- tempsptop1.hum

summary(tempsptop1.hik) #for R2 of best model

confint(tempsp.tophum)