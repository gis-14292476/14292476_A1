setwd("C:/Users/86153/Desktop/work4m/SE/ass1/14292476_A1")
library(terra)
library(sf)
library(spatstat) 
library(mlr)
#### function ####

calc_prop_focal <- function(r, radius){
  
  # 1 calculate the num of pixel
  nPix = round(radius / res(r)[1])
  nPix = (nPix * 2) + 1
  
  # 2 create matrix
  weightsMatrix = matrix(1:nPix^2, nrow=nPix, ncol=nPix)
  
  # 3 central pixel
  x = ceiling(ncol(weightsMatrix)/2)
  y = ceiling(nrow(weightsMatrix)/2)
  
  focalCell = weightsMatrix[x,y]
  indFocal = which(weightsMatrix == focalCell, arr.ind = TRUE)
  
  # 4 calculate distance
  distances = list()
  
  for(i in 1:nPix^2){
    
    ind.i = which(weightsMatrix == i, arr.ind = TRUE)
    
    diffX = abs(ind.i[1,1] - indFocal[1,1]) * res(r)[1]
    diffY = abs(ind.i[1,2] - indFocal[1,2]) * res(r)[1]
    
    dist.i = sqrt(diffX^2 + diffY^2)
    
    distances[[i]] = dist.i
  }
  
  # 5 back to matrix
  weightsMatrix[] = unlist(distances)
  
  # 6 delete outside pixel
  weightsMatrix[weightsMatrix > radius] = NA
  
  # 7 weight normalization
  weightsMatrix[!is.na(weightsMatrix)] =
    1 / length(weightsMatrix[!is.na(weightsMatrix)])
  
  # 8 focal calculate
  result = focal(r, w = weightsMatrix, fun = "sum")
  
  return(result)
}

reclass_binary <- function(lcm, rule){
  
  # transformation
  lcm = as.factor(lcm$LCMUK_1)
  
  # create reclass matrix
  RCmatrix = cbind(levels(lcm)[[1]], rule)
  
  # old class + new class
  RCmatrix = RCmatrix[,2:3]
  
  # to numeric
  RCmatrix = apply(RCmatrix, 2, as.numeric)
  
  result = classify(lcm, RCmatrix)
  
  return(result)
}



#### ----------------- 1 data processing ----------------- ####

# 1.1 read data ----
## Melesmeles spatial point data （British National Grid）----

Melesmeles = read.csv("./data/Melesmeles.csv")

## data clean
Melesmeles<-Melesmeles[!is.na(Melesmeles$Latitude),]
Melesmeles = Melesmeles[Melesmeles$Coordinate.uncertainty_m < 5000, ]

### !!! add discussion about 3000 

## convert to spatial point
Melesmeles_latlong = data.frame( x = Melesmeles$Longitude, y = Melesmeles$Latitude )
Melesmeles_sp = st_as_sf( Melesmeles_latlong, coords = c("x", "y"), crs = "epsg:4326")

## scot shapefile ----
scot=st_read('./data/scotSamp.shp')

## LCM raster (British National Grid) ----
LCM=rast("./data/LCMUK.tif")

## scotDEM ----
demScot=rast('./data/demScotland.tif')

## check projection ----
scot=st_transform(scot,crs(LCM))
Melesmeles_sp = st_transform(Melesmeles_sp, crs(LCM))


# 1.2 data processing ----
## LCM raster ----
LCM=crop(LCM,st_buffer(scot, dist= 1000))
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
# add why choose fact=4

## Cull data（england）----
MelesmelesFin=Melesmeles_sp[scot,]
LCM=crop(LCM,scot,mask=TRUE)

## Calculate environment variables ----
# Foraging resources
reclassForage = c(0,1,0,1,1,rep(0,17))
forage = reclass_binary(LCM, reclassForage)

# Disturbance
reclassDisturb = c(rep(0,20),1,1)
disturb = reclass_binary(LCM, reclassDisturb)

# Refuge habitat
reclassRefuge = c(0,1,1,rep(0,6),1,rep(0,12))
refuge = reclass_binary(LCM, reclassRefuge)


# --- Landscape proportion (moving window) ---

lcm_forage_2000  = calc_prop_focal(forage, 2000)
lcm_disturb_2300 = calc_prop_focal(disturb, 2300)
lcm_refuge_800  = calc_prop_focal(refuge, 800)

# add why choose 1800/2300

demScot=terra::resample(demScot,lcm_forage_2000)

#### ----------------- 2 create dataset -----------------  ####

# 2.1 add all environment variables ----
allEnv=c(lcm_forage_2000,lcm_disturb_2300,lcm_refuge_800,demScot)
names(allEnv)=c('food','urban','habitat','dem')

# 2.2 ‘presence’ and ‘background’ -----
## create the background points----
set.seed(33)

back = spatSample(
  allEnv,size=3000,
  as.points=TRUE,
  method="random",
  na.rm=TRUE) 
# add why choose 2000

back=back[!is.na(back$food),]
### ???

back=st_as_sf(back,crs="EPSG:27700")

## organize ‘presence’ and ‘background’----

# extract env var for presence point (MelesmelesFin)
eP=terra::extract(allEnv,MelesmelesFin)
# env var + presence point 
Pres.cov=st_as_sf(cbind(eP,MelesmelesFin))


## add text for differencr point ----

# add text for presence point (MelesmelesFin)
Pres.cov$Pres=1
# delet ID
Pres.cov=Pres.cov[,-1]

# add text for background point
Back.cov=st_as_sf(data.frame(back,Pres=0))


# 2.3 merging data ----
## merging point coordinates ----
# get coordinate
coordsPres=st_coordinates(Pres.cov)
coordsBack=st_coordinates(back)

# merging point coordinates
coords=data.frame(rbind(coordsPres,coordsBack))
colnames(coords)=c("x","y")

all.cov=rbind(Pres.cov,Back.cov)
all.cov=cbind(all.cov,coords)


## data cleaning ----
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)


####----------------- 3 define parameters ----------------- ####

# 3.1 Create category task ----
task=all.cov

task$Pres=as.factor(task$Pres)

task = makeClassifTask(
  data = task[,c('food','urban','habitat','dem',"Pres")],
  target = "Pres",
  positive = "1",
  coordinates = task[,c("x","y")])


# 3.2 Define the cross-validation strategy ----

# Randomly shuffle all points, disregarding spatial location
perf_levelCV = makeResampleDesc(
  method = "RepCV",
  predict = "test",
  folds = 5,
  reps = 5)

# Training and test points are spatially separated.
# Evaluate the model's predictive ability in "unsampled regions"
perf_level_spCV = makeResampleDesc(
  method = "SpRepCV",
  folds = 5, 
  reps = 5) 

