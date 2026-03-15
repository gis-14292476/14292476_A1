setwd("C:/Users/86153/Desktop/work4m/SE/ass1/14292476_A1")
library(terra)
library(sf)
library(spatstat) 

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

#### 1 data processing ####

# 1.1 read data
## Melesmeles spatial point data （British National Grid）----

Melesmeles = read.csv("./data/Melesmeles.csv")

## data clean
Melesmeles<-Melesmeles[!is.na(Melesmeles$Latitude),]
Melesmeles = Melesmeles[Melesmeles$Coordinate.uncertainty_m < 5000, ]

### !!! add discussion about 3000 

## convert to spatial point
Melesmeles_latlong = data.frame( x = Melesmeles$Longitude, y = Melesmeles$Latitude )
Melesmeles_sp = st_as_sf( Melesmeles_latlong, coords = c("x", "y"), crs = "epsg:4326")


## England shapefile ----
england = st_read('./data/England.shp')

## LCM raster (British National Grid) ----
LCM=rast("./data/LCMUK.tif")

## check projection ----
england=st_transform(england,crs(LCM))
Melesmeles_sp = st_transform(Melesmeles_sp, crs(LCM))

# 1.2 data processing

## LCM raster ----
LCM=crop(LCM,st_buffer(england, dist= 1000))
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
# add why choose fact=4

## Cull data（england）----
MelesmelesFin=Melesmeles_sp[england,]
LCM=crop(LCM,england,mask=TRUE)


## Calculate environment variables ----
reclassLeaf  = c(0,1,rep(0,20))
reclassUrban = c(rep(0,19),1,1)

broadleaf = reclass_binary(LCM, reclassLeaf)
urban = reclass_binary(LCM, reclassUrban)

lcm_wood_1800  = calc_prop_focal(broadleaf, 1800)
lcm_urban_2300 = calc_prop_focal(urban, 2300)

# add why choose 1800/2300

# 
#### draw picture ####

# figure 1 data show
plot(st_geometry(MelesmelesFin), add=TRUE, col="red")
