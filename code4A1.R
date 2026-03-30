#### ----------------- 0 Environment configuration ----------------- ####

setwd("C:/Users/86153/Desktop/work4m/SE/ass1/14292476_A1")

library(terra)
library(sf)
library(spatstat) 
library(mlr)
library(dismo)
library(maxnet)
library(glmnet)
library(precrec)
library(car)

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

raster.as.im = function(im) {
  # Get raster resolution (cell size)
  r = raster::res(im)[1]
  
  # Get raster origin (lower-left corner coordinates)
  orig = ext(im)[c(1,3)]
  
  # Construct x coordinates (column direction)
  xx = orig[1] + seq(from = 0, to = (ncol(im) - 1) * 100, by = r)
  
  # Construct y coordinates (row direction)
  yy = orig[2] + seq(from = 0, to = (nrow(im) - 1) * 100, by = r)
  
  # Create value matrix from raster
  # - fill by row
  # - flip vertically (because raster origin is bottom-left,
  #   while spatstat expects top-left)
  mat = matrix(raster::values(im),
               ncol = ncol(im),
               nrow = nrow(im),
               byrow = TRUE)[nrow(im):1, ]
  
  # Convert to spatstat 'im' object
  return(spatstat.geom::im(mat,
                           xcol = xx,
                           yrow = yy))
}

landBuffer <- function(speciesData, r, landcover){         
  
  speciesBuffer <- st_buffer(speciesData, dist=r)                     
  
  bufferlandcover <- crop(landcover, speciesBuffer)              
  
  masklandcover <- extract(bufferlandcover, speciesBuffer,fun="sum")      
  
  landcoverArea <- masklandcover$LCMUK_1*625  
  
  percentcover <- landcoverArea/as.numeric(st_area(speciesBuffer))*100 
  
  return(percentcover)                                       
}

Test_GLM <- function(resList,speciesData,radii){
  # list to dataframe
  resFin=do.call("cbind",resList)
  glmData=data.frame(resFin)
  #intuitive column names
  colnames(glmData)=paste("radius",radii,sep="")
  #add in the presences data
  glmData$Pres<-speciesData$Pres
  
  glmRes=data.frame(radius=NA,loglikelihood=NA)
  
  for(i in radii){
    
    #build the model formula with format "response variable ~ explanatory variable, family, data"
    #Here "binomial" is the error distribution normally used for a binary outcome (e.g. 0, 1) 
    n.i=paste0("Pres~","radius",i,sep ="")
    
    glm.i=glm(formula(n.i),family = "binomial",data = glmData)
    
    ll.i=as.numeric(logLik(glm.i))
    
    glmRes=rbind(glmRes,c(i,ll.i))
  }
  return(glmRes)
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

# draw (Distribution of Coordinate Uncertainty)----
unc <- Melesmeles$Coordinate.uncertainty_m
breaks <- c(1000, 2000, 5000, 10000, Inf) 
unc_cut <- cut(unc, breaks = breaks, right = FALSE)
percent <- table(unc_cut) / length(unc) * 100
percent_nonzero <- percent[percent > 0]

barplot(percent_nonzero,
        col = "#4C72B0",
        border = NA,
        main = "Distribution of Coordinate Uncertainty",
        ylab = "Percentage (%)",
        las = 1,
        cex.names = 0.5)

Melesmeles = Melesmeles[Melesmeles$Coordinate.uncertainty_m <= 1000, ]

# coord_dataframe
Melesmeles_latlong = data.frame(
  x = Melesmeles$Longitude, 
  y = Melesmeles$Latitude )
# point object
Melesmeles_sp = vect(
  Melesmeles_latlong,
  geom = c("x","y"),
  crs = "EPSG:4326"  
)


## scot shapefile ----
scot=st_read('./data/scotSamp.shp')

## LCM raster (British National Grid) ----
LCM=rast("./data/LCMUK.tif")

## scotDEM ----
demScot=rast('./data/demScotland.tif')


## check projection and Cull data----
scot=st_transform(scot,crs(LCM))
studyExtent <- aggregate(vect(scot))

Melesmeles_sp<-project(Melesmeles_sp,crs(LCM))
MelesmelesFin<-crop(Melesmeles_sp,studyExtent)

# 1.2 data processing ----
## LCM raster ----
LCM = LCM$LCMUK_1
# get coords
coords <- crds(MelesmelesFin)
# buffer
studyExtent_buffer <- ext(min(coords[,1]) - 5000,max(coords[,1]) + 5000,min(coords[,2]) - 5000,max(coords[,2]) + 5000)

LCM <- crop(LCM, studyExtent_buffer)

LCM=aggregate(LCM,fact=4,fun="modal")
# add why choose fact=4

## create the background points----
set.seed(11)

back = spatSample(
  studyExtent,
  size=2500,
  method="random") 

plot(LCM)                          # 先画底图（栅格）
plot(scot,col='white', add=TRUE)
plot(back, add=TRUE, col="blue", pch=16, cex=0.5)   # 再画背景点
plot(MelesmelesFin, add=TRUE, col="red", pch=16, cex=0.5) # 最后画出现点


# 1.3 choose feature ----

# extract infor
eA_1<-extract(LCM,back)
eP_1<-extract(LCM,MelesmelesFin)

# get frequency
table(eA_1[,2])
table(eP_1[,2])

# calculate density
hA <- hist(eA_1[,2], breaks = c(0:21), plot = FALSE)
hP <- hist(eP_1[,2], breaks = c(0:21), plot = FALSE)

# calculate difference
diff_density <- hP$density - hA$density

# draw ----
bar_colors <- ifelse(diff_density > 0, "#2E8B57", "#A0522D")

bp <- barplot(diff_density,
              names.arg = 1:21,
              ylim = c(-0.2, 0.3),
              col = bar_colors,
              border = NA,
              space = 0.3,
              main = "Density Difference (Presence vs Background)",
              xlab = "Land Cover Class",
              ylab = "Density Difference",
              axes = FALSE)

grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
abline(v = bp, col = "lightgray", lty = "dotted")

axis(2, las = 1, cex.axis = 0.9)
abline(h = 0, lwd = 2, col = "gray40")

legend("topright",
       legend = c("Presence > Background", "Presence < Background"),
       fill = c("#2E8B57", "#A0522D"),
       bty = "n",
       cex = 0.9)




#### ----------------- 2 create dataset -----------------  ####

Abs<-data.frame(crds(back),Pres=0)
Pres<-data.frame(crds(MelesmelesFin),Pres=1)

MelesmelesData<-rbind(Pres,Abs)
MelesmelesData_sf=st_as_sf(MelesmelesData,coords=c("x","y"),crs="EPSG:27700")

## 2.1 Calculate environment variables ----
LCM<-as.factor(LCM)

# restricted area:Mountain, heath, bog
reclassrestrict = c(rep(0,9),rep(1,4),rep(0,9))
restrict = reclass_binary(LCM, reclassrestrict)

# Broadleaved
reclassBroadleaved = c(0,1,rep(0,20))
Broadleaved = reclass_binary(LCM, reclassBroadleaved)
plot(Broadleaved)

# Coniferouswoodland
reclassConiferous = c(0,0,1,rep(0,19))
Coniferous = reclass_binary(LCM, reclassConiferous)
plot(Coniferous)

# Arable and horticulture
reclassArable = c(rep(0,3),1,rep(0,18))
Arable = reclass_binary(LCM, reclassArable)

# Improved grassland
reclassImpgrass = c(rep(0,4),1,rep(0,17))
Impgrass = reclass_binary(LCM, reclassImpgrass)

# Acid grassland
reclassAcidgrass = c(rep(0,7),1,rep(0,14))
Acidgrass = reclass_binary(LCM, reclassAcidgrass)


## 2.2 buffer ----

# Broadleaved# ----
radii_Broadleaved<-seq(300,3000,by=300)
resList_Broadleaved=list()

for(i in radii_Broadleaved){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=Broadleaved)
  resList_Broadleaved[[i/200]]=res.i
  print(i)
}
glmRes_Broadleaved <- Test_GLM(
  resList = resList_Broadleaved,
  speciesData = MelesmelesData,
  radii = radii_Broadleaved
)
plot(glmRes_Broadleaved$radius, glmRes_Broadleaved$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")
glmRes_Broadleaved


# Arable and horticulture ----
radii_Arable<-seq(500,2000,by=100)
resList_Arable=list()

for(i in radii_Arable){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=Arable)
  resList_Arable[[i/100]]=res.i
  print(i)
}
glmRes_Arable <- Test_GLM(
  resList = resList_Arable,
  speciesData = MelesmelesData,
  radii = radii_Arable
)
plot(glmRes_Arable$radius, glmRes_Arable$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")
glmRes_Arable


# Improved grassland# ----
radii_Impgrass<-seq(200,1800,by=100)
resList_Impgrass=list()

for(i in radii_Impgrass){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=Impgrass)
  resList_Impgrass[[i/100]]=res.i
  print(i)
}
glmRes_Impgrass <- Test_GLM(
  resList = resList_Impgrass,
  speciesData = MelesmelesData,
  radii = radii_Impgrass
)
plot(glmRes_Impgrass$radius, glmRes_Impgrass$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")
glmRes_Impgrass


# Acid grassland# ----
radii_Acidgrass<-seq(400,2400,by=200)
resList_Acidgrass=list()

for(i in radii_Acidgrass){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=Acidgrass)
  resList_Acidgrass[[i/200]]=res.i
  print(i)
}
glmRes_Acidgrass <- Test_GLM(
  resList = resList_Acidgrass,
  speciesData = MelesmelesData,
  radii = radii_Acidgrass
)
plot(glmRes_Acidgrass$radius, glmRes_Acidgrass$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")
glmRes_Acidgrass


# restrict# ----
radii_restrict<-seq(100,1500,by=100)
resList_restrict=list()

for(i in radii_restrict){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=restrict)
  resList_restrict[[i/100]]=res.i
  print(i)
}
glmRes_restrict <- Test_GLM(
  resList = resList_restrict,
  speciesData = MelesmelesData,
  radii = radii_restrict
)
plot(glmRes_restrict$radius, glmRes_restrict$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "buffer", ylab = "logLik")

#Landscape proportion (moving window) ----


lcm_restrict  = calc_prop_focal(restrict, glmRes_restrict$radius[which.max(glmRes_restrict$loglikelihood)])
lcm_Acidgrass  = calc_prop_focal(Acidgrass, glmRes_Acidgrass$radius[which.max(glmRes_Acidgrass$loglikelihood)])
lcm_Impgrass  = calc_prop_focal(Impgrass, glmRes_Impgrass$radius[which.max(glmRes_Impgrass$loglikelihood)])
lcm_Broadleaved  = calc_prop_focal(Broadleaved, glmRes_Broadleaved$radius[which.max(glmRes_Broadleaved$loglikelihood)])

glmRes_restrict$radius[which.max(glmRes_restrict$loglikelihood)]
glmRes_Acidgrass$radius[which.max(glmRes_Acidgrass$loglikelihood)]
glmRes_Impgrass$radius[which.max(glmRes_Impgrass$loglikelihood)]
glmRes_Broadleaved$radius[which.max(glmRes_Broadleaved$loglikelihood)]

#Arable
dist_raster <- distance(Arable, target = 0)
dist_raster <- mask(dist_raster, Arable)
Arable_dist <- exp(-dist_raster / 5000)
plot(Arable_dist)

demScot=terra::resample(demScot,lcm_restrict)

# 2.1 add all environment variables ----
allEnv=c(lcm_Acidgrass,lcm_Impgrass,lcm_restrict,lcm_Broadleaved,Arable_dist,demScot)
names(allEnv)=c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem')

# 2.2 ‘presence’ and ‘background’ -----

back_sf=st_as_sf(back,crs="EPSG:27700")
Melesmeles_sf = st_as_sf(
  MelesmelesFin,
  coords = c("x", "y"),
  crs="EPSG:27700"
)
## organize ‘presence’ and ‘background’----

# extract env var for presence point
eP=terra::extract(allEnv,Melesmeles_sf)
eB=terra::extract(allEnv,back_sf)
# env var + presence point 
Pres.cov=st_as_sf(cbind(eP,Melesmeles_sf))
Back.cov=st_as_sf(cbind(eB,back_sf))

## add text for differencr point ----

# add text for presence point (MelesmelesFin)
Pres.cov$Pres=1
# delet ID
Pres.cov=Pres.cov[,-1]

# add text for background point
Back.cov$Pres=0
Back.cov = Back.cov[,-1]

# 2.3 merging data ----

## merging point coordinates ----
# get coordinate
coordsPres=st_coordinates(Pres.cov)
coordsBack=st_coordinates(back_sf)

# merging point coordinates
coords=data.frame(rbind(coordsPres,coordsBack))
colnames(coords)=c("x","y")
names(Pres.cov)
names(Back.cov)
all.cov=rbind(Pres.cov,Back.cov)
all.cov=cbind(all.cov,coords)


## data cleaning ----
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)


# 2.4 VIF validation ----

vif_model <- lm(Pres ~ Acidgrass + Impgrass + restrict + Broadleaved + Arable + dem, data = all.cov)
vif(vif_model)
cor(all.cov[, c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem')])
str(all.cov)


# PCA
vars <- all.cov[, c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem')]
pca <- prcomp(vars, scale. = TRUE)

pc_data <- as.data.frame(pca$x)

pca$rotation
summary(pca)

####----------------- 3 define parameters ----------------- ####

# 3.1 Create category task ----
pc_data$Pres <- as.factor(all.cov$Pres)
pc_data$x <- all.cov$x
pc_data$y <- all.cov$y

task = makeClassifTask(
  data = pc_data[, c("PC1", "PC2", "PC3","Pres")],
  target = "Pres",
  positive = "1",
  coordinates = pc_data[, c("x", "y")])

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



####----------------- 4 Binomial ----------------- ####

# 4.1 define Binomial (logistic regression) ----
lrnBinomial = makeLearner(
  "classif.binomial",      #Classification model of binomial distribution
  predict.type = "prob",
  fix.factors.prediction = TRUE)


# 4.2 validation ----
cvBinomial = mlr::resample(
  learner = lrnBinomial,
  task =task,
  resampling = perf_levelCV, 
  measures = mlr::auc,
  show.info = FALSE)

print(cvBinomial)

sp_cvBinomial = mlr::resample(
  learner = lrnBinomial,
  task =task,
  resampling = perf_level_spCV, 
  measures = mlr::auc,
  models = TRUE,
  keep.pred = TRUE)
pred <- sp_cvBinomial$pred
e <- evaluate(pred$data$prob.1[pred$data$truth == 1],
              pred$data$prob.1[pred$data$truth == 0])

print(sp_cvBinomial)
plot(e, "ROC")

sp_cvBinomial = mlr::resample(
  learner = lrnBinomial,
  task =task,
  resampling = perf_level_spCV, 
  measures = mlr::auc,
  models = TRUE,
  keep.pred = TRUE)
pred <- sp_cvBinomial$pred
e <- evaluate(pred$data$prob.1[pred$data$truth == 1],
              pred$data$prob.1[pred$data$truth == 0])

print(sp_cvBinomial)
plot(e, "ROC")

roc_obj <- roc(
  response = pred$data$truth,
  predictor = pred$data$prob.1   # 注意类别名称
)

plot(roc_obj, col = "blue", lwd = 2)

print(sp_cvBinomial)


# 4.3 make partition plots----

plots = createSpatialResamplingPlots(
  task,
  resample=cvBinomial,
  crs=crs(allEnv),
  datum=crs(allEnv),
  color.test = "red",
  point.size = 1)

plotsSP = createSpatialResamplingPlots(
  task,
  resample=sp_cvBinomial,
  crs=crs(allEnv),
  datum=crs(allEnv),
  color.test = "red",
  point.size = 1)

plots = createSpatialResamplingPlots(
  task,
  resample=cvBinomial,
  crs=crs(allEnv),
  datum=crs(allEnv),
  color.test = "red",
  point.size = 1)

library(cowplot)
cowplot::plot_grid(plotlist = plots[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots[["Labels"]])





####----------------- 5 Random Forest ----------------- ####

# 5.1 define Random Forest ----
lrnRF = makeLearner(
  "classif.ranger",
  predict.type = "prob",
  fix.factors.prediction = TRUE)


# 5.2 Random Forest validation----
cvRF = mlr::resample(
  learner = lrnRF,
  task =task,
  resampling = perf_levelCV, 
  measures = mlr::auc,
  show.info = FALSE)

print(cvRF)

sp_cvRF = mlr::resample(
  learner = lrnRF, 
  task =task,
  resampling = perf_level_spCV, 
  measures = mlr::auc,
  show.info = FALSE)

print(sp_cvRF)




# 5.3 Random Forest Parameter Tuning----
getParamSet(lrnRF)
# Parameter tuning command
paramsRF = makeParamSet(
  makeIntegerParam("mtry",lower = 1,upper = 3),
  makeIntegerParam("min.node.size",lower = 1,upper = 20),
  makeIntegerParam("num.trees",lower = 100,upper = 500)
)
# Parameter Tuning Evaluation Method
tune_level = makeResampleDesc(method = "SpCV", iters = 5)
# Parameter search method
ctrl = makeTuneControlRandom(maxit = 50)

tuned_RF = tuneParams(learner = lrnRF,
                      task = task,
                      resampling = tune_level,
                      measures = mlr::auc,
                      par.set = paramsRF,
                      control = ctrl,
                      show.info = FALSE)

print(tuned_RF)



####----------------- 6 MaxEnt ----------------- ####

# 6.1 spacial grid ----
area_grid = st_make_grid(
  MelesmelesFin,
  c(50000, 50000),
  what = "polygons",
  square = T)

area_grid_sf=st_as_sf(area_grid)
area_grid_sf$grid_id=1:length(lengths(area_grid))

# plot and check
plot(area_grid_sf$x)
plot(Melesmeles_sf$geometry,add=T)

# define spacial folds
folds=area_grid_sf$grid_id

# Unified coordinate system
dataPoints=st_as_sf(all.cov,coords = c("x","y"))
st_crs(dataPoints)=crs(area_grid_sf)

# num for fold
folds=5

# presence and background
### why do this ???
Pres.cov=all.cov[all.cov$Pres==1,]
Back.cov=all.cov[all.cov$Pres==0,]

kfold_pres = kfold(Pres.cov, folds)
kfold_back = kfold(Back.cov, folds)


# 6.2 cross-validation ----
# AUC list
eMax=list()

for (i in 1:folds) {
  # get train and test data
  train = Pres.cov[kfold_pres!= i,]
  test = Pres.cov[kfold_pres == i,]
  
  backTrain=Back.cov[kfold_back!=i,]
  backTest=Back.cov[kfold_back==i,]
  
  dataTrain=rbind(train,backTrain)
  dataTest=rbind(test,backTest)
  
  # training
  maxnetMod=maxnet(dataTrain$Pres, dataTrain[,1:4])
  # evaluate on test data（AUC）
  eMax[[i]] = evaluate(p=dataTest[ which(dataTest$Pres==1),],a=dataTest[which(dataTest$Pres==0),],maxnetMod)
}

# get all AUC
aucMax = sapply(eMax, function(x){slot(x, 'auc')} )

print(mean(aucMax))


# 6.3 spatial cross-validation ----

# Use grid IDs as folds for spatial partitioning
folds = area_grid_sf$grid_id

# Initialize list to store AUC results for each fold
maxEvalList = list()

for (i in folds) {
  
  # Define training area (all grids except the current one)
  gridTrain = subset(area_grid_sf, area_grid_sf$grid_id != i)
  
  # Extract training points by spatial intersection with training grids
  # and remove geometry to create a data frame
  train = data.frame(
    st_drop_geometry(
      st_intersection(gridTrain, dataPoints)
    )
  )
  
  # Define test area (current grid only)
  gridTest = subset(area_grid_sf, area_grid_sf$grid_id == i)
  
  # Extract test points within the test grid
  test = data.frame(
    st_drop_geometry(
      st_intersection(gridTest, dataPoints)
    )
  )  
  
  # Train MaxEnt model (maxnet)
  # train$Pres: presence/absence response
  # train[1:4]: predictor variables
  # classes = "lq": linear and quadratic features
  maxnetMod = maxnet(train$Pres, train[1:4],
                     classes = "lq")   
  
  # Predict on test data (cloglog scale gives occurrence probability)
  pred = predict(maxnetMod, test, type = "cloglog")
  
  # Evaluate model performance (ROC and PR curves)
  precrec_proc = evalmod(scores = pred,
                         labels = test$Pres,
                         mode = "prcroc")
  
  # Calculate AUC value
  modauc = precrec::auc(
    precrec::evalmod(scores = pred, 
                     labels = test$Pres)
  )
  
  # Store AUC of current fold
  maxEvalList[[i]] = modauc$aucs[1]
  
  # Print current fold ID (progress tracking)
  print(i)
}

# Compute mean AUC across all spatial folds
mean(unlist(maxEvalList))




####----------------- 7. Ensemble Modeling (GLM + RF + MaxEnt) ----------------- ####

### 7.1 Calculate Model Weights (Based on Spatial CV AUC) ----
# Extract AUC from spatial cross-validation results
auc_values <- c(
  binomial = as.numeric(sp_cvBinomial$aggr),
  rf       = as.numeric(sp_cvRF$aggr),
  maxent   = mean(unlist(maxEvalList))
)

# Normalize weights so they sum to 1
weights <- auc_values / sum(auc_values)
print("Model Weights:")
print(weights)

### 7.2 Train Final Models using Full Dataset ----
# 1. Logistic Regression (GLM)
lrnBinomial <- makeLearner("classif.binomial", predict.type = "prob", fix.factors.prediction = TRUE)
fit_glm <- mlr::train(lrnBinomial, task)

# 2. Random Forest (RF)
lrnRF <- makeLearner("classif.ranger", predict.type = "prob", fix.factors.prediction = TRUE)
if(exists("tuned_RF")) lrnRF <- setHyperPars(lrnRF, par.vals = tuned_RF$x)
fit_rf <- mlr::train(lrnRF, task)

# 3. MaxEnt
pc_data$Pres <- as.numeric(as.character(pc_data$Pres))
maxnet_data <- pc_data[, c("PC1", "PC2", "PC3")]
maxnet_pres <- pc_data$Pres
fit_maxent  <- maxnet(maxnet_pres, maxnet_data, classes = "lq")
unique(pc_data$Pres)
### 7.3 Prepare Environmental Data for Prediction ----

allPC <- terra::predict(allEnv, pca, index = 1:3)
names(allPC) <- c("PC1", "PC2", "PC3")

plot(allPC)

# PC to Dataframe
env_df_all <- as.data.frame(allPC, xy = TRUE, na.rm = FALSE)
pred_vars <- c("PC1", "PC2", "PC3")

final_mask <- complete.cases(env_df_all[, pred_vars])
env_pred <- env_df_all[final_mask, pred_vars]

### 7.4 Generate Predictions ----
# Predict probability for each model
pred_glm    <- predict(fit_glm, newdata = env_pred)$data[, "prob.1"]
pred_rf     <- predict(fit_rf, newdata = env_pred)$data[, "prob.1"]
pred_maxent <- predict(fit_maxent, env_pred, type = "cloglog")

### 7.5 Weighted Ensemble Calculation ----
# Combine predictions using the AUC-based weights
ensemble_vec <- (weights[1] * pred_glm) + 
  (weights[2] * pred_rf) + 
  (weights[3] * pred_maxent)

### 7.6 Map Generation and Visualization ----
# Create an empty raster template based on the original environment
ensemble_raster <- rast(allEnv[[1]])
values(ensemble_raster) <- NA

# Fill the valid cells with ensemble prediction values
ensemble_raster[which(final_mask)] <- ensemble_vec
names(ensemble_raster) <- "Occurrence_Probability"
setMinMax(ensemble_raster)

# Final Plot
plot(ensemble_raster, 
     main = "Ensemble Suitability Map: Meles meles", 
     col = terrain.colors(100))

# Overlay presence points to verify accuracy
points(MelesmelesFin, pch = 20, col = "red", cex = 0.5)


####----------------- 8 Point Process Modelling ----------------- ####

# 8.1 Data preprocessing ----
## Prepare environmental covariates ----

# Convert raster layers to 'im' objects (required by spatstat)
foodIm    = raster.as.im(raster(allEnv$food))     
urbanIm   = raster.as.im(raster(allEnv$urban))
habitatIm = raster.as.im(raster(allEnv$habitat))
demIm     = raster.as.im(raster(allEnv$dem))  


## Define study window ----
# Create study window based on the urban raster layer
window.poly = as.owin(urbanIm)

# Inspect the study window
plot(window.poly)


## Create point pattern object ----

# Extract coordinates from occurrence data (sf object)
MelesmelesCoords = st_coordinates(MelesmelesFin)

# Create point pattern (ppp object)
pppMelesmeles = ppp(MelesmelesCoords[,1],
                    MelesmelesCoords[,2],
                    window = window.poly)

# Visual check: raster layer + point pattern
plot(allEnv$food)
plot(pppMelesmeles, add = TRUE)

# Remove points located outside the study window
pppMelesmeles = as.ppp(pppMelesmeles)


## Rescale spatial units ----

# Convert units from meters to kilometers (to improve numerical stability)
pppMelesmeles = spatstat.geom::rescale(pppMelesmeles, 1000)

# Ensure all covariate rasters are on the same spatial scale
foodIm    = spatstat.geom::rescale(foodIm, 1000)
urbanIm   = spatstat.geom::rescale(urbanIm, 1000)
habitatIm = spatstat.geom::rescale(habitatIm, 1000)
demIm     = spatstat.geom::rescale(demIm, 1000)



# 8.2 Exploratory analysis: test Complete Spatial Randomness (CSR)----

# Compute Ripley’s K-function with simulation envelopes
Kcsr = envelope(pppMelesmeles,
                Kest,
                nsim = 39,
                VARIANCE = TRUE,
                nSD = 1,
                global = TRUE)

# Plot results to assess deviation from CSR
plot(Kcsr, shade = c("hi", "lo"), legend = TRUE)



# 8.3 Select quadrature scheme (background point density)----

# Test different grid densities for quadrature scheme
ndTry = seq(100, 1000, by = 100)

for(i in ndTry){
  
  # Construct quadrature scheme (data points + dummy/background points)
  Q.i = quadscheme(pppMelesmeles,
                   method = "grid",
                   nd = i)
  
  # Fit a simple Poisson point process model
  fit.i = ppm(Q.i ~ foodIm + urbanIm + demIm + habitatIm)
  
  print(i)
  
  # Evaluate model fit using AIC
  print(AIC(fit.i))
}

# Select optimal grid density (based on AIC)
Q = quadscheme(pppMelesmeles, method = "grid", nd = 900)



# 8.4 Explore covariate-response relationships ----

# Non-parametric estimation of intensity vs environmental covariates
plot(rhohat(pppMelesmeles, foodIm))    
plot(rhohat(pppMelesmeles, demIm))     
plot(rhohat(pppMelesmeles, urbanIm))   
plot(rhohat(pppMelesmeles, habitatIm))  


# 8.5 Fit Poisson Point Process Model (PPM)----

# Fit inhomogeneous Poisson point process model
# Polynomial terms are used to capture non-linear relationships
# Spatial coordinates (x, y) account for large-scale spatial trends
firstPPMod = ppm(Q ~ poly(foodIm, 3) +
                   poly(habitatIm, 3) +
                   poly(demIm, 2) +
                   poly(urbanIm, 3) +
                   x + y)

# Model diagnostics using simulation envelope of K-function
firstModEnv = envelope(firstPPMod,
                       Kest,
                       nsim = 39,
                       VARIANCE = TRUE,
                       nSD = 1,
                       global = TRUE)

plot(firstModEnv)



# 8.6 Fit cluster point process model (kppm - Thomas process)----

# Fit Thomas cluster process model to account for spatial aggregation
thomasMod = kppm(Q ~ poly(foodIm) +
                   poly(habitatIm) +
                   poly(urbanIm) +
                   x + y,
                 "Thomas")

# Evaluate model fit using K-function envelope
thomasEnv = envelope(thomasMod,
                     Kest,
                     nsim = 39,
                     VARIANCE = TRUE,
                     nSD = 1,
                     global = TRUE)

plot(thomasEnv)



# 8.7 Model evaluation and spatial prediction ----

# Plot ROC curve
plot(roc(thomasMod))

# Calculate AUC (model discrimination ability)
auc.kppm(thomasMod)

# Predict spatial intensity (occurrence probability surface)
prPPMod = predict(thomasMod)

# Visualize prediction
plot(prPPMod)

# Convert prediction to raster format for GIS visualization
plot(rast(prPPMod))

