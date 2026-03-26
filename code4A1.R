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
  
  pointsBuffer <- st_buffer(speciesData, dist = r)                     
  
  pointsBuffer_vect <- vect(pointsBuffer)
  
  bufferlandcover <- crop(landcover, pointsBuffer)              
  
  masklandcover <- extract(landcover, pointsBuffer_vect, weights = TRUE)      
  
  percentcover <- tapply(1:nrow(masklandcover), masklandcover[,1], function(i){
    
    df <- masklandcover[i,]
    
    cover_area <- sum(df$weight[df[,2] == 1], na.rm = TRUE)
    
    total_area <- sum(df$weight, na.rm = TRUE)
    
    return(cover_area / total_area * 100)
  })
  
  return(percentcover)                                       
}

scale_optimization <- function(resList, radii, Pres){
  
  resFin <- do.call(cbind, resList)
  glmData <- data.frame(resFin)
  colnames(glmData) <- paste0("radius", radii)
  glmData$Pres <- Pres
  
  glmRes <- data.frame(radius = radii, logLik = NA)
  
  for(i in seq_along(radii)){
    
    varname <- paste0("radius", radii[i])
    
    glm.i <- glm(glmData$Pres ~ glmData[[varname]],
                 family = "binomial")
    
    glmRes$logLik[i] <- as.numeric(logLik(glm.i))
  }
  
  opt <- glmRes[which.max(glmRes$logLik), ]
  
  return(list(
    glmRes = glmRes,
    optimal = opt,
    data = glmData
  ))
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

Melesmeles = Melesmeles[Melesmeles$Coordinate.uncertainty_m == 1000, ]

## convert to spatial point
Melesmeles_latlong = data.frame(
  x = Melesmeles$Longitude, 
  y = Melesmeles$Latitude )

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
scot_vect = vect(scot)
Melesmeles_sp = project(Melesmeles_sp, "EPSG:27700")

LCM = rast("./data/LCMUK.tif")
LCM = LCM$LCMUK_1


Melesmeles_crop <- crop(Melesmeles_sp, scot_vect)
LCM <- crop(LCM, scot_vect)
LCM <- mask(LCM, scot_vect)

MelesmelesFin = project(Melesmeles_crop, crs(LCM))


# 1.2 data processing ----
## LCM raster ----
LCM=crop(LCM,st_buffer(scot, dist= 1000))
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")
# add why choose fact=4

## create the background points----
set.seed(33)

back = spatSample(
  LCM,size=1000,
  as.points=TRUE,
  method="random",
  na.rm=TRUE) 

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
unique(LCM)

# restricted area:Mountain, heath, bog
reclassrestrict = c(rep(0,9),rep(1,4),rep(0,9))
restrict = reclass_binary(LCM, reclassrestrict)

# habitat
reclasshabitat = c(0,1,1,rep(0,19))
habitat = reclass_binary(LCM, reclasshabitat)

# food
reclassfood = c(rep(0,3),rep(1,2),rep(0,2),1,rep(0,14))
food = reclass_binary(LCM, reclassfood)

# urban
reclassurban = c(rep(0,20),1,1)
urban = reclass_binary(LCM, reclassurban)

## 2.2 buffer ----

# calculate different buffer
radii<-seq(100,3000,by=200)
resList=list()

# habitat ----
for(i in radii){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,habitat)
  res.i
  resList[[i/100]]=res.i
  print(i)
  gc()
}

result_habit <- scale_optimization(resList, radii, MelesmelesData$Pres)

glmRes_habit <- result_habit$glmRes
opt_habit <- result_habit$optimal

# food ----
for(i in radii){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,food)
  res.i
  resList[[i/100]]=res.i
  print(i)
  gc()
}

result_food <- scale_optimization(resList, radii, MelesmelesData$Pres)

glmRes_food <- result_food$glmRes
opt_food <- result_food$optimal

# restrict ----
for(i in radii){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,restrict)
  res.i
  resList[[i/100]]=res.i
  print(i)
  gc()
}

result_restrict <- scale_optimization(resList, radii, MelesmelesData$Pres)

glmRes_restrict <- result_restrict$glmRes
opt_restrict <- result_restrict$optimal

print(opt_food,opt_habit,opt_restrict)

par(mar = c(5,5,4,2))

# draw habit ----
plot(glmRes_habit$radius, glmRes_habit$logLik,
     type = "b",
     pch = 19,
     lwd = 2,
     col = "#D55E00",
     cex = 1.2,
     cex.axis = 1,
     cex.lab = 1.2,
     bty = "n",
     main = "Habitat (栖息地)",
     xlab = "Buffer size (m)",
     ylab = "Log-likelihood")

grid(col = "lightgray", lty = "dotted")
opt_idx <- which.max(glmRes_habit$logLik)

points(glmRes_habit$radius[opt_idx],
       glmRes_habit$logLik[opt_idx],
       pch = 19,
       col = "blue",
       cex = 1.5)

text(glmRes_habit$radius[opt_idx],
     glmRes_habit$logLik[opt_idx],
     labels = paste0("Opt = ", glmRes_habit$radius[opt_idx]),
     pos = 3,
     col = "blue",
     cex = 0.9)


# draw food ----
plot(glmRes_food$radius, glmRes_food$logLik,
     type = "b",
     pch = 19,
     lwd = 2,
     col = "#D55E00",
     cex = 1.2,
     cex.axis = 1,
     cex.lab = 1.2,
     bty = "n",
     main = "food ()",
     xlab = "Buffer size (m)",
     ylab = "Log-likelihood")

grid(col = "lightgray", lty = "dotted")
opt_idx <- which.max(glmRes_food$logLik)

points(glmRes_food$radius[opt_idx],
       glmRes_food$logLik[opt_idx],
       pch = 19,
       col = "blue",
       cex = 1.5)

text(glmRes_food$radius[opt_idx],
     glmRes_food$logLik[opt_idx],
     labels = paste0("Opt = ", glmRes_food$radius[opt_idx]),
     pos = 3,
     col = "blue",
     cex = 0.9)


# draw restrict ----
plot(glmRes_restrict$radius, glmRes_restrict$logLik,
     type = "b",
     pch = 19,
     lwd = 2,
     col = "#D55E00",
     cex = 1.2,
     cex.axis = 1,
     cex.lab = 1.2,
     bty = "n",
     main = "restrict ()",
     xlab = "Buffer size (m)",
     ylab = "Log-likelihood")

grid(col = "lightgray", lty = "dotted")
opt_idx <- which.max(glmRes_restrict$logLik)

points(glmRes_restrict$radius[opt_idx],
       glmRes_restrict$logLik[opt_idx],
       pch = 19,
       col = "blue",
       cex = 1.5)

text(glmRes_restrict$radius[opt_idx],
     glmRes_restrict$logLik[opt_idx],
     labels = paste0("Opt = ", glmRes_restrict$radius[opt_idx]),
     pos = 3,
     col = "blue",
     cex = 0.9)


#urban
dist_raster <- distance(urban, target = 0)
dist_raster <- mask(dist_raster, urban)
urban_dist <- exp(-dist_raster / 5000)

#Landscape proportion (moving window) ----

lcm_habit = calc_prop_focal(habitat, opt_habit$radius)
lcm_food = calc_prop_focal(food, opt_food$radius)
lcm_restrict  = calc_prop_focal(restrict, opt_restrict$radius)

demScot=terra::resample(demScot,lcm_habit)

# 2.1 add all environment variables ----
allEnv=c(lcm_habit,lcm_food,lcm_restrict,urban_dist,demScot)
names(allEnv)=c('habit','food','restrict','urban','dem')

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
Back.cov = Back.cov[ , -c(1, 7)]

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

vif_model <- lm(Pres ~ habit + food + restrict + urban + dem, data = all.cov)
vif(vif_model)
cor(all.cov[, c("habit", "food", "restrict", "urban", "dem")])
str(all.cov)



# PCA
vars <- all.cov[, c("habit", "food", "restrict", "urban", "dem")]
pca <- prcomp(vars, scale. = TRUE)

pc_data <- as.data.frame(pca$x)

pca$rotation


####----------------- 3 define parameters ----------------- ####

# 3.1 Create category task ----
pc_data$Pres <- as.factor(all.cov$Pres)
pc_data$x <- all.cov$x
pc_data$y <- all.cov$y

task = makeClassifTask(
  data = pc_data[, c("PC1", "PC2", "PC3", "Pres")],
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
  show.info = FALSE)

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
maxnet_data <- pc_data[, c("PC1", "PC2", "PC3")]
maxnet_pres <- pc_data$Pres
fit_maxent  <- maxnet(maxnet_pres, maxnet_data, classes = "lq")

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

# Save the final result
writeRaster(ensemble_raster, "Ensemble_Suitability_Final.tif", overwrite = TRUE)


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

