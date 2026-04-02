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
library(ranger)
library(gstat)
library(sp)
library(raster)

set.seed(23)

#### function ####
create_dataset <- function(raster, var_names,presence_pts, background_pts){
  Env =c(raster)
  names(Env)=var_names
  
  # ‘presence’ and ‘background’ 
  back_sf=st_as_sf(background_pts,crs="EPSG:27700")
  presence_sf = st_as_sf(presence_pts,coords = c("x", "y"),crs="EPSG:27700")
  
  # extract env var for presence point ----
  e4P=terra::extract(Env,presence_sf)
  e4B=terra::extract(Env,back_sf)
  
  # env var + presence point 
  Pres.cov=st_as_sf(cbind(e4P,presence_sf))
  Back.cov=st_as_sf(cbind(e4B,back_sf))
  
  # add text for differencr point 
  Pres.cov$Pres=1
  Pres.cov=Pres.cov[,-1]
  
  Back.cov$Pres=0
  Back.cov = Back.cov[,-1]
  
  
  # merging point coordinates
  coordsPres=st_coordinates(Pres.cov)
  coordsBack=st_coordinates(back_sf)
  
  coords=data.frame(rbind(coordsPres,coordsBack))
  colnames(coords)=c("x","y")
  
  all.cov=rbind(Pres.cov,Back.cov)
  all.cov=cbind(all.cov,coords)
  
  
  # data cleaning
  all.cov=na.omit(all.cov)
  all.cov=st_drop_geometry(all.cov)
  
  return(all.cov)}

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

# 1.1 read data

Melesmeles = read.csv("./data/Melesmeles.csv")
scot=st_read('./data/scotSamp.shp')
LCM=rast("./data/LCMUK.tif")
demScot=rast('./data/demScotland.tif')

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


# chlip----
# Melesmeles coord_dataframe ----
Melesmeles_latlong = data.frame(x = Melesmeles$Longitude, 
                                y = Melesmeles$Latitude )
# point object
Melesmeles_sp = vect( Melesmeles_latlong, geom = c("x","y"), crs = "EPSG:4326")

Melesmeles_proj<-project(Melesmeles_sp,crs(LCM))
scot=st_transform(scot,crs(LCM))
studyExtent = vect(st_union(scot))

Melesmeles_sp <-crop(Melesmeles_proj,studyExtent)

plot(Melesmeles_sp)

## LCM raster ----
LCM = LCM$LCMUK_1

studyExtent_buffer <- buffer(studyExtent, width = 5000)
LCM <- crop(LCM, studyExtent_buffer)
LCM=aggregate(LCM,fact=4,fun="modal")


#### -----------------  create dataset for 1 -----------------  ####
back = spatSample(
  studyExtent,
  size=1000,
  method="random")
# draw ----
par(mar = c(5, 5, 3, 5))

plot(LCM, col = hcl.colors(20, "Viridis"))

plot(scot,
     add = TRUE,
     col = adjustcolor("white", alpha.f = 0.3), 
     border = "white",
     lty = 2,
     lwd = 2)

plot(back, add=TRUE,
     col = adjustcolor("grey", alpha.f = 0.9),
     pch = 16,
     cex = 0.4)

plot(Melesmeles_sp, add=TRUE,
     col = "orange",
     pch = 16,
     cex = 0.5)

legend(x = 220000, y = 690000, 
       legend = c("Background", "Presence"),
       col = c(adjustcolor("lightblue", alpha.f = 0.9), "orange"),
       pch = 16,
       pt.cex = c(0.6, 0.8),
       bty = "o",       
       bg = "white")   

# 标题
title("Distribution of Meles meles in Scotland")



# 1.3 choose feature ----

# extract infor
eA_1<-extract(LCM,back)
eP_1<-extract(LCM,Melesmeles_sp)

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
       cex = 0.9,
       y.intersp = 1.5)


#### ----------------- 2 create dataset -----------------  ####

## definite reclass ----
# restrict
reclassrestrict = c(rep(0,9),rep(1,4),rep(0,9))
# Broadleaved
reclassBroadleaved = c(0,1,rep(0,20))
# Coniferouswoodland
reclassConiferous = c(0,0,1,rep(0,19))
# Arable and horticulture
reclassArable = c(rep(0,3),1,rep(0,18))
# Improved grassland
reclassImpgrass = c(rep(0,4),1,rep(0,17))
# Acid grassland
reclassAcidgrass = c(rep(0,7),1,rep(0,14))

## 2.1 Calculate environment variables ----
 
restrict = reclass_binary(LCM, reclassrestrict)
Broadleaved = reclass_binary(LCM, reclassBroadleaved)
Coniferous = reclass_binary(LCM, reclassConiferous)
Arable = reclass_binary(LCM, reclassArable)
Impgrass = reclass_binary(LCM, reclassImpgrass)
Acidgrass = reclass_binary(LCM, reclassAcidgrass)

## 2.2 buffer ----

Abs<-data.frame(crds(back),Pres=0)
Pres<-data.frame(crds(Melesmeles_sp),Pres=1)

# bind the two data frames by row (both dataframes have the same column headings)
MelesmelesData<-rbind(Pres,Abs)
MelesmelesData_sf=st_as_sf(MelesmelesData,coords=c("x","y"),crs="EPSG:27700")


# Broadleaved# ----
radii_Broadleaved<-seq(100,2500,by=100)
resList_Broadleaved=list()

for(i in radii_Broadleaved){
  res.i=landBuffer(speciesData=MelesmelesData_sf,r=i,landcover=Broadleaved)
  resList_Broadleaved[[i/100]]=res.i
  print(i)}
glmRes_Broadleaved <- Test_GLM(
  resList = resList_Broadleaved,
  speciesData = MelesmelesData,
  radii = radii_Broadleaved
)


# Arable and horticulture ----
radii_Arable<-seq(1000,10000,by=1000)
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


# Improved grassland# ----
radii_Impgrass<-seq(100,1800,by=100)
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


# Acid grassland# ----
radii_Acidgrass<-seq(400,3000,by=200)
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

# draw ----
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
ylim_all <- range(c(glmRes_restrict$loglikelihood,
                    glmRes_Acidgrass$loglikelihood,
                    glmRes_Arable$loglikelihood,
                    glmRes_Broadleaved$loglikelihood))

plot_fun <- function(res, title){
  
  res <- res[is.finite(res$loglikelihood), ]
  
  best_idx <- which.max(res$loglikelihood)
  best_r   <- res$radius[best_idx]
  best_ll  <- res$loglikelihood[best_idx]
  
  plot(res$radius, res$loglikelihood,
       type = "b",
       pch = 16,
       col = "#2C7BB6",
       lwd = 2,
       xlab = "Buffer (m)",
       ylab = "logLik",
       main = title)
  
  points(best_r, best_ll, col = "red", pch = 19, cex = 1.2)
  abline(v = best_r, lty = 2, col = "red")
  abline(h = best_ll, lty = 2, col = "grey60")
  
  text(best_r, best_ll,
       labels = paste0(" ", best_r, " m"),
       pos = 4,
       cex = 0.8,
       col = "red")
}

plot_fun(glmRes_restrict, "Restrict")
plot_fun(glmRes_Acidgrass, "Acid grassland")
plot_fun(glmRes_Impgrass, "Improved grassland")
plot_fun(glmRes_Broadleaved, "Broadleaved woodland")

summary(glmRes_restrict$loglikelihood)
#Landscape proportion (moving window) ----

lcm_restrict  = calc_prop_focal(restrict, glmRes_restrict$radius[which.max(glmRes_restrict$loglikelihood)])
lcm_Acidgrass  = calc_prop_focal(Acidgrass, glmRes_Acidgrass$radius[which.max(glmRes_Acidgrass$loglikelihood)])
lcm_Impgrass  = calc_prop_focal(Impgrass, glmRes_Impgrass$radius[which.max(glmRes_Impgrass$loglikelihood)])
lcm_Broadleaved  = calc_prop_focal(Broadleaved, glmRes_Broadleaved$radius[which.max(glmRes_Broadleaved$loglikelihood)])

glmRes_restrict$radius[which.max(glmRes_restrict$loglikelihood)]
glmRes_Acidgrass$radius[which.max(glmRes_Acidgrass$loglikelihood)]
glmRes_Impgrass$radius[which.max(glmRes_Impgrass$loglikelihood)]
glmRes_Broadleaved$radius[which.max(glmRes_Broadleaved$loglikelihood)]

#Arable ----
par(mfrow = c(1,2))

target <- Arable
target[target == 0] <- NA 
dist_raster <- distance(target)
Arable_dist <- exp(-sqrt(dist_raster) / 300)
par(mfrow = c(1,2), mar = c(4,4,3,5))

plot(dist_raster,
     main = "(a) Distance to arable land",
     col = hcl.colors(50, "Viridis"),
     axes = FALSE,
     box = FALSE)

plot(Arable_dist,
     main = "(b) Distance-decay effect",
     col = hcl.colors(50, "Viridis"),
     axes = FALSE,
     box = FALSE)


demScot=terra::resample(demScot,lcm_restrict)

lcm_Acidgrass   <- mask(lcm_Acidgrass, demScot)
lcm_Impgrass    <- mask(lcm_Impgrass, demScot)
lcm_restrict    <- mask(lcm_restrict, demScot)
lcm_Broadleaved <- mask(lcm_Broadleaved, demScot)
Arable_dist     <- mask(Arable_dist, demScot)

par(mfrow = c(3,2), mar = c(3,3,3,4))

scale01 <- function(x){
  (x - min(values(x), na.rm=TRUE)) /
    (max(values(x), na.rm=TRUE) - min(values(x), na.rm=TRUE))
}

vars <- list(
  scale01(lcm_Acidgrass),
  scale01(lcm_Impgrass),
  scale01(lcm_restrict),
  scale01(lcm_Broadleaved),
  scale01(Arable_dist),
  scale01(demScot)
)

titles <- c("(a) Acid grassland",
            "(b) Improved grassland",
            "(c) Study area",
            "(d) Broadleaved woodland",
            "(e) Arable (decay)",
            "(f) Elevation")


for(i in 1:6){
  plot(vars[[i]],
       main = titles[i],
       col = hcl.colors(50, "Viridis"),
       zlim = c(0,1),
       axes = FALSE,
       box = FALSE)
}

# 2.1 add all environment variables ----

lcm_Acidgrass = scale01(lcm_Acidgrass)
lcm_Impgrass = scale01(lcm_Impgrass)
lcm_restrict = scale01(lcm_restrict)
lcm_Broadleaved = scale01(lcm_Broadleaved)
Arable_dist = scale01(Arable_dist)
demScot = scale01(demScot)

ml_data = create_dataset(c(lcm_Acidgrass,lcm_Impgrass,lcm_restrict,lcm_Broadleaved,Arable_dist,demScot),
                         c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem'),
                         Melesmeles_sp,back)

# 2.4 VIF validation ----

vif_model <- lm(Pres ~ Acidgrass + Impgrass + restrict + Broadleaved + Arable + dem, data = ml_data)
vif(vif_model)

# draw ----
dev.off()

par(mar = c(5,5,3,2))

cor_mat <- cor(ml_data[, c('Acidgrass','Impgrass','restrict',
                           'Broadleaved','Arable','dem')],
               use = "complete.obs")

col_fun <- colorRampPalette(c("#3B4CC0", "white", "#B40426"))

n <- ncol(cor_mat)

image(1:n, 1:n, cor_mat,
      col = col_fun(100),
      axes = FALSE,
      xlab = "", ylab = "",
      zlim = c(-1,1))

axis(1, at = 1:n, labels = FALSE)
text(1:n, par("usr")[3] - 0.2,  
     labels = colnames(cor_mat),
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)

axis(2, at = 1:n,
     labels = colnames(cor_mat),
     las = 2, cex.axis = 0.8)

abline(h = 1:(n+1)-0.5, col = "grey90")
abline(v = 1:(n+1)-0.5, col = "grey90")

for(i in 1:n){
  for(j in 1:n){
    val <- round(cor_mat[i,j],2)
    
    text(i, j, val,
         col = ifelse(abs(val) > 0.5, "white", "black"),
         cex = 0.8)
  }
}

title("Correlation matrix", cex.main = 1.1)


# +++++++++++++++++++ +++++++++++++++++++  +++++++++++++++++++ # +++++++++++++++++++# +++++++++++++++++++
# spacial grid ----
# make a grid for spatial blocking
scot_union <- st_union(scot)

grid_sf <- st_as_sf(st_make_grid(scot_union,
                                 #what = "polygons", square = T,
                                 n = c(2,3)))
grid_sf$grid_id <- 1:nrow(grid_sf)

area_grid_sf <- st_intersection(grid_sf,st_geometry(scot_union))

area_grid_sf <- st_as_sf(area_grid_sf)
area_grid_sf$grid_id <- 1:nrow(area_grid_sf)

# plot----
plot(st_geometry(area_grid_sf),
     col = NA,
     border = "grey40",
     lwd = 1)
centroids = st_centroid(area_grid_sf)

text(st_coordinates(centroids),labels = area_grid_sf$grid_id,
     cex = 1.5,
     col = "red")

plot(Melesmeles_sp,add = TRUE,pch = 1,         
     col = "black",
     cex = 0.7)

dataPoints=st_as_sf(ml_data,coords = c("x","y"))
st_crs(dataPoints)=crs(area_grid_sf)


folds = area_grid_sf$grid_id
vars  = c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem')

# maxtent_data----
library(precrec)
maxEvalList=list()
mm_list = list() 
for (i in folds) {
  
  
  gridTrain=subset(area_grid_sf,area_grid_sf$grid_id!=i)
  
  train=data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  

  gridTest=subset(area_grid_sf,area_grid_sf$grid_id==i)
  
  
  test=data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))  
  
  
  maxnetMod=maxnet(train$Pres, train[,vars],
                   classes="lq")   
  
  
  pred=predict(maxnetMod, test,type="cloglog")
  mm = evalmod(scores = pred, labels = test$Pres)
  mm_list[[i]] = mm
  precrec_proc=evalmod(scores = pred,labels = test$Pres,mode = "prcroc")
  
  modauc=precrec::auc(precrec::evalmod(scores = pred, 
                                       labels = test$Pres))
  
  maxEvalList[[i]]=modauc$aucs[1]
  print(i)
  
}
maxEvalList
mean(unlist(maxEvalList))

par(mfrow = c(2,3),
    mar = c(4,4,3,1),
    bg = "white")

for (i in seq_along(mm_list)) {
  
  mm <- mm_list[[i]]
  auc_val <- precrec::auc(mm)$aucs[1]
  
  x <- mm$roc[[1]]$x
  y <- mm$roc[[1]]$y
  
  plot(0, 0,
       type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste0("AUC= ", round(auc_val, 3)),
       panel.first = {
         rect(0, 0, 1, 1, col = "white", border = NA)
         abline(0, 1, col = "black", lwd = 1)
       })
  
  # 红点
  points(x, y, col = "red", pch = 1)
}

# glm ----
glmEvalList = list()
mg_list = list() 
for (i in folds) {
  
  gridTrain = subset(area_grid_sf, area_grid_sf$grid_id != i)
  train = data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  
  gridTest = subset(area_grid_sf, area_grid_sf$grid_id == i)
  test = data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))  
  
  train$Pres = as.numeric(train$Pres)
  test$Pres  = as.numeric(test$Pres)
  
  form = as.formula(
    paste("Pres ~", paste(vars, collapse = " + "))
  )
  
  glmMod = glm(
    form,
    data = train,
    family = binomial(link = "logit")
  )
  
  pred = predict(glmMod, test[, vars], type = "response")
  mg = evalmod(scores = pred, labels = test$Pres)
  mg_list[[i]] = mg
  precrec_proc = evalmod(scores = pred, labels = test$Pres, mode = "prcroc")
  
  modauc = precrec::auc(
    precrec::evalmod(scores = pred, labels = test$Pres)
  )
  
  glmEvalList[[i]] = modauc$aucs[1]
  
  print(i)
}

glmEvalList
mean(unlist(glmEvalList))

library(ranger)
par(mfrow = c(2,3),
    mar = c(4,4,3,1),
    bg = "white")

for (i in seq_along(mg_list)) {
  
  mg <- mg_list[[i]]
  auc_val <- precrec::auc(mg)$aucs[1]
  
  x <- mm$roc[[1]]$x
  y <- mm$roc[[1]]$y
  
  plot(0, 0,
       type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste0("AUC= ", round(auc_val, 3)),
       panel.first = {
         rect(0, 0, 1, 1, col = "white", border = NA)
         abline(0, 1, col = "black", lwd = 1)
       })
  
  points(x, y, col = "red", pch = 1)
}

# rf ----
rfEvalList = list()
rf_list = list() 
for (i in folds) {
  
  gridTrain = subset(area_grid_sf, area_grid_sf$grid_id != i)
  train = data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  
  gridTest = subset(area_grid_sf, area_grid_sf$grid_id == i)
  test = data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))  
  
  train$Pres = as.factor(train$Pres)
  test$Pres  = as.factor(test$Pres)
  
  mtry_grid = c(2, 3, 4, length(vars))
  min_node_grid = c(1, 5, 10)
  
  best_auc = -Inf
  best_model = NULL
  
  # ---单grid search ---
  for (m in mtry_grid) {
    for (n in min_node_grid) {
      
      rfMod = ranger(
        formula = as.formula(
          paste("Pres ~", paste(vars, collapse = " + "))
        ),
        data = train,
        probability = TRUE,
        mtry = m,
        min.node.size = n,
        num.trees = 500,
        importance = "impurity"
      )
      
      pred_train = predict(rfMod, train[, vars])$predictions[,2]
      
      auc_tmp = precrec::auc(
        precrec::evalmod(scores = pred_train, labels = as.numeric(as.character(train$Pres)))
      )$aucs[1]
      
      if (auc_tmp > best_auc) {
        best_auc = auc_tmp
        best_model = rfMod
      }
    }
  }
  
  pred = predict(best_model, test[, vars])$predictions[,2]
  rf = evalmod(scores = pred, labels = test$Pres)
  rf_list[[i]] = rf
  test_num = as.numeric(as.character(test$Pres))
  
  precrec_proc = evalmod(scores = pred, labels = test_num, mode = "prcroc")
  
  modauc = precrec::auc(
    precrec::evalmod(scores = pred, labels = test_num)
  )
  
  rfEvalList[[i]] = modauc$aucs[1]
  
  print(paste("Fold:", i, "AUC:", modauc$aucs[1]))
}

rfEvalList
mean(unlist(rfEvalList))

par(mfrow = c(2,3),
    mar = c(4,4,3,1),
    bg = "white")

for (i in seq_along(rf_list)) {
  
  rf <- rf_list[[i]]
  auc_val <- precrec::auc(rf)$aucs[1]
  
  x <- rf$roc[[1]]$x
  y <- rf$roc[[1]]$y
  
  plot(0, 0,
       type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste0("AUC= ", round(auc_val, 3)),
       panel.first = {
         rect(0, 0, 1, 1, col = "white", border = NA)
         abline(0, 1, col = "black", lwd = 1)
       })
  
  points(x, y, col = "red", pch = 1)
}


# glmpca ----
glmpcaEvalList = list()
glmpca_list = list() 
for (i in folds) {
  
  gridTrain = subset(area_grid_sf, area_grid_sf$grid_id != i)
  train = data.frame(st_drop_geometry(st_intersection(gridTrain, dataPoints)))
  
  gridTest = subset(area_grid_sf, area_grid_sf$grid_id == i)
  test = data.frame(st_drop_geometry(st_intersection(gridTest, dataPoints)))  
  
  train$Pres = as.numeric(train$Pres)
  test$Pres  = as.numeric(test$Pres)
  
  pca_mod = prcomp(train[, vars], scale. = TRUE)
  
  var_exp = summary(pca_mod)$importance[2,]
  
  k = which(cumsum(var_exp) >= 0.9)[1]
  
  train_pca = as.data.frame(pca_mod$x[, 1:k])
  
  test_pca = as.data.frame(
    predict(pca_mod, newdata = test[, vars])[, 1:k]
  )
  
  train_pca$Pres = train$Pres
  test_pca$Pres  = test$Pres
  
  form = as.formula(
    paste("Pres ~", paste(colnames(train_pca)[-ncol(train_pca)], collapse = "+"))
  )
  
  glmMod = glm(
    form,
    data = train_pca,
    family = binomial(link = "logit")
  )
  
  pred = predict(glmMod, test_pca, type = "response")
  glmpca = evalmod(scores = pred, labels = test$Pres)
  glmpca_list[[i]] = glmpca
  precrec_proc = evalmod(scores = pred, labels = test_pca$Pres, mode = "prcroc")
  
  modauc = precrec::auc(
    precrec::evalmod(scores = pred, labels = test_pca$Pres)
  )
  
  glmpcaEvalList[[i]] = modauc$aucs[1]
  
  print(paste("Fold:", i, "PCs:", k, "AUC:", modauc$aucs[1]))
}

glmpcaEvalList
mean(unlist(glmpcaEvalList))
par(mfrow = c(2,3),
    mar = c(4,4,3,1),
    bg = "white")

for (i in seq_along(glmpca_list)) {
  
  glmpca <- glmpca_list[[i]]
  auc_val <- precrec::auc(glmpca)$aucs[1]
  
  x <- glmpca$roc[[1]]$x
  y <- glmpca$roc[[1]]$y
  
  plot(0, 0,
       type = "n",
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste0("AUC= ", round(auc_val, 3)),
       panel.first = {
         rect(0, 0, 1, 1, col = "white", border = NA)
         abline(0, 1, col = "black", lwd = 1)
       })
  
  points(x, y, col = "red", pch = 1)
}




####----------------- 7. Ensemble Modeling (GLM + RF + MaxEnt) ----------------- ####

### 7.1 Calculate Model Weights (Based on Spatial CV AUC) ----
# Extract AUC from spatial cross-validation results
auc_values <- c(
  binomial = mean(unlist(glmEvalList)),
  rf       = mean(unlist(rfEvalList)),
  maxent   = mean(unlist(maxEvalList)))

weights <- auc_values / sum(auc_values)
form
names(ml_data)
vars  = c('Acidgrass','Impgrass','restrict','Broadleaved','Arable','dem')
form <- as.formula(
  paste("Pres ~", paste(vars, collapse = " + "))
)

fit_glm    <- glm(form, data = ml_data, family = binomial)
fit_rf     <- best_model
fit_maxent <- maxnet(ml_data$Pres, ml_data[,vars], classes = "lq")
allEnv <- c(
  lcm_Acidgrass,
  lcm_Impgrass,
  lcm_restrict,
  lcm_Broadleaved,
  Arable_dist,
  demScot
)
names(allEnv) <- vars
final_mask <- complete.cases(terra::values(allEnv))
env_pred <- as.data.frame(terra::values(allEnv))
env_pred <- env_pred[final_mask, ]
env_pred <- env_pred[, vars]
# 预测

pred_glm    <- predict(fit_glm, newdata = env_pred, type = "response")
pred_rf <- predict(fit_rf, data = env_pred)$predictions[,2]
pred_maxent <- predict(fit_maxent, env_pred, type = "cloglog")
env_pred <- as.data.frame(terra::values(allEnv))
# ensemble
ensemble_vec <- weights["binomial"] * pred_glm +
  weights["rf"]       * pred_rf +
  weights["maxent"]   * pred_maxent

### 7.6 Map Generation and Visualization ----
dev.off
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
points(Melesmeles_sp, pch = 20, col = "red", cex = 0.6)



