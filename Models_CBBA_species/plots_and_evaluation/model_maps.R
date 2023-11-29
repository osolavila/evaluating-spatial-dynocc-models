library(foreach)
library(doParallel)

load("./Data/Site_covariables/cat_site_covs.Rdata")


distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5)) + 0.2
cells_id_coord$lon <- as.numeric(cells_id_coord$lon)
cells_id_coord$lat <- as.numeric(cells_id_coord$lat)

cat_aonc2_list <- list()

for (study_species in species.info$spp){
  cat_aonc2 <- get_cat_sdm(study_species)
  cat_aonc2_list[[study_species]] <- cat_aonc2
}



####  BRM cat connectivity best d ####

if (file.exists(file.path(output.path, "brm_cat_connectivity_best_d.Rdata"))){
  load(file.path(output.path, "brm_cat_connectivity_best_d.Rdata"))
} else {
  n.cores <- parallel::detectCores() - 2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  brm_cat_connectivity_best_d <- foreach(
    sp = 1:length(species.info$spp),
    .packages = "raster"
  ) %dopar% {
  
  
    # print(paste("Calculating connectivity:", sp, species.info$spp[sp]) )
    brm_best_d <- which.max(metrics2$auc$abs$occu_BRM_colext[[sp]] )
    max_dist <- distances_analysed[ min( brm_best_d, 8 ) ] ##20km max
  
    ##SpeciesAONC2 SDM (categories) for all Catalonia
    cat_aonc_cat <- as.integer(as.character(
      cat_aonc2_list[[sp]]$aonc2_sdm_cat))
  
    sp.autocovariate <- vector(length=nrow(cells_id_coord))
    for (cell_i in 1:nrow(cells_id_coord)) {
      ##Neighbours
      points_distances <- pointDistance(cells_id_coord[cell_i ,c("lon", "lat")],
                                        cells_id_coord[,c("lon", "lat")],
                                        lonlat=F)
      points_distances <- points_distances/1000
      ##Eliminate same cell distances
      points_distances[cell_i] <- NA
  
      index_max_dist <- which(points_distances<max_dist)
      # a <- which(points_distances<max_dist)
      # b <- points_distances[a]
  
      ##Calculate neighbours info (whithout scaling)
      cell.autocovariate <- sum(cat_aonc_cat[index_max_dist], na.rm=T)/
        length(index_max_dist)
      sp.autocovariate[cell_i] <- cell.autocovariate
  
    }
  
    study.autocovariate <- sp.autocovariate[index_study_in_cat]
    # c <- sapply(index_max_dist[[3]],
    #             function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x))
    # scaled.study.autocovariate <- scale(study.autocovariate)
    # cor(scaled.study.autocovariate, site_covs$sp.autocovariate)
    ## Values are not exactly the same that connectivity calculated for the models
    ## Since coordinates study cells are different grid centroid
    ## But correlation is 0.999 so I should not worry
    mean.autocovariate <- mean(study.autocovariate)
    sd.autocovariate <- sd(study.autocovariate)
  
    scaled.sp.autocovariate <- (sp.autocovariate - mean.autocovariate) /
      sd.autocovariate
  
    return( list(
      con = data.frame( cbind(sp.autocovariate, scaled.sp.autocovariate) ),
      d = max_dist,
      sp = species.info$spp[sp]) )
  
  }
  
  parallel::stopCluster(cl = my.cluster)
  
  names(brm_cat_connectivity_best_d) <- species.info$spp
  
  save(brm_cat_connectivity_best_d,
       file=file.path(output.path, "brm_cat_connectivity_best_d.Rdata"))
}


####  BRM-habitat cat connectivity best d ####


n.cores <- parallel::detectCores() - 2

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

brm_habitat_cat_connectivity_best_d <- foreach(
  sp = 1:length(species.info$spp),
  .packages = "raster"
) %dopar% {
  
  
  # print(paste("Calculating connectivity:", sp, species.info$spp[sp]) )
  brm_best_d <- which.max(metrics2$auc$abs$occu_BRM_habitat_colext[[sp]] )
  max_dist <- distances_analysed[ min( brm_best_d, 8 ) ] ##20km max
  
  ##SpeciesAONC2 SDM (categories) for all Catalonia
  cat_aonc_cat <- as.integer(as.character(
    cat_aonc2_list[[sp]]$aonc2_sdm_cat))
  
  sp.autocovariate <- vector(length=nrow(cells_id_coord))
  for (cell_i in 1:nrow(cells_id_coord)) {
    ##Neighbours
    points_distances <- pointDistance(cells_id_coord[cell_i ,c("lon", "lat")],
                                      cells_id_coord[,c("lon", "lat")],
                                      lonlat=F)
    points_distances <- points_distances/1000
    ##Eliminate same cell distances
    points_distances[cell_i] <- NA
    
    index_max_dist <- which(points_distances<max_dist)
    # a <- which(points_distances<max_dist)
    # b <- points_distances[a]
    
    ##Calculate neighbours info (whithout scaling)
    cell.autocovariate <- sum(cat_aonc_cat[index_max_dist], na.rm=T)/
      length(index_max_dist)
    sp.autocovariate[cell_i] <- cell.autocovariate
    
  }
  
  study.autocovariate <- sp.autocovariate[index_study_in_cat]
  # c <- sapply(index_max_dist[[3]],
  #             function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x))
  # scaled.study.autocovariate <- scale(study.autocovariate)
  # cor(scaled.study.autocovariate, site_covs$sp.autocovariate)
  ## Values are not exactly the same that connectivity calculated for the models
  ## Since coordinates study cells are different grid centroid
  ## But correlation is 0.999 so I should not worry
  mean.autocovariate <- mean(study.autocovariate)
  sd.autocovariate <- sd(study.autocovariate)
  
  scaled.sp.autocovariate <- (sp.autocovariate - mean.autocovariate) /
    sd.autocovariate
  
  return( list(
    con = data.frame( cbind(sp.autocovariate, scaled.sp.autocovariate) ),
    d = max_dist,
    sp = species.info$spp[sp]) )
  
}

parallel::stopCluster(cl = my.cluster)

names(brm_habitat_cat_connectivity_best_d) <- species.info$spp

save(brm_habitat_cat_connectivity_best_d,
     file=file.path(output.path, "brm_habitat_cat_connectivity_best_d.Rdata") )




#### Predict colext models state variables  ####

predict_colext_models <- function(model.fit, site_covs,
                                  psi0_formula = formula(" ~ aonc2_sdm"),
                                  col_formula,
                                  ext_formula){
  
  ##Occu t0 
  psi.params <- coef(model.fit@estimates$psi)
  ##Gamma
  gamma.params <- coef(model.fit@estimates$col)
  ##Epsilon
  epsilon.params <- coef(model.fit@estimates$ext)
  
  ##Create model matrix
  psi0_mf <- model.frame(psi0_formula, site_covs)
  psi0_mm <- model.matrix(psi0_formula, psi0_mf)
  ##Calculate psi0
  psi0 <- plogis(psi0_mm %*% psi.params)
  # psi0 <- predict(model, type="psi")[,1] ##Equivalent for model sites
  
  ##Create model matrix
  col_mf <- model.frame(col_formula, site_covs)
  col_mm <- model.matrix(col_formula, col_mf)
  ##Calculate col
  col <- plogis(col_mm %*% gamma.params)
  # col <- predict(model, type="col")[,1] ##Equivalent for model sites
  
  
  ##Create model matrix
  ext_mf <- model.frame(ext_formula, site_covs)
  ext_mm <- model.matrix(ext_formula, ext_mf)
  
  ##Calculate ext
  ext <- plogis(ext_mm %*% epsilon.params)
  # ext <- predict(model, type="ext")[,1] ##Equivalent for model sites
  
  psi1 <- psi0*(1-ext) + (1-psi0)*col
  # psi1 <- model@projected[2,2,]
  
  return( data.frame(
    psi0=psi0,
    psi1=psi1,
    diff=psi1-psi0,
    p.col=col,
    p.ext=ext
    
  ))
}

#### BRM  ####
load(file.path(output.path, "fits_BRM_colext.Rdata") )
load( file=file.path(output.path, "brm_cat_connectivity_best_d.Rdata") )


brm_predictions_cat <- list()


for (sp in 1:length(species.info$spp) ){
  # print(paste("Calculating connectivity:", sp, species.info$spp[sp]) )
  brm_best_d <- which.max(metrics2$auc$abs$occu_BRM_colext[[sp]] )
  best_model <-  min( brm_best_d, 8 )
  
  
  
  models <- fits_BRM_colext[[sp]][[best_model]]
  ##Site covariables catalonia
  aonc2_sdm <- cat_aonc2_list[[sp]]$aonc2_sdm
  sp.autocovariate <- brm_cat_connectivity_best_d[[sp]][[
    1]]$scaled.sp.autocovariate
  site_covs <- data.frame(aonc2_sdm= aonc2_sdm,
                          sp.autocovariate=sp.autocovariate)
  
  trains_predictions <- list()
  for (i in 1:length(models)){
    model <- models[[i]][[2]]
  
    
    trains_predictions[[i]] <- predict_colext_models(model, site_covs, 
                                                     psi0_formula = formula(" ~ aonc2_sdm"),
                                                     col_formula = formula(" ~ sp.autocovariate"),
                                                     ext_formula = formula(" ~ sp.autocovariate"))
  }
  ##calculate mean
  brm_predictions_cat[[species.info$spp[sp]]] <- 
    Reduce(`+`, trains_predictions) / length(trains_predictions)
  
  
}
#### Habitat  ####
load(file.path(output.path, "fits_habitat.Rdata") )


habitat_predictions_cat <- list()


for (sp in 1:length(species.info$spp) ){
  
  models <- fits_habitat[[sp]]
  ##Site covariables catalonia
  aonc2_sdm <- cat_aonc2_list[[sp]]$aonc2_sdm
  site_covs <- cbind(cat_site_covs, aonc2_sdm)
  
  col_formula <- as.formula(col.equations[sp])
  ext_formula <- as.formula(ext.equations[sp])
  
  trains_predictions <- list()
  for (i in 1:length(models)){
    model <- models[[i]][[2]]
    
    
    trains_predictions[[i]] <- predict_colext_models(model, site_covs, 
                                                     psi0_formula = formula(" ~ aonc2_sdm"),
                                                     col_formula=col_formula,
                                                     ext_formula=ext_formula)
  }
  ##calculate mean
  habitat_predictions_cat[[species.info$spp[sp]]] <- 
    Reduce(`+`, trains_predictions) / length(trains_predictions)
  
  
}
#### BRM-habitat  ####
load(file.path(output.path, "fits_habitat_BRM.Rdata") )
load( file=file.path(output.path, "brm_habitat_cat_connectivity_best_d.Rdata") )


brm_habitat_predictions_cat <- list()


for (sp in 1:length(species.info$spp) ){
  # print(paste("Calculating connectivity:", sp, species.info$spp[sp]) )
  brm_best_d <- which.max(metrics2$auc$abs$occu_BRM_habitat_colext[[sp]] )
  best_model <-  min( brm_best_d, 8 )
  
  models <- fits_habitat_BRM[[sp]][[best_model]]
  ##Site covariables catalonia
  aonc2_sdm <- cat_aonc2_list[[sp]]$aonc2_sdm
  sp.autocovariate <- brm_habitat_cat_connectivity_best_d[[sp]][[
    1]]$scaled.sp.autocovariate
  site_covs <- data.frame(aonc2_sdm= aonc2_sdm,
                          sp.autocovariate=sp.autocovariate)
  site_covs <- cbind(cat_site_covs, site_covs)
  
  col_formula <- as.formula(paste0(col.equations[sp], " + sp.autocovariate"))
  ext_formula <- as.formula(paste0(ext.equations[sp], " + sp.autocovariate"))
  
  
  trains_predictions <- list()
  for (i in 1:length(models)){
    model <- models[[i]][[2]]
    
    
    trains_predictions[[i]] <- predict_colext_models(model, site_covs, 
                                                     psi0_formula = formula(" ~ aonc2_sdm"),
                                                     col_formula=col_formula,
                                                     ext_formula=ext_formula)
  }
  ##calculate mean
  brm_habitat_predictions_cat[[species.info$spp[sp]]] <- 
    Reduce(`+`, trains_predictions) / length(trains_predictions)
  
  
}

####  Get prediction rasters ####
colext_models_predictions_cat <- list(
  BRM = brm_predictions_cat,
  Habitat = habitat_predictions_cat,
  Habitat_BRM = brm_habitat_predictions_cat
)

save(colext_models_predictions_cat,
     file=file.path(output.path,"colext_models_predictions_cat.Rdata") )

load(file=file.path(output.path,"colext_models_predictions_cat.Rdata"))


load("./Data/cells_1x1_id_coord.Rdata")
load("./Data/raster_malla_1km.Rdata")


list_prediction_raster <- list()

for (pred_type in names(colext_models_predictions_cat[[1]][[1]]) ){
  species_list <- list()
  for (sp in species.info$spp){
    
    model_type_list <- list()
    for (model_type in names(colext_models_predictions_cat) ){
      ##Get predictions all Catalonia cells 1x1
      model_prediction <- colext_models_predictions_cat[[
        model_type]][[sp]][[pred_type]] 
      model_prediction <- data.frame(cell_1x1 = cells_id_coord$cell_1x1,
                                     model_prediction)
      model_prediction <- cells_1x1_id_coord %>% 
        left_join(model_prediction, by="cell_1x1")
      
      prediction_raster <- raster_malla_1km
      prediction_raster[] <- model_prediction$model_prediction[raster_malla_1km[]]
      
      model_type_list[[model_type]] <- prediction_raster
    }
    
    species_list[[sp]] <- model_type_list
  }
  
  list_prediction_raster[[pred_type]] <- species_list
}


#### Change maps 3.0 ####
windows(11.2,4.1)
par(mfrow=c(1,3), oma=c(0,0,3,2))
par(mar = c(1, 1, 1, 1))

for (sp in 1:nrow(species.info)){
  
  prediction_raster <- list_prediction_raster[["diff"]][[sp]]
  
  for (i in 1:length(prediction_raster)){
    
    if (i!=3){
      plot(prediction_raster[[i]], col=colors_canvi, zlim=c(-1,1),
           bty="n", box=FALSE, axes=F, legend=F,
           main=names(prediction_raster)[i] )
    } else{
      plot(prediction_raster[[i]], col=colors_canvi, zlim=c(-1,1),
           bty="n", box=FALSE, axes=F, legend=T,
           main=names(prediction_raster)[i] )
    }
    
  }
  
  mtext(substitute(italic(x), list(x=species.info$spp[sp])), side = 3, line = 1, outer = TRUE)
  
  
  savePlot(filename=file.path(plots.path, "change_maps_png", paste0(sp, ".png")),
           type="png")
}