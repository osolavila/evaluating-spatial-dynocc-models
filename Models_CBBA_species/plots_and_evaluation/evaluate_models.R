#### Load detection data and analysis functions ####
Dir.Base <- "~/evaluating-spatial-dynamic-occupancy-models"
setwd(Dir.Base)

##Load detection data
source(file.path(Dir.Base, "Paper_scripts", "data_preparation",
                 "Prepare_species_detecion_data.R"))
##Load functions for the analysis
source(file.path(Dir.Base, "Paper_scripts", "Functions",
                 "get_colext_umf.R"))
source(file.path(Dir.Base, "Paper_scripts", "data_preparation",
                 "cat_site_covs.R"))
source(file.path(Dir.Base, "Paper_scripts", "Functions",
                 "data_analysis_utils.R"))

output.path <- file.path(Dir.Base, "colext", "Analysis", "Rdata",
                         "revision", "cross_validation")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}

plots.path <- file.path(Dir.Base, "colext", "Analysis","plots", "revision",
                        "cross_validation")
if (!file.exists(plots.path)){
  dir.create(file.path(plots.path))
}

library(Metrics) ##AUC and rmse
library(reshape2) ##Melt
library(viridis)## Color palette
library(IFMcpp)
library(ParallelSpatialColext)
library(writexl)


#### Gather species info and select study species ####

## DF: Species|#AONC2|#Col|#Ext|%change|spatial_Group
species.occurance.change <- db_aonc_change %>% 
  filter(cell_1x1 %in% cells_study$cell_1x1,
         spp %in% species.groups$spp) %>%
  group_by(spp, .drop=FALSE) %>%
  summarise( aonc2 = length(which(Canvi==0)) + length(which(Canvi==-1)),
             col= length(which(Canvi==1)),
             ext= length(which(Canvi==-1)),
  ) %>%
  # filter(spp %in% species.groups$spp) %>%
  mutate(perc.change= round((col+aonc2-ext)/(aonc2)*100-100,1),
         group= species.groups$group) %>% 
  filter(aonc2>100)

##DF Detection info: detection rates and best model
load(file="./colext/Analysis/Rdata/species_detection_models_info.Rdata") 

##Join species occurrance DF and species detection DF
species.info <- cbind(species.occurance.change, 
                      species_detection_models_info[,-1])


##Species UMF -> Covariates and detection data
species_UMFs<- list()
##SpeciesAONC2 SDM (categories) for all Catalonia
cat_aonc2_list <- list()

for (study_species in species.info$spp){
  simUMF <- get_UMF(study_species, LiDAR = F, habitat_aonc3 = F,
                    SDM = T)
  species_UMFs[[study_species]] <- simUMF
  
  cat_aonc2 <- get_cat_sdm(study_species)
  cat_aonc2_list[[study_species]] <- cat_aonc2
}

study_species_vector <- species.info$spp
##Get detectpion data
species_detections <- list()
for (study_species in study_species_vector){
  simUMF <- species_UMFs[[study_species]]
  
  
  detection_data <- array(simUMF@y, dim = c(nrow(simUMF@y),2,2) )
  ysum <- apply(detection_data, c(1,3), sum)
  ysum[ysum!=0] <- 1 
  
  species_detections[[study_species]] <- ysum
}

#### Gather neighbours distances info ####
##Neighbours
points_distances <- pointDistance(cells_study[,c("longitude", "latitude")],
                                  cells_id_coord[,c("lon", "lat")],
                                  lonlat=F)
points_distances <- points_distances/1000

index_study_in_cat <- match(cells_study$cell_1x1, cells_id_coord$cell_1x1)
##Eliminate same cell distances
for (i in 1:nrow(points_distances)){
  points_distances[i, index_study_in_cat[i]] <- NA
}

index_cells_30km <- list()
distance_cells_30km <- list()
for (i in 1:nrow(cells_study)){
  neigh <- which(points_distances[i, ] < 30)
  index_cells_30km[[i]] <- neigh
  distance_cells_30km[[i]] <- points_distances[i, neigh]
}

distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
index_max_dist <- list()
names_dist <- paste0("max_dist_", c(1:3, seq(5,10,2.5), seq(15,30,5)) )
i <- 1
for (max_dist in c(1:3, seq(5,10,2.5), seq(15,30,5))+ 0.2 ){
  index_max_dist[[names_dist[i]]] <- apply(points_distances,1,function(x)
    which(x<max_dist))
  i<-i+1
}

rm(points_distances)
gc()
#### Habitat formulas ####
##Best habitat variables
variables.path <- file.path(Dir.Base, "colext", "Analysis", "Rdata",
                            "revision", "more_species")
load(file.path(variables.path, "best_habitat_variables_col_models.Rdata"))
load(file.path(variables.path, "best_habitat_variables_ext_models.Rdata"))

col.equations <- vector(length = length(species.info$spp))
names(col.equations) <- species.info$spp

ext.equations <- vector(length = length(species.info$spp))
names(ext.equations) <- species.info$spp

for (study_species in species.info$spp){
  ##Col formula
  col_predictors <- best_habitat_variables_col_models[[study_species]]$variable
  quadratic <- best_habitat_variables_col_models[[study_species]]$Quadratic_term
  quadratic <- quadratic[!is.na(col_predictors)]
  col_predictors <- col_predictors[!is.na(col_predictors)]
  col_predictors[quadratic] <- sapply(col_predictors[quadratic], function (x)
    paste0(x , " + I(", x, "^2)"))
  col_predictors <- paste0(col_predictors, collapse = " + ")
  col_formula <- paste0("~", col_predictors)
  col.equations[study_species] <- col_formula
  
  ##Ext formula
  ext_predictors <- best_habitat_variables_ext_models[[study_species]]$variable
  quadratic <- best_habitat_variables_ext_models[[study_species]]$Quadratic_term
  quadratic <- quadratic[!is.na(ext_predictors)]
  ext_predictors <- ext_predictors[!is.na(ext_predictors)]
  ext_predictors[quadratic] <- sapply(ext_predictors[quadratic], function (x)
    paste0(x , " + I(", x, "^2)"))
  ext_predictors <- paste0(ext_predictors, collapse = " + ")
  ext_formula <- paste0("~", ext_predictors)
  ext.equations[study_species] <- ext_formula
}

study_species <- study_species_vector[1]
simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
                  SDM = F)
site_covs <- simUMF@siteCovs
#### Predict functions  ####
predict_det_colext_model <- function(model.fit, study_species,
                                     col_formula=~1, ext_formula=~1,
                                     spatial.model=NULL,
                                     spatial.col=F, spatial.ext=F,
                                     autolog=NULL,
                                     alpha=NULL, div=NULL){
  
  if (class(model.fit)=="unmarkedFitColExt"){
    ##Occu t0 
    psi.params <- coef(model.fit, type="psi")
    ##Gamma
    gamma.params <- coef(model.fit, type="col")
    ##Epsilon
    epsilon.params <- coef(model.fit, type="ext")
    
  } else if (class(model.fit)=="unmarkedEstimateList"){
    ##Occu t0 
    psi.params <- coef(model.fit@estimates$psi)
    ##Gamma
    gamma.params <- coef(model.fit@estimates$col)
    ##Epsilon
    epsilon.params <- coef(model.fit@estimates$ext)
    
  }
  
  z <- matrix(nrow = nrow(cells_study), ncol=2)
  ##Occu time 0
  
  ##Prob occu t1
  psi <- species_UMFs[[study_species]]@siteCovs$aonc2_sdm
  psi.landscape <- cat_aonc2_list[[study_species]]$aonc2_sdm
  z[,1] <- plogis(cbind(1,psi) %*% psi.params)
  
  ##Create model matrix
  col_mf <- model.frame(col_formula, site_covs)
  col_mm <- model.matrix(col_formula, col_mf)
  
  ##Create model matrix
  ext_mf <- model.frame(ext_formula, site_covs)
  ext_mm <- model.matrix(ext_formula, ext_mf)
  
  if (!spatial.col){
    
    ##Calculate col
    prob.col <- plogis(col_mm %*% gamma.params)
    
  } else if (spatial.model=="IFM"){
    
    ##Calculate col
    prob.col <- plogis(col_mm %*% gamma.params)
    
    base_delta <- prob.col
    
    # prob.col <- rcpp_phi(distance_cells_30km, index_cells_30km,
    #                      base_delta, psi.landscape, alpha, div)
    
    exp_disp_times_aonc <- list()
    for (i in 1:length(distance_cells_30km)){
      exp_times_aonc <-
        exp(-distance_cells_30km[[i]]/alpha)*
        psi.landscape[index_cells_30km[[i]]]/div
      exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
    }
    
    prob.col <- rcpp_y(exp_disp_times_aonc, base_delta)
    
  } else if (spatial.model=="BRM"){
    
    col_mm <- cbind(col_mm, autolog)
    prob.col <- plogis(col_mm %*% gamma.params)
  }
  
  
  ##Ext probability
  if (!spatial.ext){
    ##Ext probability
    prob.ext <- plogis(ext_mm %*% epsilon.params)
  } else if (spatial.model=="IFM"){
    
    ##Calculate ext
    prob.ext <- plogis(ext_mm %*% epsilon.params)
    
    ##Ext probability
    base_delta <- prob.ext
    
    exp_disp_times_aonc <- list()
    for (i in 1:length(distance_cells_30km)){
      exp_times_aonc <-
        exp(-distance_cells_30km[[i]]/alpha)*
        psi.landscape[index_cells_30km[[i]]]/div
      exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
    }
    
    prob.ext <- 1 - rcpp_y(exp_disp_times_aonc, base_delta)
    
  } else if (spatial.model=="BRM"){
    
    ext_mm <- cbind(ext_mm, autolog)
    prob.ext <- plogis(ext_mm %*% epsilon.params)
  }
  
  
  z[,2] <- z[,1]*(1-prob.ext)+(1-z[,1])*prob.col
  
  
  return(z)
  
}


predict.deterministic.base <- function(params_models){
  
  predictions <- list()
  for (study_species in names(params_models)){
    occu_trains <- list()
    for (train in 1:5){
      ##extract params
      model.fit <- params_models[[study_species]][[train]]
      
      z <- predict_det_colext_model(model.fit, study_species)
      occu_trains[[train]] <- z
      
    }
    
    predictions[[study_species]] <- occu_trains
  }
  
  return(predictions)
  
}


predict.deterministic.habitat <- function(params_models){
  
  predictions_habitat <- list()
  for (study_species in names(params_models)){
    col_formula <- as.formula(col.equations[study_species])
    ext_formula <- as.formula(ext.equations[study_species])
    occu_trains <- list()
    for (train in 1:5){
      ##extract params
      model.fit <- params_models[[study_species]][[train]]
      
      z <- predict_det_colext_model(model.fit, study_species,
                                    col_formula=col_formula,
                                    ext_formula=ext_formula)
      occu_trains[[train]] <- z
      
    }
    
    predictions_habitat[[study_species]] <- occu_trains
  }
  
  return(predictions_habitat)
  
}


predict.deterministic.brm <- function(params_models, spatial.col, spatial.ext){
  
  predictions_BRM <- list()
  for (study_species in names(params_models)){
    cat_aonc_cat <- as.integer(as.character(
      cat_aonc2_list[[study_species]]$aonc2_sdm_cat))
    predictions <- list()
    
    for (dist in names(index_max_dist)) {
      
      ##Calculate neighbours info
      sp.autocovariate <- scale( sapply(index_max_dist[[dist]],
                                        function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
      # cat(paste0("Analysing ", study_species, ", dist: ", alpha, "\n"), 
      #     file="log_3.txt", append=TRUE)
      occu_trains <- list()
      for (train in 1:5){
        ##extract params
        model.fit <- params_models[[study_species]][[dist]][[train]]
        
        z <- predict_det_colext_model(model.fit, study_species,
                                      col_formula=~1, ext_formula=~1,
                                      spatial.model= "BRM",
                                      spatial.col= spatial.col,
                                      spatial.ext= spatial.ext,
                                      autolog=sp.autocovariate)
        occu_trains[[train]] <- z
        
      }
      predictions[[dist]] <- occu_trains
      
    }
    
    predictions_BRM[[study_species]] <- predictions
  }
  return(predictions_BRM)
  
}


predict.deterministic.brm.habitat <- function(params_models, spatial.col, spatial.ext){
  
  predictions_BRM <- list()
  for (study_species in names(params_models)){
    col_formula <- as.formula(col.equations[study_species])
    ext_formula <- as.formula(ext.equations[study_species])
    cat_aonc_cat <- as.integer(as.character(
      cat_aonc2_list[[study_species]]$aonc2_sdm_cat))
    predictions <- list()
    
    for (dist in names(index_max_dist)) {
      
      ##Calculate neighbours info
      sp.autocovariate <- scale( sapply(index_max_dist[[dist]],
                                        function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
      # cat(paste0("Analysing ", study_species, ", dist: ", alpha, "\n"), 
      #     file="log_3.txt", append=TRUE)
      occu_trains <- list()
      for (train in 1:5){
        ##extract params
        model.fit <- params_models[[study_species]][[dist]][[train]]
        
        z <- predict_det_colext_model(model.fit, study_species,
                                      col_formula=col_formula,
                                      ext_formula=ext_formula,
                                      spatial.model= "BRM",
                                      spatial.col= spatial.col,
                                      spatial.ext= spatial.ext,
                                      autolog=sp.autocovariate)
        occu_trains[[train]] <- z
        
      }
      predictions[[dist]] <- occu_trains
      
    }
    
    predictions_BRM[[study_species]] <- predictions
  }
  return(predictions_BRM)
  
}

predict.deterministic.ifm <- function(params_models, spatial.col, spatial.ext,
                                      best_div_factors, div_factors){
  
  predictions_IFM <- list()
  for (study_species in names(params_models)){
    
    predictions <- list()
    
    for (alpha in c(1,2.5,5,10,20) ){
      # cat(paste0("Analysing ", study_species, ", dist: ", alpha, "\n"), 
      #     file="log_3.txt", append=TRUE)
      occu_trains <- list()
      for (train in 1:5){
        ##extract params
        model.fit <- params_models[[study_species]][[as.character(alpha)]][[train]]
        div <- div_factors[
          best_div_factors[[
            study_species]][[as.character(alpha)]][[train]] ]
        
        z <- predict_det_colext_model(model.fit, study_species,
                                      col_formula=~1, ext_formula=~1,
                                      spatial.model= "IFM",
                                      spatial.col= spatial.col,
                                      spatial.ext= spatial.ext,
                                      autolog=NULL,
                                      alpha=alpha, div=div)
        occu_trains[[train]] <- z
        
      }
      predictions[[as.character(alpha)]] <- occu_trains
      
    }
    
    predictions_IFM[[study_species]] <- predictions
  }
  return(predictions_IFM)
  
}




predict.deterministic.ifm.habitat <- function(params_models,
                                              best_div_factors, div_factors,
                                              spatial.col=T, spatial.ext=T){
  
  
  predictions_IFM <- list()
  for (study_species in names(params_models)){
    
    predictions <- list()
    col_formula <- as.formula(col.equations[study_species])
    ext_formula <- as.formula(ext.equations[study_species])
    
    for (alpha in c(1,2.5,5,10,20) ){
      # cat(paste0("Analysing ", study_species, ", dist: ", alpha, "\n"), 
      #     file="log_3.txt", append=TRUE)
      occu_trains <- list()
      for (train in 1:5){
        ##extract params
        model.fit <- params_models[[study_species]][[as.character(alpha)]][[train]]
        div <- div_factors[
          best_div_factors[[
            study_species]][[as.character(alpha)]][[train]] ]
        
        z <- predict_det_colext_model(model.fit, study_species,
                                      col_formula=col_formula, ext_formula=ext_formula,
                                      spatial.model= "IFM",
                                      spatial.col= spatial.col,
                                      spatial.ext= spatial.ext,
                                      autolog=NULL,
                                      alpha=alpha, div=div)
        occu_trains[[train]] <- z
        
      }
      predictions[[as.character(alpha)]] <- occu_trains
      
    }
    
    predictions_IFM[[study_species]] <- predictions
  }
  return(predictions_IFM)
  
}




#### Models occu predictions ####
load(file.path(output.path, "train_test_splits.Rdata"))

if (file.exists(file.path(output.path,
                          "models_det_occu_predictions.Rdata"))){
  load(file.path(output.path,
                 "models_det_occu_predictions.Rdata"))
  rm(distance_cells_30km, index_cells_30km, cat_aonc2_list, cells_1x1)
  gc()
} else{
  ## Base models
  load(file.path(output.path, "fits_base.Rdata") )
  estimates_base <- lapply(fits_base, function(x)lapply(x, function(y) y[[2]]))
  base_colext <- predict.deterministic.base(estimates_base)
  range(base_colext$`Aegithalos caudatus`[[4]][-test_ids[[4]],] - fits_base$`Aegithalos caudatus`[[4]][[1]])
  rm(fits_base)

  ## Habitat models
  load(file.path(output.path, "fits_habitat.Rdata") )
  estimates_habitat <- lapply(fits_habitat, function(x)lapply(x, function(y) y[[2]]))
  habitat_colext <- predict.deterministic.habitat(estimates_habitat)
  range(habitat_colext$`Aegithalos caudatus`[[4]][-test_ids[[4]],] - fits_habitat$`Aegithalos caudatus`[[4]][[1]])
  rm(fits_habitat)
  
  ## BRM models
  load(file.path(output.path, "fits_BRM_col.Rdata") )
  estimates_brm_col <- lapply(fits_BRM_col, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  brm_col <- predict.deterministic.brm(estimates_brm_col, T, F)
  range(brm_col$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_BRM_col$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_BRM_col)
  
  load(file.path(output.path, "fits_BRM_ext.Rdata") )
  estimates_brm_ext <- lapply(fits_BRM_ext, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  brm_ext <- predict.deterministic.brm(estimates_brm_ext, F, T)
  range(brm_ext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_BRM_ext$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_BRM_ext)
  
  load(file.path(output.path, "fits_BRM_colext.Rdata") )
  estimates_brm_colext <- lapply(fits_BRM_colext, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  brm_colext <- predict.deterministic.brm(estimates_brm_colext, T, T)
  range(brm_colext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_BRM_colext$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_BRM_colext)
  
  load(file.path(output.path, "fits_habitat_BRM.Rdata") )
  estimates_habitat_brm_colext <- lapply(fits_habitat_BRM, function(x)
    lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  brm_habitat_colext <- predict.deterministic.brm.habitat(
    estimates_habitat_brm_colext, T, T)
  range(brm_habitat_colext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_habitat_BRM$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_habitat_BRM)
  
  gc()
  
  ## IFM models
  div_factors <- c(1,10,100)
  load(file.path(output.path, "fits_IFM_col.Rdata") )
  load(file=file.path(output.path, "div_factor_IFM_col.Rdata") )
  estimates_ifm_col <- lapply(fits_IFM_col, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  ifm_col <- predict.deterministic.ifm(estimates_ifm_col, T, F,
                                       best_div_factors = div_factor_IFM_col,
                                       div_factors = div_factors)
  range(ifm_col$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_IFM_col$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_IFM_col)
  
  div_factors <- c(1,10,100)
  load(file.path(output.path, "fits_IFM_ext.Rdata") )
  load(file=file.path(output.path, "div_factor_IFM_ext.Rdata") )
  estimates_ifm_ext <- lapply(fits_IFM_ext, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  ifm_ext <- predict.deterministic.ifm(estimates_ifm_ext, F, T,
                                       best_div_factors = div_factor_IFM_ext,
                                       div_factors = div_factors)
  range(ifm_ext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_IFM_ext$`Aegithalos caudatus`[[3]][[4]][[1]]) ##It is not the same because on original fitting function exp_dist <0.001 is filtered out
  rm(fits_IFM_ext)
  
  div_factors <- c(1, 2.5, 5, 10, 25, 50, 100, 200)
  load(file.path(output.path, "fits_IFM_colext.Rdata") )
  load(file=file.path(output.path, "div_factor_IFM_colext.Rdata") )
  estimates_ifm_colext <- lapply(fits_IFM_colext, function(x)lapply(x, function(y)
    lapply(y, function(z) z[[2]])) )
  ifm_colext <- predict.deterministic.ifm(estimates_ifm_colext, T, T,
                                       best_div_factors = div_factor_IFM_colext,
                                       div_factors = div_factors)
  range(ifm_colext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_IFM_colext$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_IFM_colext)
  
  
  div_factors <- c(1, 2.5, 5, 10, 25, 50, 100, 200)
  load(file.path(output.path, "fits_habitat_IFM.Rdata") )
  load(file=file.path(output.path, "div_factor_IFM_colext.Rdata") )
  estimates_ifm_habitat_colext <- lapply(fits_habitat_IFM, function(x)lapply(x,
    function(y) lapply(y, function(z) z[[2]])) )
  ifm_habitat_colext <- predict.deterministic.ifm.habitat(estimates_ifm_habitat_colext, T, T,
                                          best_div_factors = div_factor_IFM_colext,
                                          div_factors = div_factors)
  range(ifm_habitat_colext$`Aegithalos caudatus`[[3]][[4]][-test_ids[[4]],] - fits_habitat_IFM$`Aegithalos caudatus`[[3]][[4]][[1]])
  rm(fits_habitat_IFM)
  
  
  models_det_occu_predictions <- list(
    "occu_base_colext" = base_colext,
    "occu_BRM_col" = brm_col,
    "occu_BRM_ext" = brm_ext,
    "occu_BRM_colext" = brm_colext,
    "occu_BRM_habitat_colext" = brm_habitat_colext,
    "occu_habitat_colext" = habitat_colext,
    "occu_IFM_col" = ifm_col,
    "occu_IFM_ext" = ifm_ext,
    "occu_IFM_colext" = ifm_colext,
    "occu_IFM_habitat_colext" = ifm_habitat_colext
    
  )
  
  save(models_det_occu_predictions, file=file.path(output.path,
                                                   "models_det_occu_predictions.Rdata"))
  
}

#### Evaluation functions ####
eval.metrics <- function(y,z){
  auc <- auc(y[,2], z[,2])
  rmse <- rmse(y[,2], z[,2])
  
  index_change <- y[,1] != y[,2]
  
  auc_change <- auc(y[index_change, 2], z[index_change, 2])
  rmse_change <- rmse(y[index_change, 2], z[index_change, 2])
  
  obs_change <- y[index_change, 2] - y[index_change, 1]
  pred_change <- z[index_change, 2] - z[index_change, 1]
  
  corr <- cor(obs_change, pred_change)
  
  eval.metrics <- c("auc"=auc, "auc_change"=auc_change,
                    "rmse"=rmse, "rmse_change"=rmse_change,
                    "corr"=corr)
  return(eval.metrics)
}


eval.trains <- function(obs, trains.occu){
  eval.list <- list()
  for ( i in 1:length(train_ids)){
    y <- obs[test_ids[[i]],]
    z <- trains.occu[[i]][test_ids[[i]],]
    eval.list[[i]] <- eval.metrics(y,z)
  }
  return(do.call(cbind, eval.list) )
}

eval.occu.predictions <- function(predicted.occu){
  eval.list <- list()
  for (study_species in study_species_vector){
    y <- species_detections[[study_species]]
    if (typeof(predicted.occu[[study_species]][[1]])=="list" ){
      eval.list[[study_species]] <- lapply(predicted.occu[[study_species]], 
                                           function (z) eval.trains(y,z))
    } else{
      eval.list[[study_species]] <- 
        eval.trains(y, predicted.occu[[study_species]])
    }
    
  }
  return(eval.list)
}

compare_metrics <- function(metrics1, metrics2){
  eval.list <- list()
  for (study_species in study_species_vector){
    y1 <- metrics1[[study_species]]
    y2 <- metrics2[[study_species]]
    if (typeof(y1) =="list" & typeof(y2) =="list"){
      eval.list[[study_species]] <- mapply(function(x,z) x - z, y1, y2)
    } else if (typeof(y1) =="list" & typeof(y2) =="double") {
      eval.list[[study_species]] <- lapply(y1, function(x) x-y2)
      
    } else if (typeof(y1) =="double" & typeof(y2) =="double") {
      eval.list[[study_species]] <- y1-y2
    } 
  }
  return(eval.list)
}

select.best.model <- function(model.metrics, metric){
  eval.list <- vector("numeric", length(study_species_vector))
  names(eval.list) <- study_species_vector
  for (study_species in study_species_vector){
    y <- model.metrics[[study_species]]
    if (typeof(y)=="list" ){
      metrics.distances <- do.call(rbind, y)
      eval.list[study_species] <- max(metrics.distances[,metric])
    } else{
      eval.list[study_species] <- max(y[metric])
    }
  }
  return(eval.list)
}


select.metric <- function(model.metrics, metric){
  eval.list <- list()
  for (study_species in study_species_vector){
    y <- model.metrics[[study_species]]
    if (typeof(y)=="list" ){
      metrics.distances <- do.call(rbind, y)
      eval.list[[study_species]] <- metrics.distances[,metric]
    } else{
      eval.list[[study_species]] <- y[metric]
    }
  }
  return(eval.list)
}

####  Evaluate models ####
base.metrics <- eval.occu.predictions(models_det_occu_predictions$occu_base_colext)

metrics.names <- c("auc", "auc_change", "rmse", "rmse_change", "corr")
# load(file= file.path(output.path, "occu_BRM_colext.Rdata") )
# BRM.colext.metrics <- eval.occu.predictions(occu_BRM_colext)
# BRM.colext.metrics.vs.base <- compare_metrics(BRM.colext.metrics, base.metrics)
metrics <- list()
metrics.vs.base <- list()
for (model.occu in names(models_det_occu_predictions)){
  occu <- models_det_occu_predictions[[model.occu]]
  metrics[[model.occu]] <- eval.occu.predictions(occu)
  metrics.vs.base[[model.occu]] <- compare_metrics(
    metrics[[model.occu]], base.metrics)
}

metrics.mean <- rapply(metrics, function(x) apply(x, 1, mean), how="list")
metrics.sd <- rapply(metrics, function(x) apply(x, 1, sd), how="list")

metrics.vs.base.mean <- rapply(metrics.vs.base, function(x) apply(x, 1, mean), how="list")
metrics.vs.base.sd <- rapply(metrics.vs.base, function(x) apply(x, 1, sd), how="list")


##Metrics best performing models
best.metrics <- list()
for (m in metrics.names){
  ##AUC difference vs base
  best.absolute <- list()
  best.diff.base <- list()
  for (model.occu in names(models_det_occu_predictions)){
    best.absolute[[model.occu]] <- select.best.model(
      metrics.mean[[model.occu]], m)
    best.diff.base[[model.occu]] <- select.best.model(
      metrics.vs.base.mean[[model.occu]], m)
  }
  best.diff.base <- as.data.frame(do.call(cbind, best.diff.base) )
  best.absolute <- as.data.frame(do.call(cbind, best.absolute) )
  best.metrics[[m]] <- list("abs"= best.absolute, "rel"= best.diff.base)
}


metrics2 <- list()
for (m in metrics.names){
  model.metrics.abs <- list()
  model.metrics.rel <- list()
  for (model.occu in names(models_det_occu_predictions)){
    model.metrics.rel[[model.occu]] <- select.metric(
      metrics.vs.base.mean[[model.occu]], m)
    model.metrics.abs[[model.occu]] <- select.metric(
      metrics.mean[[model.occu]], m)
  }
  metrics2[[m]] <- list("abs"= model.metrics.abs, "rel"= model.metrics.rel)
}
####  Figure 1: AUC boxplots ####
model_names <- c(
  "fixed",
  "spatial BRM col" ,
  "spatial BRM ext",
  "spatial BRM",
  "habitat-spatial BRM",
  "habitat colext" ,
  "spatial IFM col" ,
  "spatial IFM ext",
  "spatial IFM",
  "habitat-spatial IFM")

#model_colors <- viridis(10)[c(1,3,7,9,10)]
model_colors <- c("gray", turbo(20)[c(11,17)], viridis(20)[c(2,12,20)],
                  turbo(20)[c(12,19 )], viridis(20)[c(5,15)] )
index_all_models <-  c(1,4,9,6,5,10)

par(mar = c(2, 5, 2, 1))

# for (m in metrics.names){
  m="auc"
  boxplot(best.metrics[[m]]$abs[index_all_models], 
          names=model_names[index_all_models],
          las=2, col=model_colors[index_all_models], 
          ylab="AUC", xaxt='n')
  abline(h=median(best.metrics[[m]]$abs$occu_base_colext) )
  
  boxplot(best.metrics[[m]]$rel[index_all_models][-1], 
          names=model_names[index_all_models][-1],
          las=2, col=model_colors[index_all_models][-1], 
          ylab=expression(paste(Delta, "AUC")),
          xaxt='n')
  
  abline(h=0)
# }
  
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =model_names[index_all_models],
       pch=15,
       pt.cex=2, cex=1, bty='n',
       col = model_colors[index_all_models])

#### Figure S1: d/alpha vs AUC models ####
buffer_distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
inf_fun_distances_analysed <- c(1,2.5,5,10,20)


#model_colors <- viridis(10)[c(1,3,7,9,10)]
model_colors <- viridis(20)[c(2,5,12,15,20)]

metrics.list <- metrics2$auc$rel
var <- "cor"
max_metrics <- max(unlist(metrics.list))
min_metrics <- min(unlist(metrics.list))

#specify path to save PDF to
destination <- file.path(plots.path, "auc_negative.pdf")
#open PDF
pdf(file=destination, width = 8.3, height = 11.7)

layout(matrix(c(1:10), 5, 2, byrow = TRUE), widths = lcm(c(9,9)),
       heights = lcm(rep(5,5)))

#adjust plot margins
par(mar = c(4, 4, 2, 1))

for (i in 1:length(study_species_vector)){
  
  ##Spatial buffer
  plot( buffer_distances_analysed, metrics.list$occu_BRM_colext[[i]],
        xlim = c(1,max(inf_fun_distances_analysed)) ,
        # xlim = c(1,max(buffer_distances_analysed)) ,
        # ylim = c(0,max(meta_col_diff_aic, auto_col_diff_aic)),
        ylim = c(min_metrics*1.1,max_metrics*1.1),
        col=model_colors[1], 
        main = substitute(italic(x), list(x=study_species_vector[i])),
        type = "l", lty=2, lwd=1,
        xlab=ifelse(i %in% c(9,10,19,20,29,30,39,40,45,46),
                    expression(paste("d or ",alpha, " (km)")), "") ,
        ylab= ifelse(i%%2==0, "", expression(paste(Delta, "AUC change")) ) )
  points(buffer_distances_analysed, metrics.list$occu_BRM_colext[[i]],
         pch= 22, cex=1.5, bg=model_colors[1], col="white", lwd=0.5)
  
  ##Spatial inference function
  lines(inf_fun_distances_analysed, metrics.list$occu_IFM_colext[[i]],
        lty=2, lwd=1, col=model_colors[2])
  points(inf_fun_distances_analysed, metrics.list$occu_IFM_colext[[i]],
         pch= 23, cex=1.5, bg=model_colors[2], col="white", lwd=0.5)
  
  ##Habitat
  lines(buffer_distances_analysed, 
        rep(metrics.list$occu_habitat_colext[[i]],
            length(buffer_distances_analysed)),
        lty=1, lwd=3, col=model_colors[5])
  
  ##Habitat+Spatial brm
  lines(buffer_distances_analysed,
        metrics.list$occu_BRM_habitat_colext[[i]],
        lty=2, lwd=1, col=model_colors[3])      
  points(buffer_distances_analysed,
         metrics.list$occu_BRM_habitat_colext[[i]],
         pch= 22, cex=1.5, bg=model_colors[3], col="white", lwd=0.5)
  
  ##Habitat+Spatial inference function
  lines(inf_fun_distances_analysed,
        metrics.list$occu_IFM_habitat_colext[[i]],
        lty=2, lwd=1, col=model_colors[4])
  points(inf_fun_distances_analysed, 
         metrics.list$occu_IFM_habitat_colext[[i]],
         pch= 23, cex=1.5, bg=model_colors[4], col="white", lwd=0.5)
  
  # legend("topleft", col = c("blue", "blue", "yellow", "green", "green"),
  #        legend = c("spatial_buffer", "spatial_inf_fun", "habitat",
  #                   "habitat_spatial_buffer", "habitat_spatial_inf_fun"), 
  #        ncol = 2, 
  #        lty=c(1,2,1,1,2), lwd=3)
  abline(h=0)
  
}

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c("spatial BRM", "spatial IFM", "habitat",
                           "habitat-spatial BRM", "habitat-spatial IFM"),
       pch=c(15,18, NA, 15,18), lty = c(NA,NA,1,NA,NA),
       pt.cex=2, cex=1, bty='n', lwd= 3,
       col = model_colors[c(1,2,5,3,4)])
mtext("Legend", line=-1.5, at=0.5, cex=1)
dev.off()


#### Figure S2: Connectivity effects ####
load(file.path(output.path, "fits_BRM_colext.Rdata") )
estimates_brm_colext <- lapply(fits_BRM_colext, function(x)lapply(x, function(y)
  lapply(y, function(z) z[[2]])) )

cat_aonc2_list <- list()

for (study_species in species.info$spp){
  cat_aonc2 <- get_cat_sdm(study_species)
  cat_aonc2_list[[study_species]] <- cat_aonc2
}

mean_col_prob_sp <- list()
mean_ext_prob_sp <- list()
connectivity_sp <- list()

for (sp in 1:length(study_species_vector)){
  
  best_distance <- which.max(metrics2$auc$abs$occu_BRM_colext[[sp]] )
  
  estimates <- estimates_brm_colext[[sp]][[best_distance]]
  
  
  ##Species max and minimum connectivity
  cat_aonc_cat <- as.integer(as.character(
    cat_aonc2_list[[sp]]$aonc2_sdm_cat))
  ##Calculate neighbours info
  connectivity <- scale( sapply(index_max_dist[[best_distance]],
                                function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
  
  min_connectivity <- min(connectivity)
  max_connectivity <- max(connectivity)
  
  ##Scale connectivity 0-1
  diff_range_connectivity <- max_connectivity - min_connectivity
  ##Make connectivity values between 0 and 1
  scaled_connectivity <- (connectivity-min_connectivity)/diff_range_connectivity
  ##Get index of ordered connectivity values
  ordered_index <- order(connectivity)
  
  
  
  col_prob_trains <- matrix(nrow=length(scaled_connectivity),
                            ncol= length(estimates))
  ext_prob_trains <- matrix(nrow=length(scaled_connectivity),
                            ncol= length(estimates))
  
  for (train in 1:length(estimates)){
    col_estimates <- estimates[[train]]@estimates$col@estimates
    col_prob <- plogis(col_estimates[1] +
                         col_estimates[2]*connectivity)
    col_prob_trains[,train] <- col_prob[ordered_index]
    
    ext_estimates <- estimates[[train]]@estimates$ext@estimates
    ext_prob <- plogis(ext_estimates[1] +
                         ext_estimates[2]*connectivity)
    ext_prob_trains[,train] <- ext_prob[ordered_index]
    
  }
  
  
  mean_col_prob_sp[[sp]] <- apply(col_prob_trains, 1, mean)
  mean_ext_prob_sp[[sp]] <- apply(ext_prob_trains, 1, mean)
  connectivity_sp[[sp]] <- scaled_connectivity[ordered_index]
  
}

conn_model_colors <- turbo(20)[c(11,17)]
#adjust plot margins
par(mar = c(4, 4, 2, 1))

destination = file.path(plots.path, "col_and_ext_vs_connectivity.pdf")
pdf(file=destination, width = 8.3, height = 11.7)

layout(matrix(c(1:10), 5, 2, byrow = TRUE), widths = lcm(c(9,9)),
       heights = lcm(rep(5,5)))

for (sp in 1:length(species.info$spp)){
  sp_name <- species.info$spp[sp]
  plot(connectivity_sp[[sp]], mean_col_prob_sp[[sp]], 
       type="l", ylim=c(0,1), 
       xlab="",
       ylab= "",
       col=conn_model_colors[1], lwd=3,
       main = substitute(italic(x), list(x=sp_name)) )
  lines(connectivity_sp[[sp]], mean_ext_prob_sp[[sp]],
        col=conn_model_colors[2], lwd=3)
  legend("topleft", col = conn_model_colors,
         legend = c("col", "ext") , ncol = 2, lty=1, lwd=3)
  title(xlab="neighborhood connectivity", ylab= "predicted probability",
        line=2.2, cex.lab=1.2)
}

dev.off() 
#### Final results manuscript ####

##Mean metrics of each method range:
##AUC
range(sapply(best.metrics$auc$abs, mean) )
##AUC change
range(sapply(best.metrics$auc_change$abs, mean) )
##Correlation
range(sapply(best.metrics$corr$abs, mean) )

##Mean metrics compared to the base colext models range:
##Spatial colext
##BRM
range(best.metrics$auc$rel$occu_BRM_colext)
mean(best.metrics$auc$rel$occu_BRM_colext)
length(which(best.metrics$auc$rel$occu_BRM_colext>0))

range(best.metrics$auc_change$rel$occu_BRM_colext)
mean(best.metrics$auc_change$rel$occu_BRM_colext)
length(which(best.metrics$auc_change$rel$occu_BRM_colext>0))

range(best.metrics$corr$rel$occu_BRM_colext)
mean(best.metrics$corr$rel$occu_BRM_colext)
length(which(best.metrics$corr$rel$occu_BRM_colext>0))
##IFM
range(best.metrics$auc$rel$occu_IFM_colext)
mean(best.metrics$auc$rel$occu_IFM_colext)
length(which(best.metrics$auc$rel$occu_IFM_colext>0))

range(best.metrics$auc_change$rel$occu_IFM_colext)
mean(best.metrics$auc_change$rel$occu_IFM_colext)
length(which(best.metrics$auc_change$rel$occu_IFM_colext>0))

range(best.metrics$corr$rel$occu_IFM_colext)
mean(best.metrics$corr$rel$occu_IFM_colext)
length(which(best.metrics$corr$rel$occu_IFM_colext>0))

##BRM vs IFM
mean(abs(best.metrics$auc$abs$occu_BRM_colext - best.metrics$auc$abs$occu_IFM_colext))
range(abs(best.metrics$auc$abs$occu_BRM_colext - best.metrics$auc$abs$occu_IFM_colext))
sd(abs(best.metrics$auc$abs$occu_BRM_colext - best.metrics$auc$abs$occu_IFM_colext))
length(which((best.metrics$auc$abs$occu_BRM_colext - best.metrics$auc$abs$occu_IFM_colext)>0))

##Habitat
range(best.metrics$auc$rel$occu_habitat_colext)
mean(best.metrics$auc$rel$occu_habitat_colext)
length(which(best.metrics$auc$rel$occu_habitat_colext>0))

range(best.metrics$auc_change$rel$occu_habitat_colext)
mean(best.metrics$auc_change$rel$occu_habitat_colext)
length(which(best.metrics$auc_change$rel$occu_habitat_colext>0))

range(best.metrics$corr$rel$occu_habitat_colext)
mean(best.metrics$corr$rel$occu_habitat_colext)
length(which(best.metrics$corr$rel$occu_habitat_colext>0))


##Habitat vs spatial colext
##auc
best.spatial.auc <- pmax(best.metrics$auc$rel$occu_IFM_colext,
                         best.metrics$auc$rel$occu_BRM_colext)
best.spatial.auc[best.spatial.auc<0]<-0

best.habitat.auc <- best.metrics$auc$rel$occu_habitat_colext
best.habitat.auc[best.habitat.auc<0]<-0

length(which(best.habitat.auc>best.spatial.auc))

hab_div_spatial <- best.habitat.auc/best.spatial.auc

#### Table S2: habitat vs spatial species classification  ####

colext_effects <- data.frame(spp=species.info$spp, type=NA)
##Mainly habitat
colext_effects$type[which(hab_div_spatial>5)] <- "mainly habitat"
##Stronger habitat than connectivity
colext_effects$type[which(hab_div_spatial<5 & hab_div_spatial>2)] <- " both stronger habitat"
##Similar habitat and connectivity
colext_effects$type[which(hab_div_spatial<2 & hab_div_spatial>0.5)] <- "both similar"
##Stronger connectivity than habitat
colext_effects$type[which(hab_div_spatial<0.5 & hab_div_spatial>0.2)] <- "both stronger connectivity"
##Mainly connectivity
colext_effects$type[which(hab_div_spatial<0.2)] <- "mainly connectivity"

table(colext_effects$type)/46



# hist(best.habitat.auc/best.spatial.auc)
# index_hab_best <- which(best.habitat.auc>best.spatial.auc)
# index_spa_best <- which(best.habitat.auc<best.spatial.auc)
# 
# range((best.habitat.auc[index_hab_best]/best.spatial.auc[index_hab_best]) )
# hab_div_spatial <- best.habitat.auc[index_hab_best]/best.spatial.auc[index_hab_best]
# length(which(hab_div_spatial<3))/37
# (length(which(hab_div_spatial<9)) - length(which(hab_div_spatial<3)) )/37
# length(which(hab_div_spatial>9))/37
# 
# 
# spa_div_habitat <- best.spatial.auc[index_spa_best]/best.habitat.auc[index_spa_best]


##auc_change
best.spatial.auc_change <- pmax(best.metrics$auc_change$rel$occu_IFM_colext,
                         best.metrics$auc_change$rel$occu_BRM_colext)
length(which(best.metrics$auc_change$rel$occu_habitat_colext>best.spatial.auc_change))

##corr
best.spatial.corr <- pmax(best.metrics$corr$rel$occu_IFM_colext,
                         best.metrics$corr$rel$occu_BRM_colext)
length(which(best.metrics$corr$rel$occu_habitat_colext>best.spatial.corr))




##Habitat-spatial model
##auc
best.spatial_habitat.auc <- pmax(best.metrics$auc$rel$occu_BRM_habitat_colext,
                                 best.metrics$auc$rel$occu_IFM_habitat_colext)
best.spatial.or.habitat.auc <- pmax(best.spatial.auc,
                                    best.metrics$auc$rel$occu_habitat_colext)

##Absolute AUC increase
length(which(best.spatial_habitat.auc> best.spatial.or.habitat.auc))
range(best.spatial_habitat.auc - best.spatial.or.habitat.auc)
mean(best.spatial_habitat.auc - best.spatial.or.habitat.auc)
hist(best.spatial_habitat.auc - best.spatial.or.habitat.auc)

##Relative AUC increase
# range((best.spatial_habitat.auc - best.spatial.or.habitat.auc)/
#         best.spatial.or.habitat.auc)
# mean((best.spatial_habitat.auc - best.spatial.or.habitat.auc)/
#         best.spatial.or.habitat.auc)
# hist((best.spatial_habitat.auc - best.spatial.or.habitat.auc)/
#        best.spatial.or.habitat.auc, breaks=20)
# which(((best.spatial_habitat.auc - best.spatial.or.habitat.auc)/
#   best.spatial.or.habitat.auc)>0.5)

hab.spa.rel.increase <- best.spatial_habitat.auc/best.spatial.or.habitat.auc
range(hab.spa.rel.increase)
mean(hab.spa.rel.increase)


best.spatial.auc <- pmax(best.metrics$auc$rel$occu_IFM_colext,
                         best.metrics$auc$rel$occu_BRM_colext)
best.habitat.auc <- best.metrics$auc$rel$occu_habitat_colext

colext_effects$AUC_habitat <- round(best.habitat.auc, 3)
colext_effects$AUC_spatial <- round(best.spatial.auc, 3)

colext_effects$AUC_habitat_spatial <- round(best.spatial_habitat.auc, 3)
colext_effects$hab.spa.rel.increase <- round(hab.spa.rel.increase, 2)

colext_effects$hab.spa.abs.increase <- round(best.spatial_habitat.auc - best.spatial.or.habitat.auc, 3)

# names(colext_effects) <- c("species", "habitat vs spatial effects",
#                            "AUC difference", "prop. AUC increase")

colext_effects$type <- factor(colext_effects$type, 
                              levels = c("mainly habitat",
                                         " both stronger habitat",
                                         "both similar",
                                         "both stronger connectivity",
                                         "mainly connectivity" ) )

colext_effects <- colext_effects[order(colext_effects$type),]


write_xlsx(colext_effects, "habitat_vs_spatial.xlsx")

group_by(colext_effects, type) %>% summarise(mean=mean(hab.spa.rel.increase))

boxplot(data=colext_effects,hab.spa.rel.increase~type, las=2,
        ylab=expression(paste(Delta, "AUC")), xlab=NULL)

boxplot(data=colext_effects,hab.spa.abs.increase~type, las=2,
        ylab=expression(paste(Delta, "AUC")), xlab=NULL)


##auc_change
best.spatial_habitat.auc_change <- pmax(best.metrics$auc_change$rel$occu_BRM_habitat_colext,
                                 best.metrics$auc_change$rel$occu_IFM_habitat_colext)
best.spatial.or.habitat.auc_change <- pmax(best.spatial.auc_change,
                                    best.metrics$auc_change$rel$occu_habitat_colext)

length(which(best.spatial_habitat.auc_change> best.spatial.or.habitat.auc_change))
range(best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)
mean(best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)
hist(best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)

range((best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)/
        best.spatial.or.habitat.auc_change)
mean((best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)/
       best.spatial.or.habitat.auc_change)
hist((best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)/
       best.spatial.or.habitat.auc_change, breaks=20)
which(((best.spatial_habitat.auc_change - best.spatial.or.habitat.auc_change)/
         best.spatial.or.habitat.auc_change)>0.5)





best.spatial_habitat.corr <- pmax(best.metrics$corr$rel$occu_BRM_habitat_colext,
                                 best.metrics$corr$rel$occu_IFM_habitat_colext)
best.spatial.or.habitat.corr <- pmax(best.spatial.corr,
                                    best.metrics$corr$rel$occu_habitat_colext)


length(which(best.spatial_habitat.corr> best.spatial.or.habitat.corr))
range(best.spatial_habitat.corr - best.spatial.or.habitat.corr)
mean(best.spatial_habitat.corr - best.spatial.or.habitat.corr)
hist(best.spatial_habitat.corr - best.spatial.or.habitat.corr)

range((best.spatial_habitat.corr - best.spatial.or.habitat.corr)/
        best.spatial.or.habitat.corr)
mean((best.spatial_habitat.corr - best.spatial.or.habitat.corr)/
       best.spatial.or.habitat.corr)
hist((best.spatial_habitat.corr - best.spatial.or.habitat.corr)/
       best.spatial.or.habitat.corr, breaks=20)
which(((best.spatial_habitat.corr - best.spatial.or.habitat.corr)/
         best.spatial.or.habitat.corr)>0.5)

#### Table S1: Best selected variables  ####
##Habitat IFM vs Habitat BRM
study_species <- species.info$spp[1]
simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
                  SDM = T)

predictors <- colnames(simUMF@siteCovs)[
  -c(1, grep("sdm", colnames(simUMF@siteCovs)) )] ##Exclude cell_id and sdm



# a <- sub("^[^_]*", "", predictors)
# gsub('[_]', '', a)

pred_names <- c("forest cover",                             
                "basal area"     ,                                 
                "tree cover"      ,                                
                "tree density"     ,                                
                "dbh"               ,                     
                "tree height"        ,                        
                "closed forest" ,
                "canopy heigh avg","distance big cities",
                "coastal rocks", "coastal sands",
                "wetlands", "alpine ponds", "riparian habitats",                  
                "grasslands alpine",                 "grasslands mediterranean"           ,
                "grasslands montane",                  "shrublands alpine"                  ,
                "shrublands mediterranean",            "shrublands montane"                 ,
                "forests coniferous alpine" ,           "forests coniferous mediterranean"    ,
                "forests coniferous montane" ,          "forests deciduous alpine"            ,
                "forests deciduous montane"   ,         "forests sclerophyllous mediterranean",
                "forests sclerophyllous montane",       "riparian forests"                   ,
                "rocky mountains",                     "bare areas"                         ,
                "crops flooded"   ,                    "crops irrigated"                    ,
                "crops rainfed"    ,                   "forest plantations"                 ,
                "suburban areas"    ,                  "tree crops irrigated"                ,
                "tree crops rainfed"  ,                 "urban areas"                        ,
                "shannon index"       ,                "burnt areas"                        ,
                "dhi min"              ,               "dhi season"                            ,
                "dhi sum")         

names(pred_names) <- predictors


# load("./colext/Analysis/Rdata/best_habitat_variables_col_models.Rdata")
# load("./colext/Analysis/Rdata/best_habitat_variables_ext_models.Rdata")
# 

habitat.vars <- data.frame(matrix(ncol=3, nrow=0))
colnames(habitat.vars) <- c("species", "Colonisation", "Extinction")
same.vars.num <- c()

for (i in 1:length(best_habitat_variables_col_models)){
  sp <- names(best_habitat_variables_col_models)[i]
  col <- best_habitat_variables_col_models[[i]]$variable
  ext <- best_habitat_variables_ext_models[[i]]$variable
  
  same.vars.num <- c(same.vars.num, length(which(col %in% ext)) )
  sp.habitat.vars <- data.frame(species=sp, Colonisation=col, Extinction=ext)
  
  habitat.vars <- rbind(habitat.vars, sp.habitat.vars)
}

habitat.vars$Colonisation <- pred_names[habitat.vars$Colonisation]
habitat.vars$Extinction <- pred_names[habitat.vars$Extinction]

write_xlsx(habitat.vars, "./hab_vars.xlsx")
