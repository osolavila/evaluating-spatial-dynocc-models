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

library(foreach)
library(doParallel)
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

#### Model parameters ####
intercept <- "~ 1"
atlas <- "~ atlas"
quadratic_date <- "~ julian_day_std + I(julian_day_std^2)"
quadratic_date_atlas <- "~ atlas + julian_day_std + I(julian_day_std^2)"
interaction_atlas_quadratic_date <- "~ atlas + julian_day_std + 
I(julian_day_std^2) + atlas:julian_day_std + atlas:I(julian_day_std^2)"

det_model_formulas <- c(intercept, atlas, quadratic_date, quadratic_date_atlas, 
                        interaction_atlas_quadratic_date)
names(det_model_formulas) <- c("intercept", "atlas", "quadratic_date",
                               "quadratic_date_atlas",
                               "interaction_atlas_quadratic_date")

psi_formula <- formula(" ~ aonc2_sdm")

load(file.path(output.path, "train_test_splits.Rdata"))
#### Fit Metapopulation models####
source("./colext/Functions/spatial_colext_v2.R")

environment(colext.fit2) <- asNamespace('unmarked')
assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')

environment(colext2) <- asNamespace('unmarked')
assignInNamespace("colext", colext2, ns = 'unmarked')

#### Fit col spatial IFM ####

parallel::detectCores()

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

writeLines(c(""), "log.txt")

fits_IFM_col <- foreach(
  study_species = species.info$spp, 
  .packages = c("unmarked", "ParallelSpatialColext"),
  .noexport = 'rcpp_y'
) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  cat(paste("Analysing species", study_species,"\n"), file="log.txt", append=TRUE)
  
  environment(colext.fit2) <- asNamespace('unmarked')
  assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')
  
  environment(colext2) <- asNamespace('unmarked')
  assignInNamespace("colext", colext2, ns = 'unmarked')
  

  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  ##SDM AONC2 
  cat_sdm_spp <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  
  nP <- get.nP(psiformula = psi_formula, # First-year occupancy
               gammaformula = ~ 1, # Colonization
               epsilonformula = ~ 1, # Extinction
               pformula = det.formula)
  
  for (alpha in c(1,2.5,5,10,20) ){
    
    print(paste0("Analysing ", study_species, ", dist: ", alpha))
    
    train.fits <- list()
    for (t in 1:length(test_ids) ){
      
      test <- test_ids[[t]]
      simUMF<- simUMF0
      simUMF@y[test,] <- NA
      
      best_model <- list()
      best_starts <- list()
      div_i <-1
      for (div in c(1,10,100)){
        ##Calculates exp kernel multiplied by psi of neighbours
        exp_disp_times_aonc <- list()
        for (i in 1:length(distance_cells_30km)){
          exp_times_aonc <-
            exp(-distance_cells_30km[[i]]/alpha)*
            cat_sdm_spp$aonc2_sdm[index_cells_30km[[i]]]/div
          exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
        }
        
        
        ##Inference function for colonization
        col_exp_disp_times_aonc <- exp_disp_times_aonc
        
        ##fit models
        spatial.model <- try(unmarked::colext(
          psiformula = psi_formula, # First-year occupancy
          gammaformula = ~ 1, # Colonization
          epsilonformula = ~ 1, # Extinction
          pformula = det.formula, # Detection
          data = simUMF,
          spatial.col = T,
          col_exp_disp_times_aonc= col_exp_disp_times_aonc) )
        
        starts_list <- list()
        starts_list[[1]] <- rep(0, nP)
        
        spatial.models <- list()
        model.list <- list()
        
        if (class(spatial.model) == "unmarkedFitColExt"){
          model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
          model.list[[2]] <- spatial.model@estimates
          model.list[[3]] <- spatial.model@AIC
          spatial.models[[1]] <- model.list
        } else{
          spatial.models[[1]] <- list(NA,NA,NA)
        }
        
        
        for (i in 2:3){
          
          # if (is.null(starts)) 
          if (i==1){
            starts <- rep(0, nP)
          } else{
            starts <- runif(nP, min=-1.5, max=1.5)
          }
          
          ##Model fitting
          starts_list[[i]] <- starts
          
          spatial.model <- try(unmarked::colext(
            psiformula = psi_formula, # First-year occupancy
            gammaformula = ~ 1, # Colonization
            epsilonformula = ~ 1, # Extinction
            pformula = det.formula, # Detection
            data = simUMF,
            starts = starts,
            spatial.col = T,
            col_exp_disp_times_aonc= col_exp_disp_times_aonc) )
          
          if (class(spatial.model) == "unmarkedFitColExt"){
            model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
            model.list[[2]] <- spatial.model@estimates
            model.list[[3]] <- spatial.model@AIC
            spatial.models[[i]] <- model.list
          } else{
            spatial.models[[i]] <- list(NA,NA,NA)
          }
        }
        
        index_best <- which.min(sapply(spatial.models, function(x) x[[3]] ) )
        
        best_model[[div_i]] <- spatial.models[[index_best]]
        best_starts[[div_i]] <- starts_list[[index_best]]
        
        div_i <- div_i + 1
      }
      
      index_best <- which.min(sapply(best_model, function(x) x[[3]] ) )
      
      spatial.model <- best_model[[index_best]]
      starts <- best_starts[[index_best]]
      
      train.fits[[t]] <- list(spatial.model, index_best, starts)
    
    }
    
    dist.models[[as.character(alpha)]] <- train.fits
    
  }
  
  
  
  return(list(dist.models, study_species))
  
  
}

parallel::stopCluster(cl = my.cluster)


sp.names <- sapply(fits_IFM_col, function(x) x[[2]])

fits_IFM_col <- lapply(fits_IFM_col, function(x) x[[1]])
names(fits_IFM_col) <- sp.names

starts_IFM_col <- lapply(fits_IFM_col, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[3]])) )

div_factor_IFM_col <- lapply(fits_IFM_col, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )

fits_IFM_col <- lapply(fits_IFM_col, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

save(starts_IFM_col, 
     file=file.path(output.path, "starts_IFM_col.Rdata"))

save(div_factor_IFM_col,
     file=file.path(output.path, "div_factor_IFM_col.Rdata"))

save(fits_IFM_col, 
     file=file.path(output.path, "fits_IFM_col.Rdata"))


#### Fit ext spatial IFM ####

parallel::detectCores()

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

writeLines(c(""), "log.txt")

fits_IFM_ext <- foreach(
  study_species = species.info$spp, 
  .packages = c("unmarked", "ParallelSpatialColext"),
  .noexport = 'rcpp_y'
) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  cat(paste("Analysing species", study_species,"\n"), file="log.txt", append=TRUE)
  
  environment(colext.fit2) <- asNamespace('unmarked')
  assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')
  
  environment(colext2) <- asNamespace('unmarked')
  assignInNamespace("colext", colext2, ns = 'unmarked')
  
  
  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  ##SDM AONC2 
  cat_sdm_spp <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  
  nP <- get.nP(psiformula = psi_formula, # First-year occupancy
               gammaformula = ~ 1, # Colonization
               epsilonformula = ~ 1, # Extinction
               pformula = det.formula)
  
  for (alpha in c(1,2.5,5,10,20) ){
    
    print(paste0("Analysing ", study_species, ", dist: ", alpha))
    
    train.fits <- list()
    for (t in 1:length(test_ids) ){
      
      test <- test_ids[[t]]
      simUMF<- simUMF0
      simUMF@y[test,] <- NA
    
      best_model <- list()
      best_starts <- list()
      div_i <-1
      for (div in c(1,10,100)){
        ##Calculates exp kernel multiplied by psi of neighbours
        exp_disp_times_aonc <- list()
        for (i in 1:length(distance_cells_30km)){
          exp_times_aonc <-
            exp(-distance_cells_30km[[i]]/alpha)*
            cat_sdm_spp$aonc2_sdm[index_cells_30km[[i]]]/div
          exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
        }
        
        
        ##fit models
        spatial.model <- try(unmarked::colext(
          psiformula = psi_formula, # First-year occupancy
          gammaformula = ~ 1, # Colonization
          epsilonformula = ~ 1, # Extinction
          pformula = det.formula, # Detection
          data = simUMF,
          spatial.ext = T,
          ext_exp_disp_times_aonc= exp_disp_times_aonc) )
        
        starts_list <- list()
        starts_list[[1]] <- rep(0, nP)
        
        spatial.models <- list()
        model.list <- list()
        
        if (class(spatial.model) == "unmarkedFitColExt"){
          model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
          model.list[[2]] <- spatial.model@estimates
          model.list[[3]] <- spatial.model@AIC
          spatial.models[[1]] <- model.list
        } else{
          spatial.models[[1]] <- list(NA,NA,NA)
        }
        
        
        for (i in 2:3){
          
          # if (is.null(starts)) 
          if (i==1){
            starts <- rep(0, nP)
          } else{
            starts <- runif(nP, min=-1.5, max=1.5)
          }
          
          ##Model fitting
          starts_list[[i]] <- starts
          
          spatial.model <- try(unmarked::colext(
            psiformula = psi_formula, # First-year occupancy
            gammaformula = ~ 1, # Colonization
            epsilonformula = ~ 1, # Extinction
            pformula = det.formula, # Detection
            data = simUMF,
            starts = starts,
            spatial.ext = T,
            ext_exp_disp_times_aonc= exp_disp_times_aonc) )
          
          if (class(spatial.model) == "unmarkedFitColExt"){
            model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
            model.list[[2]] <- spatial.model@estimates
            model.list[[3]] <- spatial.model@AIC
            spatial.models[[i]] <- model.list
          } else{
            spatial.models[[i]] <- list(NA,NA,NA)
          }
        }
        
        index_best <- which.min(sapply(spatial.models, function(x) x[[3]] ) )
        
        best_model[[div_i]] <- spatial.models[[index_best]]
        best_starts[[div_i]] <- starts_list[[index_best]]
        
        div_i <- div_i + 1
      }
      
      index_best <- which.min(sapply(best_model, function(x) x[[3]] ) )
      
      spatial.model <- best_model[[index_best]]
      starts <- best_starts[[index_best]]
      
      train.fits[[t]] <- list(spatial.model, index_best, starts)
      
    }
    
    dist.models[[as.character(alpha)]] <- train.fits
    
  }
  
  
  
  return(list(dist.models, study_species))
  
  
}

parallel::stopCluster(cl = my.cluster)


sp.names <- sapply(fits_IFM_ext, function(x) x[[2]])

fits_IFM_ext <- lapply(fits_IFM_ext, function(x) x[[1]])
names(fits_IFM_ext) <- sp.names

starts_IFM_ext <- lapply(fits_IFM_ext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[3]])) )

div_factor_IFM_ext <- lapply(fits_IFM_ext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )

fits_IFM_ext <- lapply(fits_IFM_ext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

save(starts_IFM_ext, 
     file=file.path(output.path, "starts_IFM_ext.Rdata"))

save(div_factor_IFM_ext,
     file=file.path(output.path, "div_factor_IFM_ext.Rdata"))


save(fits_IFM_ext, 
     file=file.path(output.path, "fits_IFM_ext.Rdata"))


#### Fit colext spatial IFM ####

parallel::detectCores()

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

writeLines(c(""), "log.txt")

fits_IFM_colext <- foreach(
  study_species = species.info$spp, 
  .packages = c("unmarked", "ParallelSpatialColext"),
  .noexport = 'rcpp_y'
) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  cat(paste("Analysing species", study_species,"\n"), file="log.txt", append=TRUE)
  
  environment(colext.fit2) <- asNamespace('unmarked')
  assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')
  
  environment(colext2) <- asNamespace('unmarked')
  assignInNamespace("colext", colext2, ns = 'unmarked')
  
  
  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  ##SDM AONC2 
  cat_sdm_spp <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  
  nP <- get.nP(psiformula = psi_formula, # First-year occupancy
               gammaformula = ~ 1, # Colonization
               epsilonformula = ~ 1, # Extinction
               pformula = det.formula)
  
  for (alpha in c(1,2.5,5,10,20) ){
    
    print(paste0("Analysing ", study_species, ", dist: ", alpha))
    
    train.fits <- list()
    for (t in 1:length(test_ids) ){
      
      test <- test_ids[[t]]
      simUMF<- simUMF0
      simUMF@y[test,] <- NA
    
      best_model <- list()
      best_starts <- list()
      div_i <-1
      for (div in c(1, 2.5, 5, 10, 25, 50, 100, 200)){
        ##Calculates exp kernel multiplied by psi of neighbours
        exp_disp_times_aonc <- list()
        for (i in 1:length(distance_cells_30km)){
          exp_times_aonc <-
            exp(-distance_cells_30km[[i]]/alpha)*
            cat_sdm_spp$aonc2_sdm[index_cells_30km[[i]]]/div
          exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
        }
        
        
        ##fit models
        spatial.model <- try(unmarked::colext(
          psiformula = psi_formula, # First-year occupancy
          gammaformula = ~ 1, # Colonization
          epsilonformula = ~ 1, # Extinction
          pformula = det.formula, # Detection
          data = simUMF,
          spatial.col = T,
          spatial.ext = T,
          col_exp_disp_times_aonc= exp_disp_times_aonc,
          ext_exp_disp_times_aonc = exp_disp_times_aonc) )
        
        starts_list <- list()
        starts_list[[1]] <- rep(0, nP)
        
        spatial.models <- list()
        model.list <- list()
        
        if (class(spatial.model) == "unmarkedFitColExt"){
          model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
          model.list[[2]] <- spatial.model@estimates
          model.list[[3]] <- spatial.model@AIC
          spatial.models[[1]] <- model.list
        } else{
          spatial.models[[1]] <- list(NA,NA,NA)
        }
        
        
        for (i in 2:3){
          
          # if (is.null(starts)) 
          if (i==1){
            starts <- rep(0, nP)
          } else{
            starts <- runif(nP, min=-1.5, max=1.5)
          }
          
          ##Model fitting
          starts_list[[i]] <- starts
          
          ##fit models
          spatial.model <- try(unmarked::colext(
            psiformula = psi_formula, # First-year occupancy
            gammaformula = ~ 1, # Colonization
            epsilonformula = ~ 1, # Extinction
            pformula = det.formula, # Detection
            data = simUMF,
            starts = starts,
            spatial.col = T,
            spatial.ext = T,
            col_exp_disp_times_aonc= exp_disp_times_aonc,
            ext_exp_disp_times_aonc = exp_disp_times_aonc) )
          
          if (class(spatial.model) == "unmarkedFitColExt"){
            model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
            model.list[[2]] <- spatial.model@estimates
            model.list[[3]] <- spatial.model@AIC
            spatial.models[[i]] <- model.list
          } else{
            spatial.models[[i]] <- list(NA,NA,NA)
          }
        }
        
        index_best <- which.min(sapply(spatial.models, function(x) x[[3]] ) )
        
        best_model[[div_i]] <- spatial.models[[index_best]]
        best_starts[[div_i]] <- starts_list[[index_best]]
        
        div_i <- div_i + 1
      }
      
      index_best <- which.min(sapply(best_model, function(x) x[[3]] ) )
      
      spatial.model <- best_model[[index_best]]
      starts <- best_starts[[index_best]]
      
      train.fits[[t]] <- list(spatial.model, index_best, starts)
      
    }
    
    dist.models[[as.character(alpha)]] <- train.fits
    
  }

  return(list(dist.models, study_species))
  
}

parallel::stopCluster(cl = my.cluster)


sp.names <- sapply(fits_IFM_colext, function(x) x[[2]])

fits_IFM_colext <- lapply(fits_IFM_colext, function(x) x[[1]])
names(fits_IFM_colext) <- sp.names

starts_IFM_colext <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[3]])) )

div_factor_IFM_colext <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )

fits_IFM_colext <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

save(starts_IFM_colext, 
     file=file.path(output.path, "starts_IFM_colext.Rdata"))

save(div_factor_IFM_colext,
     file=file.path(output.path, "div_factor_IFM_colext.Rdata"))


save(fits_IFM_colext, 
     file=file.path(output.path, "fits_IFM_colext.Rdata"))