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
  cat_aonc_cat <- as.integer(as.character(cat_aonc2$aonc2_sdm_cat))
  cat_aonc2_list[[study_species]] <- cat_aonc_cat
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

distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
index_max_dist <- list()
names_dist <- paste0("max_dist_", c(1:3, seq(5,10,2.5), seq(15,30,5)) )
i <- 1
for (max_dist in c(1:3, seq(5,10,2.5), seq(15,30,5))+ 0.2 ){
  index_max_dist[[names_dist[i]]] <- apply(points_distances,1,function(x)
    which(x<max_dist))
  i<-i+1
}
rm(i)
gc()

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

autocovariate_formula <- formula(" ~ sp.autocovariate")

load(file.path(output.path, "train_test_splits.Rdata"))
#### Fit brm col ####

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


writeLines(c(""), "log.txt")

fits_BRM_col <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  
  
  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  cat_aonc_cat <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  for (dist in names(index_max_dist)) {
    
    ##Calculate neighbours info
    sp.autocovariate <- scale( sapply(index_max_dist[[dist]],
                                      function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
    
    ##Add neighbour info to colext fitting data
    siteCovs <- cbind(simUMF0@siteCovs, sp.autocovariate)
    simUMF <- unmarkedMultFrame(
      y = simUMF0@y,
      siteCovs = siteCovs,
      obsCovs = simUMF0@obsCovs,
      yearlySiteCovs=simUMF0@yearlySiteCovs,
      numPrimary=simUMF0@numPrimary)
    
    train.fits <- list()
    for (i in 1:length(test_ids) ){
      
      test <- test_ids[[i]]
      simUMF_train<- simUMF
      simUMF_train@y[test,] <- NA
      ##fit models
      spatial.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                                  gammaformula = autocovariate_formula, # Colonization
                                  epsilonformula = ~ 1, # Extinction
                                  pformula = det.formula, # Detection
                                  data = simUMF_train) )
      
      model.list <- list()
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        train.fits[[i]] <- model.list
      } else{
        train.fits[[i]] <- list(NA, NA, NA)
      }
      
    }
    dist.models[[dist]] <- train.fits
  }
  return(dist.models)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_BRM_col) <- species.info$spp

save(fits_BRM_col,
     file=file.path(output.path, "fits_BRM_col.Rdata"))


#### Fit brm ext  ####

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


writeLines(c(""), "log.txt")

fits_BRM_ext <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  cat_aonc_cat <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  for (dist in names(index_max_dist)) {
    print(paste0("Analysing ", study_species, ", dist: ", dist))
    
    ##Calculate neighbours info
    sp.autocovariate <- scale( sapply(index_max_dist[[dist]],
                                      function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
    
    ##Add neighbour info to colext fitting data
    siteCovs <- cbind(simUMF0@siteCovs, sp.autocovariate)
    simUMF <- unmarkedMultFrame(
      y = simUMF0@y,
      siteCovs = siteCovs,
      obsCovs = simUMF0@obsCovs,
      yearlySiteCovs=simUMF0@yearlySiteCovs,
      numPrimary=simUMF0@numPrimary)
    
    train.fits <- list()
    for (i in 1:length(test_ids) ){
      
      test <- test_ids[[i]]
      simUMF_train<- simUMF
      simUMF_train@y[test,] <- NA
    
      spatial.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                                  gammaformula = ~ 1 , # Colonization
                                  epsilonformula = autocovariate_formula, # Extinction
                                  pformula = det.formula, # Detection
                                  data = simUMF_train) )
      
      model.list <- list()
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        train.fits[[i]] <- model.list
      } else{
        train.fits[[i]] <- list(NA, NA, NA)
      }
    
    }
    dist.models[[dist]] <- train.fits
  }
  
  return(dist.models)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_BRM_ext) <- species.info$spp

save(fits_BRM_ext,
     file=file.path(output.path, "fits_BRM_ext.Rdata"))



####  Fit brm colext  ####
n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


writeLines(c(""), "log.txt")

fits_BRM_colext <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {

  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  cat_aonc_cat <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  for (dist in names(index_max_dist)) {
    print(paste0("Analysing ", study_species, ", dist: ", dist))
    
    ##Calculate neighbours info
    sp.autocovariate <- scale( sapply(index_max_dist[[dist]],
                                      function(x) sum(cat_aonc_cat[x], na.rm=T)/length(x)) )
    
    ##Add neighbour info to colext fitting data
    siteCovs <- cbind(simUMF0@siteCovs, sp.autocovariate)
    simUMF <- unmarkedMultFrame(
      y = simUMF0@y,
      siteCovs = siteCovs,
      obsCovs = simUMF0@obsCovs,
      yearlySiteCovs=simUMF0@yearlySiteCovs,
      numPrimary=simUMF0@numPrimary)
    
    train.fits <- list()
    for (i in 1:length(test_ids) ){
      
      test <- test_ids[[i]]
      simUMF_train<- simUMF
      simUMF_train@y[test,] <- NA
      
      spatial.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                                  gammaformula = autocovariate_formula, # Colonization
                                  epsilonformula = autocovariate_formula, # Extinction
                                  pformula = det.formula, # Detection
                                  data = simUMF_train) )
      
      model.list <- list()
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        train.fits[[i]] <- model.list
      } else{
        train.fits[[i]] <- list(NA, NA, NA)
      }
      
    }
    dist.models[[dist]] <- train.fits
  }
  
  return(dist.models)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_BRM_colext) <- species.info$spp

save(fits_BRM_colext,
     file=file.path(output.path, "fits_BRM_colext.Rdata"))
