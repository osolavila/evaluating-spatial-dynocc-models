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
##SpeciesAONC2 SDM for all Catalonia
cat_aonc2_list <- list()

for (study_species in species.info$spp){
  simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
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

##Div factor for neigh col prob 
##To find optimal solution with starts range (-1.5,1.5)
load(file.path(output.path, "div_factor_IFM_colext.Rdata"))
div_factors <- c(1, 2.5, 5, 10, 25, 50, 100, 200)

#### Fit Metapopulation models####
source("./colext/Functions/spatial_colext_v2.R")

environment(colext.fit2) <- asNamespace('unmarked')
assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')

environment(colext2) <- asNamespace('unmarked')
assignInNamespace("colext", colext2, ns = 'unmarked')

####  Fit colext habitat ifm models ####
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

fits_habitat_IFM <- foreach(
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
  
  epsilon <- as.formula(ext.equations[study_species])
  gamma <- as.formula(col.equations[study_species])
  
  ##SDM AONC2 
  cat_sdm_spp <- cat_aonc2_list[[study_species]]
  
  dist.models <- list()
  
  nP <- get.nP(psiformula = psi_formula, # First-year occupancy
               gammaformula = gamma, # Colonization
               epsilonformula = epsilon, # Extinction
               pformula = det.formula)
  
  for (alpha in c(1,2.5,5,10,20) ){
    
    print(paste0("Analysing ", study_species, ", alpha: ", alpha))
    
    train.fits <- list()
    for (t in 1:length(test_ids) ){
      
      test <- test_ids[[t]]
      simUMF<- simUMF0
      simUMF@y[test,] <- NA
      
      model.list <- list()
      
      ##Calculates exp kernel multiplied by psi of neighbours
      exp_disp_times_aonc <- list()
      div <- div_factors[
        div_factor_IFM_colext[[study_species]][[
        as.character(alpha)]][[t]]
      ]
      for (i in 1:length(distance_cells_30km)){
        exp_times_aonc <-
          exp(-distance_cells_30km[[i]]/alpha)*
          cat_sdm_spp$aonc2_sdm[index_cells_30km[[i]]]/div
        exp_disp_times_aonc[[i]] <- exp_times_aonc[exp_times_aonc>0.001]
      }
      
      ##fit models
      spatial.model <- try(unmarked::colext(
        psiformula = psi_formula, # First-year occupancy
        gammaformula = gamma, # Colonization
        epsilonformula = epsilon, # Extinction
        pformula = det.formula, # Detection
        data = simUMF,
        spatial.col = T,
        spatial.ext = T,
        col_exp_disp_times_aonc= exp_disp_times_aonc,
        ext_exp_disp_times_aonc = exp_disp_times_aonc) )
      
      starts_list <- list()
      starts_list[[1]] <- rep(0, nP)
      spatial.models <- list()
      
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        spatial.models[[1]] <- model.list
        
      } else{
        spatial.models[[1]] <- list(NA,NA,NA)
      }
      
      
      for (i in 2:10){
        
        # if (is.null(starts)) 
        if (i==1){
          starts <- rep(0, nP)
        } else{
          starts <- runif(nP, min=-2, max=2)
        }
        
        ##Model fitting
        starts_list[[i]] <- starts
        
        ##fit models
        spatial.model <- try(unmarked::colext(
          psiformula = psi_formula, # First-year occupancy
          gammaformula = gamma, # Colonization
          epsilonformula = epsilon, # Extinction
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
      
      spatial.model <- spatial.models[[index_best]]
      starts <- starts_list[[index_best]]
      
      train.fits[[t]] <- list(spatial.model, index_best, starts)
      
    }
    
    dist.models[[as.character(alpha)]] <- train.fits
    
  }
  

  return(list(dist.models, study_species))
  
  
}

parallel::stopCluster(cl = my.cluster)


sp.names <- sapply(fits_habitat_IFM, function(x) x[[2]])

fits_habitat_IFM <- lapply(fits_habitat_IFM,
                           function(x) x[[1]])
names(fits_habitat_IFM) <- sp.names

starts_habitat_IFM <- lapply(fits_habitat_IFM, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[3]])) )

fits_habitat_IFM <- lapply(fits_habitat_IFM, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

save(starts_habitat_IFM, 
     file=file.path(output.path, "starts_habitat_IFM.Rdata"))

save(fits_habitat_IFM, 
     file=file.path(output.path, "fits_habitat_IFM.Rdata"))

