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
##Habitat models covariables
study_species <- species.info$spp[1]
simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
                  SDM = T)

predictors <- colnames(simUMF@siteCovs)[
  -c(1, grep("sdm", colnames(simUMF@siteCovs)) )] ##Exclude cell_id and sdm


num_vars <- 4


species_UMFs<- list()
for (study_species in species.info$spp){
  simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
                    SDM = T)
  species_UMFs[[study_species]] <- simUMF
}


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


#### Fit colext models with best habitat variables  ####

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

writeLines(c(""), "log.txt")
fits_habitat <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  cat(paste("Analysing species", study_species,"\n"), file="log.txt", append=TRUE)
  simUMF0 <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  gamma <- as.formula(col.equations[study_species])
  epsilon <- as.formula(ext.equations[study_species])
  
  train.fits <- list()
  for (t in 1:length(test_ids) ){
    
    test <- test_ids[[t]]
    simUMF<- simUMF0
    simUMF@y[test,] <- NA
  
  
    habitat.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                                 gammaformula = gamma, # Colonization
                                 epsilonformula = epsilon, # Extinction
                                 pformula = det.formula, # Detection
                                 data = simUMF), silent = T )
    
    nP <- (habitat.model@AIC - habitat.model@negLogLike*2)/2
    
    starts_list <- list()
    starts_list[[1]] <- rep(0, nP)
    habitat.models <- list()
    model.list <- list()
    if (class(habitat.model) == "unmarkedFitColExt"){
      model.list[[1]] <- aperm(habitat.model@projected[2,,], c(2,1)) ##matrix sites - seasons
      model.list[[2]] <- habitat.model@estimates
      model.list[[3]] <- habitat.model@AIC
      habitat.models[[1]] <- model.list
    } else{
      habitat.models[[1]] <- list(NA,NA,NA)
    }
    
     
    for (i in 2:10){
      print(i)
      # if (is.null(starts)) 
      if (i==1){
        starts <- rep(0, nP)
      } else{
        starts <- runif(nP, min=-1.5, max=1.5)
      }
      
      ##Model fitting
      starts_list[[i]] <- starts
      
      habitat.model <- try(
        colext(psiformula = psi_formula, # First-year occupancy
               gammaformula = gamma, # Colonization
               epsilonformula = epsilon, # Extinction
               pformula = det.formula, # Detection
               data = simUMF,
               starts = starts) )
      
      if (class(habitat.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(habitat.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- habitat.model@estimates
        model.list[[3]] <- habitat.model@AIC
        habitat.models[[i]] <- model.list
      } else{
        habitat.models[[i]] <- list(NA,NA,NA)
      }
    }
    
    index_best <- which.min(sapply(habitat.models, function(x) x[[3]] ) )
    
    habitat.model <- habitat.models[[index_best]]
    starts <- starts_list[[index_best]]
    
    train.fits[[t]] <- list(habitat.model, study_species, starts)
  
  }
  
  return(train.fits)
  
  
}

parallel::stopCluster(cl = my.cluster)

sp.names <- species.info$spp

starts_habitat <- lapply(fits_habitat,
                             function(x) 
                               lapply(x, function(y) y[[3]]))

names(starts_habitat) <- sp.names
save(starts_habitat, 
     file=file.path(output.path, "starts_habitat.Rdata"))


fits_habitat <- lapply(fits_habitat, function(x) 
  lapply(x, function(y) y[[1]]))

names(fits_habitat) <- sp.names

save(fits_habitat, 
     file=file.path(output.path, "fits_habitat.Rdata"))

