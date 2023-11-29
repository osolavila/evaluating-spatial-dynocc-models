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


#### Fit base ####
load(file.path(output.path, "train_test_splits.Rdata"))

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


writeLines(c(""), "log.txt")

fits_base <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {

  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  simUMF <- species_UMFs[[study_species]]
  
  train.fits <- list()
  for (i in 1:length(test_ids) ){
    
    test <- test_ids[[i]]
    simUMF0<- simUMF
    simUMF0@y[test,] <- NA
    ##fit models
    spatial.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                                gammaformula = ~ 1, # Colonization
                                epsilonformula = ~ 1, # Extinction
                                pformula = det.formula, # Detection
                                data = simUMF0) )
    
    model.list <- list()
    if (class(spatial.model) == "unmarkedFitColExt"){
      model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
      model.list[[2]] <- spatial.model@estimates
      model.list[[3]] <- spatial.model@AIC
      # model.list[[4]] <- spatial.model@formula
    }
    
    train.fits[[i]] <- model.list
  
  }

  return(train.fits)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_base) <- species.info$spp

save(fits_base,
     file=file.path(output.path, "fits_base.Rdata"))


