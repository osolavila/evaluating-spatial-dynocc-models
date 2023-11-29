#### Load detection data and analysis functions ####
Dir.Base <- "~/evaluating-spatial-dynamic-occupancy-models"
setwd(Dir.Base)

##Load detection data
source(file.path(Dir.Base, "Paper_scripts", "data_preparation",
                 "Prepare_species_detecion_data.R"))
##Load useful functions
source(file.path(Dir.Base, "Paper_scripts", "Functions",
                 "get_colext_umf.R"))

##Load libraries to parallelise code
library(foreach)
library(doParallel)


#### Gather species info and select study species ####

## DF: Species|#AONC2|#Col|#Ext|%change|Habitat_Group
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

##Select species with detection >0.5 in both sampling periods
species.index <- which(species.info$Atlas_2>0.5 &
                         species.info$Atlas_3>0.5)
species.info <- species.info[species.index,]


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

#### Fit colonization habitat models ####


parallel::detectCores()

n.cores <- parallel::detectCores() - 2

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()
  

clusterExport(cl=my.cluster, c('species_UMFs', 'predictors', 'num_vars',
                               'species.info', 'det_model_formulas',
                               'psi_formula'))

writeLines(c(""), "log.txt")

best_habitat_variables_col_models <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  
  
  sink("log.txt", append=TRUE)
  cat(paste("Analysing species", study_species,"\n"))
  
  simUMF <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  ##Select best habitat variables
  cols.names <- c("iteration", "variable", "AIC", "Quadratic_term")
  vars.selected.AIC <- data.frame(matrix(nrow=num_vars, ncol=length(cols.names)))
  colnames(vars.selected.AIC) <- cols.names
  
  
  
  for (i in 1:num_vars){
    
    print(paste0("Analysing ", study_species, ", loop: ", i))
    
    if (i == 1) {
      selected_formula <- " ~ "
      vars_selected_index <- c()
      }
    
    preds_to_analyse <- setdiff(1:length(predictors), vars_selected_index)
    l.models.aic <- c()
    q.models.aic <- c()
    for (pred in predictors[preds_to_analyse]) {
      if (i== 1){
        l_form <- formula(paste0(selected_formula, pred ))
        q_form <- formula(paste0(
          selected_formula, pred , " + I(", pred, "^2)"))
      } else{
        l_form <- formula(paste0(selected_formula, " + ", pred ))
        q_form <- formula(paste0(
          selected_formula, " + ", pred , " + I(", pred, "^2)"))
      }
      
      l.models.aic <- c(l.models.aic,
                        try(colext(psiformula = psi_formula, # First-year occupancy
                                   gammaformula = l_form, # Colonization
                                   epsilonformula = ~ 1, # Extinction
                                   pformula = det.formula, # Detection
                                   data = simUMF)@AIC ))
      
      q.models.aic <- c(q.models.aic,
                        try(colext(psiformula = psi_formula, # First-year occupancy
                                   gammaformula = q_form, # Colonization
                                   epsilonformula = ~ 1, # Extinction
                                   pformula = det.formula, # Detection
                                   data = simUMF)@AIC ))
      
    }
    
    best_aic <- min(mapply(
      function(x,y) min(x,y), l.models.aic, q.models.aic))
    index_best <- which.min(mapply(
      function(x,y) min(x,y), l.models.aic, q.models.aic))
    vars_selected_index <- c(vars_selected_index, preds_to_analyse[index_best]) ##Keep track of variables selected
    quadratic_term <- l.models.aic[index_best] > q.models.aic[index_best] ##T/F
    
    if (i>1){
      if (best_aic > vars.selected.AIC$AIC[i-1]){
        break
      }
    }
    
    pred <- predictors[preds_to_analyse[index_best]]
    vars.selected.AIC$iteration[i] <- i
    vars.selected.AIC$variable[i] <- pred
    vars.selected.AIC$AIC[i] <- best_aic
    vars.selected.AIC$Quadratic_term[i] <- quadratic_term
    
    
    ##Update formula
    if (i == 1){
      if (quadratic_term) {
        selected_formula <- paste0(selected_formula, pred , " + I(", pred, "^2)")
      } else{
        selected_formula <- paste0(selected_formula, pred )
      }
    }  else{ ## loop 2 to num_vars
      if (quadratic_term) {
        selected_formula <- paste0(
          selected_formula, " + ", pred , " + I(", pred, "^2)") 
      } else{
        selected_formula <- paste0(selected_formula, " + ", pred )
      }
    }
    
    
  }
  
  return(list(vars.selected.AIC, study_species) )
  
  
}

parallel::stopCluster(cl = my.cluster)



sp.names <- sapply(best_habitat_variables_col_models, function(x) x[[2]])
best_habitat_variables_col_models <- lapply(best_habitat_variables_col_models,
                                            function(x) x[[1]])
names(best_habitat_variables_col_models) <- sp.names

save(best_habitat_variables_col_models, 
     file="./colext/Analysis/Rdata/best_habitat_variables_col_models.Rdata")

#### Fit extinction habitat models ####


parallel::detectCores()

n.cores <- parallel::detectCores() - 6

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()



species_UMFs<- list()
for (study_species in species.info$spp){
  simUMF <- get_UMF(study_species, LiDAR = T, habitat_aonc3 = T,
                    SDM = T)
  species_UMFs[[study_species]] <- simUMF
}


clusterExport(cl=my.cluster, c('species_UMFs', 'predictors', 'num_vars',
                               'species.info', 'det_model_formulas',
                               'psi_formula'))

writeLines(c(""), "log.txt")

best_habitat_variables_ext_models <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  
  
  sink("log.txt", append=TRUE)
  cat(paste("Analysing species", study_species,"\n"))
  
  simUMF <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  ##Select best habitat variables
  cols.names <- c("iteration", "variable", "AIC", "Quadratic_term")
  vars.selected.AIC <- data.frame(matrix(nrow=num_vars, ncol=length(cols.names)))
  colnames(vars.selected.AIC) <- cols.names
  
  l_missing_preds <- list()
  q_missing_preds <- list()
  
  for (i in 1:num_vars){
    
    print(paste0("Analysing ", study_species, ", loop: ", i))
    
    if (i == 1) {
      selected_formula <- " ~ "
      vars_selected_index <- c()
    }
    
    preds_to_analyse <- setdiff(1:length(predictors), vars_selected_index)
    l.models.aic <- c()
    q.models.aic <- c()
    for (pred in predictors[preds_to_analyse]) {
      if (i== 1){
        l_form <- formula(paste0(selected_formula, pred ))
        q_form <- formula(paste0(
          selected_formula, pred , " + I(", pred, "^2)"))
      } else{
        l_form <- formula(paste0(selected_formula, " + ", pred ))
        q_form <- formula(paste0(
          selected_formula, " + ", pred , " + I(", pred, "^2)"))
      }
      
      l.models.aic <- c(l.models.aic,
                        try(colext(psiformula = psi_formula, # First-year occupancy
                                   gammaformula = ~ 1, # Colonization
                                   epsilonformula = l_form, # Extinction
                                   pformula = det.formula, # Detection
                                   data = simUMF)@AIC, silent = T ))
      
      q.models.aic <- c(q.models.aic,
                        try(colext(psiformula = psi_formula, # First-year occupancy
                                   gammaformula = ~ 1, # Colonization
                                   epsilonformula = q_form, # Extinction
                                   pformula = det.formula, # Detection
                                   data = simUMF)@AIC, silent = T ))
      
    }
    
    ##Register missed predictors
    if(typeof(l.models.aic)=="character"){
      index_error_preds <- grep("Error", l.models.aic)
      l_missing_preds[[i]] <- index_error_preds
    } else{
      l_missing_preds[[i]] <- NA
    }
    if(typeof(q.models.aic)=="character"){
      index_error_preds <- grep("Error", q.models.aic)
      q_missing_preds[[i]] <- index_error_preds
    } else{
      q_missing_preds[[i]] <- NA
    }
    
    ##Register best predictor
    best_aic <- min(mapply(
      function(x,y) min(x,y), l.models.aic, q.models.aic))
    index_best <- which.min(mapply(
      function(x,y) min(x,y), l.models.aic, q.models.aic))
    vars_selected_index <- c(vars_selected_index, preds_to_analyse[index_best]) ##Keep track of variables selected
    quadratic_term <- l.models.aic[index_best] > q.models.aic[index_best] ##T/F
    
    ##if AIC worse than previous iteration -> break
    if (i>1){
      if (best_aic > vars.selected.AIC$AIC[i-1]){
        break
      }
    }
    
    pred <- predictors[preds_to_analyse[index_best]]
    vars.selected.AIC$iteration[i] <- i
    vars.selected.AIC$variable[i] <- pred
    vars.selected.AIC$AIC[i] <- best_aic
    vars.selected.AIC$Quadratic_term[i] <- quadratic_term
    
    
    ##Update formula
    if (i == 1){
      if (quadratic_term) {
        selected_formula <- paste0(selected_formula, pred , " + I(", pred, "^2)")
      } else{
        selected_formula <- paste0(selected_formula, pred )
      }
    }  else{ ## loop 2 to num_vars
      if (quadratic_term) {
        selected_formula <- paste0(
          selected_formula, " + ", pred , " + I(", pred, "^2)") 
      } else{
        selected_formula <- paste0(selected_formula, " + ", pred )
      }
    }
    
    
  }
  
  missing_preds <- list()
  missing_preds[["l"]] <- l_missing_preds 
  missing_preds[["q"]] <- q_missing_preds
  
  return(list(vars.selected.AIC, study_species, missing_preds) )
  
  
}

parallel::stopCluster(cl = my.cluster)


sp.names <- sapply(best_habitat_variables_ext_models, function(x) x[[2]])
missed_habitat_variables_ext_models <- lapply(best_habitat_variables_ext_models,
                                              function(x) x[[3]])
names(missed_habitat_variables_ext_models) <- sp.names
save(missed_habitat_variables_ext_models, 
     file="./colext/Analysis/Rdata/missed_habitat_variables_ext_models.Rdata")


best_habitat_variables_ext_models <- lapply(best_habitat_variables_ext_models,
                                            function(x) x[[1]])
names(best_habitat_variables_ext_models) <- sp.names
save(best_habitat_variables_ext_models, 
     file="./colext/Analysis/Rdata/best_habitat_variables_ext_models.Rdata")

#### Fit colext models with best habitat variables  ####
load("./colext/Analysis/Rdata/best_habitat_variables_col_models.Rdata")
load("./colext/Analysis/Rdata/best_habitat_variables_ext_models.Rdata")

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



parallel::detectCores()

n.cores <- parallel::detectCores() - 2

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()


clusterExport(cl=my.cluster, c('species_UMFs', 'col.equations', 'ext.equations',
                               'species.info', 'det_model_formulas',
                               'psi_formula'))



habitat_colext_models <- foreach(
  study_species = species.info$spp, 
  .packages = "unmarked"
) %dopar% {
  simUMF <- species_UMFs[[study_species]]
  
  ##Species best detection model
  det.model <- species.info$best_model[species.info$spp == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  
  gamma <- as.formula(col.equations[study_species])
  epsilon <- as.formula(ext.equations[study_species])
  
  habitat.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                               gammaformula = gamma, # Colonization
                               epsilonformula = epsilon, # Extinction
                               pformula = det.formula, # Detection
                               data = simUMF), silent = T )
  
  nP <- (habitat.model@AIC - habitat.model@negLogLike*2)/2
  
  starts_list <- list()
  starts_list[[1]] <- rep(0, nP)
  habitat.models <- list()
  
  if (class(habitat.model) == "unmarkedFitColExt"){
    habitat.models[[1]] <- habitat.model
    aic <- habitat.model@AIC
  } else{
    habitat.models[[i]] <- NA
    aic <- NA
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
      habitat.models[[i]] <- habitat.model
      aic <- c(aic, habitat.model@AIC)
    } else{
      habitat.models[[i]] <- NA
      aic <- c(aic, NA)
    }
  }
  
  index_best <- which.min(aic)
  
  habitat.model <- habitat.models[[index_best]]
  starts <- starts_list[[index_best]]
  
  
  
  return(list(habitat.model, study_species, starts))
  
  
}

parallel::stopCluster(cl = my.cluster)

sp.names <- sapply(habitat_colext_models, function(x) x[[2]])

starts_habitat_colext_models <- lapply(habitat_colext_models, function(x) x[[3]])
names(starts_habitat_colext_models) <- sp.names
save(starts_habitat_colext_models, 
     file="./colext/Analysis/Rdata/starts_habitat_colext_models.Rdata")

habitat_colext_models <- lapply(habitat_colext_models, function(x) x[[1]])
names(habitat_colext_models) <- sp.names
save(habitat_colext_models, 
     file="./colext/Analysis/Rdata/habitat_colext_models.Rdata")

