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


##Clean species occurrence data
species.occurance.change <- db_aonc_change %>% 
  filter(cell_1x1 %in% cells_study$cell_1x1,
         spp %in% species.groups$spp) %>%
  group_by(spp, .drop=FALSE) %>%
  summarise( aonc2 = length(which(Canvi==0)) + length(which(Canvi==-1)),
             col= length(which(Canvi==1)),
             ext= length(which(Canvi==-1)),
  ) %>%
  # filter(spp %in% species.groups$spp) %>%
  mutate(perc.change= round((col+aonc2-ext)/(aonc2)*100-100,1)) %>% 
  filter(aonc2>100)




####  Fit colext model with different detection models####
intercept <- "~ 1"
atlas <- "~ atlas"
quadratic_date <- "~ julian_day_std + I(julian_day_std^2)"
quadratic_date_atlas <- "~ atlas + julian_day_std + I(julian_day_std^2)"
interaction_atlas_quadratic_date <- "~ atlas + julian_day_std + 
I(julian_day_std^2) + atlas:julian_day_std + atlas:I(julian_day_std^2)"


det_model_formulas <- c(intercept, atlas, quadratic_date, quadratic_date_atlas, 
                        interaction_atlas_quadratic_date)
det_model_names <- c("intercept", "atlas", "quadratic_date",
                     "quadratic_date_atlas",
                     "interaction_atlas_quadratic_date")
# model_formulas <- c(quadratic_date,
#                     paste(quadratic_date, "+ elev_std"),
#                     paste(quadratic_date, "+ lat_std"),
#                     paste(quadratic_date, "+ elev_std + lat_std"))


study_species_list <- species.occurance.change$spp

models <- list()
study_species <- "Turdus viscivorus"
for (study_species in study_species_list){
  
  simUMF <- get_UMF(study_species)
  # summary(simUMF)
  print(paste0("Analysing species: ", study_species))
  
  #### run models ####
  
  psi_formula <- formula(" ~ aonc2_sdm")
  
  det_models <- list()
  
  for (i in 1:length(det_model_formulas)){
    m1 <- try(colext(psiformula = psi_formula, # First-year occupancy
                     gammaformula = ~ 1, # Colonization
                     epsilonformula = ~ 1, # Extinction
                     pformula = as.formula(det_model_formulas[i]), # Detection
                     data = simUMF) )
    
    det_models[[det_model_names[i]]] <- m1
  }
  
  models[[study_species]] <- det_models
}

#### Analys best detection model ####

##DF: rows -> species, columns -> det models AIC
best_det_model <- data.frame(matrix(nrow= length(study_species_list),
                                    ncol= length(det_model_formulas)+2 ))

colnames(best_det_model) <- c("study_species", "best_model",
                              det_model_names)
rownames(best_det_model) <- study_species_list

best_det_model$study_species <- study_species_list

for (study_species in study_species_list){
  for (det_m in det_model_names){
    best_det_model[study_species, det_m] <- models[[study_species]][[det_m]]@AIC
  }
  best_det_model[study_species, "best_model"] <- det_model_names[
    which.min(best_det_model[study_species, det_model_names])]
}

table(best_det_model$best_model)

#### Mean detection rates per atlas period ####
det_rates <- data.frame(matrix(nrow= length(study_species_list),
                                    ncol= 4 ))

colnames(det_rates) <- c("study_species", "intercept",
                              "Atlas_2", "Atlas_3")
rownames(det_rates) <- study_species_list

det_rates$study_species <- study_species_list

for (study_species in study_species_list){
  
  det_rates[study_species, "intercept"] <- unmarked::predict(
    models[[study_species]][[1]], type="det")$Predicted[1]
  
  det_rates[study_species, c("Atlas_2", "Atlas_3")] <- unmarked::predict(
    models[[study_species]][[2]], type="det")$Predicted[c(1,3)]
  
}
length(which(det_rates$intercept>0.5))

length(which(det_rates$Atlas_2>0.5 & det_rates$Atlas_3>0.5))

table(best_det_model$best_model[
  which(det_rates$Atlas_2>0.5 & det_rates$Atlas_3>0.5)])
#### Join detection data and save as Rdata  ####

species_detection_models_info <- det_rates %>% left_join(best_det_model[,1:2],
                                             by="study_species")

save(species_detection_models_info,
     file="./colext/Analysis/Rdata/species_detection_models_info.Rdata")
load("./colext/Analysis/Rdata/species_detection_models_info.Rdata")
#### Species selection according detection models ####
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
# model_formulas <- c(quadratic_date,
#                     paste(quadratic_date, "+ elev_std"),
#                     paste(quadratic_date, "+ lat_std"),
#                     paste(quadratic_date, "+ elev_std + lat_std"))


study_species_list <- species.occurance.change$spp

det_models <- list()

study_species <- "Turdus viscivorus"
for (study_species in study_species_list){
  
  simUMF <- get_UMF(study_species)
  # summary(simUMF)
  print(paste0("Analysing species: ", study_species))
  
  #### run models ####
  
  det.model <- species_detection_models_info$best_model[
    species_detection_models_info$study_species == study_species]
  det.formula <- as.formula( det_model_formulas[det.model] )
  

  
  
  m1 <- try(colext(psiformula = psi_formula, # First-year occupancy
                   gammaformula = ~ 1, # Colonization
                   epsilonformula = ~ 1, # Extinction
                   pformula = det.formula, # Detection
                   data = simUMF) )
    
    
  
  
  det_models[[study_species]] <- m1
}

save(det_models,
     file="./colext/Analysis/Rdata/det_models.Rdata")


load("./colext/Analysis/Rdata/det_models.Rdata")

det.predictions <- list()
for (study_species in study_species_list){
  det.predictions[[study_species]] <-
    predict(det_models[[study_species]], type="det")[,1]
}


perc.cells.above.x.detectability <- function(x, det.predictions){
  perc.df <- data.frame( matrix(ncol=3, nrow=length(study_species_list)) )
  colnames(perc.df) <- c("spp", "perc.aonc2", "perc.aonc3")
  
  i<- 1
  for (study_species in study_species_list) {
   det.mat <- matrix( det.predictions[[study_species]], ncol=4, byrow = T)
   perc.df$spp[i] <- study_species
   perc.df$perc.aonc2[i] <- length(which(apply(det.mat[,1:2], 1, max)>x))/
     nrow(det.mat)
   perc.df$perc.aonc3[i] <- length(which(apply(det.mat[,3:4], 1, max)>x))/
     nrow(det.mat)
   i<-i+1
  }
  
  return(perc.df)                        
}

a <- perc.cells.above.x.detectability(0.7, det.predictions)
#### Study correlation between %colext and detection rate ####
load("./colext/Analysis/Rdata/species_detection_models_info.Rdata")

perc.col <- species.occurance.change$col / species.occurance.change$aonc2
perc.ext <- species.occurance.change$ext / species.occurance.change$aonc2
perc.colext  <- (species.occurance.change$col + species.occurance.change$ext) /
  species.occurance.change$aonc2

corr_det_rate_obs_colext <- as.data.frame(
  cbind(perc.col, perc.ext, perc.colext, det_rates$intercept))
colnames(corr_det_rate_obs_colext)[4] <- "detection.rate"

par(mar=c(5, 5, 4, 2) + 0.1)
##Col
plot(det_rates$intercept, perc.col, xlab = "Detection rate", 
     ylab = "% Colonized cells", 
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.col ~ detection.rate, data = corr_det_rate_obs_colext),
       col = "blue")
##Ext
plot(det_rates$intercept, perc.ext, xlab = "Detection rate", 
     ylab = "% Extinct cells",
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.ext ~ detection.rate, data = corr_det_rate_obs_colext),
       col = "blue")

##Ext+col
plot(det_rates$intercept, perc.colext)

par(mar=c(5, 4, 4, 2) + 0.1)

##Different detection rates between atles (atles 3 – atlas 2)
##is predictor of colonizations / extinctions?

corr_diff_det_rate_obs_colext <- as.data.frame(
  cbind(perc.col, perc.ext, 
        perc.col - perc.ext,
        species_detection_models_info$Atlas_3 - 
          species_detection_models_info$Atlas_2))

colnames(corr_diff_det_rate_obs_colext)[3:4] <- c("diff_colext",
                                                  "diff_det_rate")

plot(corr_diff_det_rate_obs_colext$diff_det_rate,
     corr_diff_det_rate_obs_colext$diff_colext, 
     xlab = "Diff fitted detection rate (AONC3-AONC2)", 
     ylab = "% Col - % Ext cells"
     )
abline(lm(diff_colext ~ diff_det_rate, data = corr_diff_det_rate_obs_colext),
       col = "blue")



#### Observed detection rates only presence  ####

obs_det_rates <- data.frame(matrix(nrow= length(study_species_list),
                                   ncol= 4 ))

colnames(obs_det_rates) <- c("study_species", "intercept",
                             "Atlas_2", "Atlas_3")
rownames(obs_det_rates) <- study_species_list

obs_det_rates$study_species <- study_species_list
for (study_species in study_species_list){
  
  simUMF <- get_UMF(study_species)
  # summary(simUMF)
  print(paste0("Analysing species: ", study_species))
  
  
  y <- simUMF@y
  detection_cells <- which(apply(y,1,sum)>=1)
  detection_cells_aonc2 <- which(apply(y[,1:2],1,sum)>=1)
  detection_cells_aonc3 <- which(apply(y[,3:4],1,sum)>=1)
 
  obs_det_rates[study_species, "intercept"] <- sum(y[detection_cells,])/
    length(y[detection_cells,])
  obs_det_rates[study_species, "Atlas_2"] <- sum(y[detection_cells_aonc2,1:2])/
    length(y[detection_cells_aonc2,1:2])
  obs_det_rates[study_species, "Atlas_3"] <- sum(y[detection_cells_aonc3,3:4])/
    length(y[detection_cells_aonc3,3:4])
    
  
}

perc.col <- species.occurance.change$col / species.occurance.change$aonc2
perc.ext <- species.occurance.change$ext / species.occurance.change$aonc2

corr_det_rate_obs_colext <- as.data.frame(
  cbind(perc.col, perc.ext,  obs_det_rates))


par(mar=c(5, 5, 4, 2) + 0.1)
##Col
plot(corr_det_rate_obs_colext$Atlas_2,
     corr_det_rate_obs_colext$perc.col,
     xlab = "Obs detection rate surveys with presence AONC2", 
     ylab = "% Colonized cells", 
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.col ~ Atlas_2, data = corr_det_rate_obs_colext),
       col = "blue")

plot(corr_det_rate_obs_colext$Atlas_3,
     corr_det_rate_obs_colext$perc.col,
     xlab = "Obs detection rate surveys with presence AONC3", 
     ylab = "% Colonized cells", 
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.col ~ Atlas_3, data = corr_det_rate_obs_colext),
       col = "blue")


##Ext
plot(corr_det_rate_obs_colext$Atlas_2,
     corr_det_rate_obs_colext$perc.ext,
     xlab = "Obs detection rate surveys with presence AONC2", 
     ylab = "% Extinct cells", 
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.ext ~ Atlas_2, data = corr_det_rate_obs_colext),
       col = "blue")

plot(corr_det_rate_obs_colext$Atlas_3,
     corr_det_rate_obs_colext$perc.ext,
     xlab = "Obs detection rate surveys with presence AONC3", 
     ylab = "% Extinct cells", 
     ylim= c(0,max(corr_det_rate_obs_colext[,1:2])))
abline(lm(perc.ext ~ Atlas_3, data = corr_det_rate_obs_colext),
       col = "blue")


##Fitted vs observed correlation rate
plot(species_detection_models_info$Atlas_2,
     corr_det_rate_obs_colext$Atlas_2,
     xlab = "Fitted detection rate AONC2", 
     ylab = "Obs detection rate surveys with presence AONC2", 
     ylim= c(0.5,1))

atlas3_fitted_vs_obs <- as.data.frame(cbind(
  species_detection_models_info$Atlas_3,
  corr_det_rate_obs_colext$Atlas_3))
colnames(atlas3_fitted_vs_obs) <- c("fitted", "obs")


plot(species_detection_models_info$Atlas_3,
    corr_det_rate_obs_colext$Atlas_3,
    xlab = "Fitted detection rate AONC3", 
    ylab = "Obs detection rate surveys with presence AONC3", 
    ylim= c(0.5,1))
abline(lm(obs ~ fitted, data = atlas3_fitted_vs_obs),
       col = "blue")


plot(species_detection_models_info$Atlas_3 - 
       species_detection_models_info$Atlas_2,
     corr_det_rate_obs_colext$Atlas_3 -
       corr_det_rate_obs_colext$Atlas_2,
     xlab = "Fitted detection rate AONC3-AONC2", 
     ylab = "Obs detection AONC3-AONC2", 
)


##Analyse outliers
atlas3_fitted_vs_obs$lm <-
  lm(obs ~ fitted, data = atlas3_fitted_vs_obs)$fitted.values
atlas3_fitted_vs_obs$res <- atlas3_fitted_vs_obs$lm - atlas3_fitted_vs_obs$obs

atlas3_fitted_vs_obs$perc.ext <- corr_det_rate_obs_colext$perc.ext
atlas3_fitted_vs_obs$perc.col <- corr_det_rate_obs_colext$perc.col


plot(atlas3_fitted_vs_obs$res,
     atlas3_fitted_vs_obs$perc.ext,
     xlab = "Residual fitted vs observed", 
     ylab = "% Extinct cells", 
     ylim= c(0,1))

plot(atlas3_fitted_vs_obs$res,
     atlas3_fitted_vs_obs$perc.col,
     xlab = "Residual fitted vs observed", 
     ylab = "% Colonized cells", 
     ylim= c(0,1))

plot(atlas3_fitted_vs_obs$res,
     atlas3_fitted_vs_obs$perc.col-atlas3_fitted_vs_obs$perc.ext,
     xlab = "Residual fitted vs observed", 
     ylab = "% Ext/col cells", 
     )

##Different detection rates between atles (atles 3 – atlas 2)
##is predictor of colonizations / extinctions?

corr_diff_det_rate_obs_colext <- as.data.frame(
  cbind(perc.col, perc.ext, 
        perc.col - perc.ext,
        obs_det_rates$Atlas_3 - 
          obs_det_rates$Atlas_2))

colnames(corr_diff_det_rate_obs_colext)[3:4] <- c("diff_colext",
                                                  "diff_det_rate")

plot(corr_diff_det_rate_obs_colext$diff_det_rate,
     corr_diff_det_rate_obs_colext$diff_colext, 
     xlab = "Diff obs detection rate (AONC3-AONC2)", 
     ylab = "% Col - % Ext cells"
)
abline(lm(diff_colext ~ diff_det_rate, data = corr_diff_det_rate_obs_colext),
       col = "blue")
#### Independent analysis detection models with only presence ####
glm_models <- list()
glm_det_rates <- data.frame(matrix(nrow= length(study_species_list),
                               ncol= 4 ))

colnames(glm_det_rates) <- c("study_species", "intercept",
                         "Atlas_2", "Atlas_3")
rownames(glm_det_rates) <- study_species_list

glm_det_rates$study_species <- study_species_list
for (study_species in study_species_list){
  
  simUMF <- get_UMF(study_species)
  # summary(simUMF)
  print(paste0("Analysing species: ", study_species))
  
  
  y <- simUMF@y
  detection_cells <- which(apply(y,1,sum)>=1)
  y <- as.vector(t(y[detection_cells,]))
  julian_day_std <- matrix(simUMF@obsCovs$julian_day_std,ncol=4, byrow = T)
  julian_day_std <- as.vector(t(julian_day_std[detection_cells,]))
  atlas <- matrix(
    simUMF@yearlySiteCovs$atlas[rep(1:nrow(simUMF@yearlySiteCovs), each = 2)],
    ncol=4, byrow=T)
  atlas <- as.vector(t(atlas[detection_cells,]))
  
  det_data <- data.frame(y=y, julian_day_std=julian_day_std, atlas=atlas)
  det_models <- list()
  for (i in 1:length(det_model_formulas)){
    log_reg <- glm(formula = as.formula(paste0("y",det_model_formulas[i])),
                   family=binomial(link = "logit"), data=det_data)
    
    
    det_models[[det_model_names[i]]] <- log_reg
  }
  glm_det_rates[study_species, "intercept"] <- predict(
    det_models[[1]], type="response")[1]
  glm_det_rates[study_species, c("Atlas_2", "Atlas_3")] <- predict(
    det_models[[2]], type="response")[c(1,3)]
  glm_models[[study_species]] <- det_models
}

##DF: rows -> species, columns -> det models AIC
best_glm_model <- data.frame(matrix(nrow= length(study_species_list),
                                    ncol= length(det_model_formulas)+2 ))

colnames(best_glm_model) <- c("study_species", "best_model",
                              det_model_names)
rownames(best_glm_model) <- study_species_list

best_glm_model$study_species <- study_species_list

for (study_species in study_species_list){
  for (det_m in det_model_names){
    best_glm_model[study_species, det_m] <- glm_models[[study_species]][[det_m]]$aic
  }
  best_glm_model[study_species, "best_model"] <- det_model_names[
    which.min(best_glm_model[study_species, det_model_names])]
}

##For how many species best model coincides
length(which(best_det_model$best_model==best_glm_model$best_model))
##35/46 agreement
table(best_det_model$best_model, best_glm_model$best_model)
## 4/11: glm -> intercept, colext <- other
## 3/11 colext <- interaction, glm <- date_atlas
## 4/11: colext -> date_atlas glm -> 2 atlas, 2 interaction

##Difference between colext estimated det rates and glm
boxplot(det_rates$intercept - glm_det_rates$intercept)
