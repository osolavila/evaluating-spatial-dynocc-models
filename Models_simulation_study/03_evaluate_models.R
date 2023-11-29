
output.path <- file.path(Dir.Base, "colext", "Analysis","Rdata", "revision", "simulation")
plots.path <- file.path(Dir.Base, "colext", "Analysis","plots", "revision", "simulation")

library(Metrics) ##AUC and rmse
library(reshape2) ##Melt
library(viridis)## Color palette

#### Gather simulated data ####
load(file=file.path(output.path, paste0("surveyed.cells.array.Rdata") ) )
load(file.path(output.path, "landscape.Rdata"))
load(file.path(output.path, "occu_sim.Rdata"))
load(file=file.path(output.path, paste0("y.array.Rdata")) )

# prop.surveyed <- 0.05
# nsurveyed <- nsites*prop.surveyed
# nyears <- 2
# nsurveys <- 2
# nscenarios <- dim(y.array)[length(dim(y.array))-1]
scenario_names <- c("habitat", "brm", "habitat_brm_1", "habitat_brm_2")
study_species_vector <- scenario_names
dynamics <- c("col", "ext", "colext")

coord <- landscape_cells[,c("x","y")]*1000 ## km to m
names(coord) <- c("lon", "lat")
psi0_var <- scale(psi0_var)
site_covs <- as.data.frame(cbind(psi0_var, col_var, ext_var, coord))



species_detections <- list()
for (d in 1:length(dynamics)){
  dumf <-list()
  for (i in 1:nscenarios){
    reps_list <- list()
    for ( rep in 1:nreps){
      ysum <- apply(y.array[rep,,,,i,d], c(1,3), sum)
      ysum[ysum!=0] <- 1 
      reps_list[[rep]]<- ysum
    }
    dumf[[i]] <- reps_list
  }
  names(dumf) <- study_species_vector
  species_detections[[d]] <- dumf
}
all_species_detections <- species_detections
names(all_species_detections) <- dynamics

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

eval.reps <- function(obs, reps.occu){
  eval.list <- list()
  for ( i in 1:nreps){
    y <- obs[[i]]
    z <- reps.occu[[i]]
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
                                           function (z) eval.reps(y,z))
    } else{
      eval.list[[study_species]] <- 
        eval.reps(y, predicted.occu[[study_species]])
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
#### True models  ####
library(raster)

coord <- landscape_cells[,c("x","y")] ##km
disp.par <- 10

index_neighs <- list()
for(i in 1:nrow(coord) ){
  ##Neighbours distances
  points_distances <- pointDistance(coord[i ,],
                                    coord[,],
                                    lonlat=F)
  ##Eliminate same cell distances
  points_distances[i] <- NA
  
  
  neigh <- which(points_distances<disp.par)
  
  index_neighs[[i]] <- neigh
}

##Params
##psi0
mean.psi1 <- -2
beta.psi1 <- 3
##hab
mean.gamma <- -8
beta.gamma <- 10
mean.epsilon <- -7
beta.epsilon <- 10
##brm
gamma.params.brm <- c(-4,7)
epsilon.params.brm <- c(4,-7)
##brm habitat
gamma.params.habitat.brm.1 <- c(-8,5,5)
gamma.params.habitat.brm.2 <- c(-8,7,3)
epsilon.params.habitat.brm.1 <- c(-7,10,-5)
epsilon.params.habitat.brm.2 <- c(-7,10,-3)

psis.col <- list()
psis.ext <- list()
psis.colext <- list()

psi0 <- plogis(mean.psi1+beta.psi1*scale(psi0_var))
sp_autocovariate <- sapply(index_neighs,
                           function(x) sum(psi0[x], na.rm=T)/length(x) ) 
##habitat
gamma.habitat <- plogis(mean.gamma+beta.gamma*col_var)
epsilon.habitat <- plogis(mean.epsilon+beta.epsilon*ext_var)

psis.col[["habitat"]] <- psi0+(1-psi0)*gamma.habitat
psis.ext[["habitat"]] <- psi0*(1-epsilon.habitat)
psis.colext[["habitat"]] <- psi0*(1-epsilon.habitat)+(1-psi0)*gamma.habitat

##brm
gamma.brm <-  plogis(cbind(rep(1, nsites), sp_autocovariate)
                     %*% gamma.params.brm)
epsilon.brm <-  plogis(cbind(rep(1, nsites), sp_autocovariate)
                       %*% epsilon.params.brm)

psis.col[["brm"]] <- psi0+(1-psi0)*gamma.brm
psis.ext[["brm"]] <- psi0*(1-epsilon.brm)
psis.colext[["brm"]] <- psi0*(1-epsilon.brm)+(1-psi0)*gamma.brm

##habitat_brm_1
gamma.habitat.brm.1 <- plogis(cbind(rep(1, nsites), col_var, sp_autocovariate)
                              %*% gamma.params.habitat.brm.1)
epsilon.habitat.brm.1 <- plogis(cbind(rep(1, nsites), ext_var, sp_autocovariate)
                                %*% epsilon.params.habitat.brm.1)
psis.col[["habitat_brm_1"]] <- psi0+(1-psi0)*gamma.habitat.brm.1
psis.ext[["habitat_brm_1"]] <- psi0*(1-epsilon.habitat.brm.1)
psis.colext[["habitat_brm_1"]] <- psi0*(1-epsilon.habitat.brm.1)+(1-psi0)*gamma.habitat.brm.1


gamma.habitat.brm.2 <- plogis(cbind(rep(1, nsites), col_var, sp_autocovariate)
                              %*% gamma.params.habitat.brm.2)
epsilon.habitat.brm.2 <- plogis(cbind(rep(1, nsites), ext_var, sp_autocovariate)
                                %*% epsilon.params.habitat.brm.2)
psis.col[["habitat_brm_2"]] <- psi0+(1-psi0)*gamma.habitat.brm.2
psis.ext[["habitat_brm_2"]] <- psi0*(1-epsilon.habitat.brm.2)
psis.colext[["habitat_brm_2"]] <- psi0*(1-epsilon.habitat.brm.2)+(1-psi0)*gamma.habitat.brm.2

true_models_landscape <- list(psis.col, psis.ext, psis.colext)
names(true_models_landscape) <- dynamics

get_true_model_replicates <- function(occu){
  reps_list <- list()
  for (i in 1:nreps){
    s <- surveyed.cells.array[,i]
    reps_list[[i]] <- cbind(psi0[s], occu[s])
  }
  return(reps_list)
}

occu_true_models <- lapply(true_models_landscape, function(x) 
  lapply(x, get_true_model_replicates))


#### Load predicted occu data and model names ####
fits.files <- list.files(output.path)
fits.files <- fits.files[grep("fits", fits.files)]
for (file in fits.files){
  load(file.path(output.path, file) )
}

fits.vars.names <- sapply(fits.files, function(x) substr(x, 1, nchar(x)-6))

model.names <- fits.vars.names
model.names[grep("BRM", fits.vars.names)] <- "BRM"
model.names[grep("IFM", fits.vars.names)] <- "IFM"
model.names[grep("habitat", fits.vars.names)] <- "habitat"
model.names[grep("BRM_habitat", fits.vars.names)] <- "BRM_habitat"
model.names[grep("IFM_habitat", fits.vars.names)] <- "IFM_habitat"


# model.type <- fits.vars.names
# model.type[grep("col", fits.vars.names)] <- "col"
# model.type[grep("ext", fits.vars.names)] <- "ext"
# model.type[grep("colext", fits.vars.names)] <- "colext"

change_list_order <- function(occus){
  model <- occus
  dist_list <- list()
  for (d in names(model[[1]]) ){
    reps_list <- list()
    for (r in 1:nreps){
      reps_list[[r]] <- model[[r]][[d]]
    }
    dist_list[[d]] <- reps_list
  }
  return(dist_list)
}
fits_BRM_colext <- lapply(fits_BRM_colext, function(x) 
  lapply(x, function (y) change_list_order(y)))
fits_IFM_colext <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function (y) change_list_order(y)))
fits_BRM_habitat <- lapply(fits_BRM_habitat, function(x) 
  lapply(x, function (y) change_list_order(y)))
fits_IFM_habitat <- lapply(fits_IFM_habitat, function(x) 
  lapply(x, function (y) change_list_order(y)))

type="colext"

occu_predictions_list <- list()



for (model.occu in fits.vars.names){
  fit_m <- get(model.occu)[[type]]
  occu_predictions_list[[model.occu]] <- 
    lapply(fit_m, function(z) lapply(z, function (x){ ##loop over scenarios and over replicates
    
    if (length(x)==3){
      c <- x[[1]]
    } else {
      c <- lapply(x, function(y) y[[1]])
    }
    return(c)
  } ))
  
}
occu_predictions_list[["true_model"]] <- occu_true_models[[type]]
species_detections <- all_species_detections[[type]]


#### Spatial models evaluation  ####
base.metrics <- eval.occu.predictions(occu_predictions_list$fits_base)

metrics.names <- c("auc", "auc_change", "rmse", "rmse_change", "corr")

# load(file= file.path(output.path, "occu_BRM_colext.Rdata") )
# BRM.colext.metrics <- eval.occu.predictions(occu_BRM_colext)
# BRM.colext.metrics.vs.base <- compare_metrics(BRM.colext.metrics, base.metrics)
metrics <- list()
metrics.vs.base <- list()
for (model.occu in names(occu_predictions_list)){
  occu <- occu_predictions_list[[model.occu]]
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
  for (model.occu in names(occu_predictions_list)){
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
  for (model.occu in names(occu_predictions_list)){
    model.metrics.rel[[model.occu]] <- select.metric(
      metrics.vs.base.mean[[model.occu]], m)
    model.metrics.abs[[model.occu]] <- select.metric(
      metrics.mean[[model.occu]], m)
  }
  metrics2[[m]] <- list("abs"= model.metrics.abs, "rel"= model.metrics.rel)
}

####  Plot correlation  ####
model_names <- c(
  "base",
  "BRM_col" ,
  "BRM_ext",
  "BRM_colext" ,
  "BRM_habitat_colext",
  "habitat_colext" ,
  "IFM_col" ,
  "IFM_ext",
  "IFM_colext" ,
  "IFM_habitat_colext")[c(1,2,4,3,6,5,10,7,9,8)]

#model_colors <- viridis(10)[c(1,3,7,9,10)]
model_colors <- c("gray", turbo(20)[c(11,17)], viridis(20)[c(2,12,20)],
                  turbo(20)[c(12,19 )], viridis(20)[c(5,15)] )[c(1,2,4,3,6,5,10,7,9,8)]

index_all_models <- c(1,2,8,4,10,3,9,5,6,7)

par(mar = c(10, 4, 2, 1))

for (m in metrics.names){
  
  boxplot(best.metrics[[m]]$abs[index_all_models], 
          names=model_names[index_all_models],
          las=2, col=model_colors[index_all_models], 
          main=m)
  abline(h=median(best.metrics[[m]]$abs$fits_base) )
  
  boxplot(best.metrics[[m]]$rel[index_all_models][-1], 
          names=model_names[index_all_models][-1],
          las=2, col=model_colors[index_all_models][-1], 
          main=m)
  
  abline(h=0)
}

#### AIC  ####

aic_predictions_list <- list()


for (model.occu in fits.vars.names){
  
  aic_predictions_list[[model.occu]] <- lapply(get(model.occu), function(x) {
    
    if (length(x)==3){
      c <- x[[3]]
    } else {
      c <- lapply(x, function(y) y[[3]])
    }
    return(c)
  } )
  
  
  
}



compare_aic <- function(metrics1, metrics2){
  eval.list <- list()
  for (study_species in study_species_vector){
    y1 <- metrics1[[study_species]]
    y2 <- metrics2[[study_species]]
    if (typeof(y1) =="list" & typeof(y2) =="list"){
      eval.list[[study_species]] <- mapply(function(x,z) 1-x/z, y1, y2)
    } else if (typeof(y1) =="list" & typeof(y2) =="double") {
      eval.list[[study_species]] <- lapply(y1, function(x) 1-x/y2)
      
    } else if (typeof(y1) =="double" & typeof(y2) =="double") {
      eval.list[[study_species]] <- 1-y1/y2
    } 
  }
  return(eval.list)
}

select.best.model.aic <- function(model.metrics){
  eval.list <- vector("numeric", length(study_species_vector))
  names(eval.list) <- study_species_vector
  for (study_species in study_species_vector){
    y <- model.metrics[[study_species]]
    if (typeof(y)=="list" ){
      metrics.distances <- do.call(rbind, y)
      eval.list[study_species] <- max(metrics.distances)
    } else{
      eval.list[study_species] <- max(y)
    }
  }
  return(eval.list)
}

aic.vs.base <- list()
best.aic.vs.base <- list()
for (model.aic in names(aic_predictions_list)){
  aic <- aic_predictions_list[[model.aic]]
  
  aic.vs.base[[model.aic]] <- compare_aic(
    aic, aic_predictions_list$fits_base)
  best.aic.vs.base[[model.aic]] <- select.best.model.aic(aic.vs.base[[model.aic]])
}
best.aic.vs.base <- as.data.frame( do.call(cbind, best.aic.vs.base) )


boxplot(best.aic.vs.base[index_all_models][-1], 
        names=model_names[index_all_models][-1],
        las=2, col=model_colors[index_all_models][-1], 
        main="AIC_gain")

abline(h=0)


#### Plot boxplots  ####
metric <- "auc"
for (i in 1:length(study_species_vector)){
  df.abs <- data.frame(matrix(nrow = nreps,
                          ncol=length(names(occu_predictions_list))))
  names(df.abs) <- names(occu_predictions_list)
  
  df.rel <- data.frame(matrix(nrow = nreps,
                              ncol=length(names(occu_predictions_list))-1))
  names(df.rel) <- names(occu_predictions_list)[-1]
  for (model.occu in names(occu_predictions_list)){
    
    if (grepl("BRM", model.occu)){
      df.abs[,model.occu] <- metrics[[model.occu]][[i]][["max_dist_10"]][metric,]
    } else if (grepl("IFM", model.occu)){
      df.abs[,model.occu] <- metrics[[model.occu]][[i]][["10"]][metric,]
    } else{
      df.abs[,model.occu] <- metrics[[model.occu]][[i]][metric,]
    }
  }
}


#### Plot AUC - d/alpha colext models (Figure 2) ####
buffer_distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
inf_fun_distances_analysed <- c(1,2.5,5,10,20)


#model_colors <- viridis(10)[c(1,3,7,9,10)]
model_colors <- viridis(20)[c(2,5,12,15,20)]

metrics.list <- metrics2$auc$rel

max_metrics <- max(unlist(metrics.list))
diff.scale <- T
#specify path to save PDF to
# destination <- file.path(plots.path, "auc_gain.pdf")
# #open PDF
# pdf(file=destination, width = 8.3, height = 11.7)

# layout(matrix(c(1:10), 5, 2, byrow = TRUE), widths = lcm(c(9,9)),
       # heights = lcm(rep(5,5)))

#adjust plot margins
par(mar = c(4, 4, 2, 1))

for (i in 1:length(study_species_vector)){
  
  if (diff.scale){
    max_metrics <- sapply(metrics.list, function(x) x[[study_species_vector[[i]]]])
    max_metrics <- max(unlist(max_metrics))
  }
  
  ##Spatial buffer
  plot( buffer_distances_analysed, metrics.list$fits_BRM_colext[[i]],
        xlim = c(1,max(inf_fun_distances_analysed)) ,
        # xlim = c(1,max(buffer_distances_analysed)) ,
        # ylim = c(0,max(meta_col_diff_aic, auto_col_diff_aic)),
        ylim = c(0,max_metrics*1.1),
        col=model_colors[1], 
        # main = substitute(italic(x), list(x=study_species_vector[i])),
        type = "l", lty=2, lwd=1,
        xlab="",
        ylab= "",
        cex.axis=1.2)
  points(buffer_distances_analysed, metrics.list$fits_BRM_colext[[i]],
         pch= 22, cex=1.5, bg=model_colors[1], col="white", lwd=0.5)
  title(xlab=expression(paste("d or ",alpha, " (km)")), ylab= expression(paste(Delta,"AUC") ),
        line=2.2, cex.lab=1.4)
  
  ##Spatial inference function
  lines(inf_fun_distances_analysed, metrics.list$fits_IFM_colext[[i]],
        lty=2, lwd=1, col=model_colors[2])
  points(inf_fun_distances_analysed, metrics.list$fits_IFM_colext[[i]],
         pch= 23, cex=1.5, bg=model_colors[2], col="white", lwd=0.5)
  
  ##Habitat
  lines(buffer_distances_analysed, 
        rep(metrics.list$fits_habitat[[i]],
            length(buffer_distances_analysed)),
        lty=1, lwd=3, col=model_colors[5])
  
  ##Habitat+Spatial brm
  lines(buffer_distances_analysed,
        metrics.list$fits_BRM_habitat[[i]],
        lty=2, lwd=1, col=model_colors[3])      
  points(buffer_distances_analysed,
         metrics.list$fits_BRM_habitat[[i]],
         pch= 22, cex=1.5, bg=model_colors[3], col="white", lwd=0.5)
  
  ##Habitat+Spatial inference function
  lines(inf_fun_distances_analysed,
        metrics.list$fits_IFM_habitat[[i]],
        lty=2, lwd=1, col=model_colors[4])
  points(inf_fun_distances_analysed, 
         metrics.list$fits_IFM_habitat[[i]],
         pch= 23, cex=1.5, bg=model_colors[4], col="white", lwd=0.5)
  
  # legend("topleft", col = c("blue", "blue", "yellow", "green", "green"),
  #        legend = c("spatial_buffer", "spatial_inf_fun", "habitat",
  #                   "habitat_spatial_buffer", "habitat_spatial_inf_fun"), 
  #        ncol = 2, 
  #        lty=c(1,2,1,1,2), lwd=3)
  
}

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c("spatial BRM colext", "spatial IFM colext", "habitat colext",
                           "habitat-spatial BRM colext", "habitat-spatial IFM colext"),
       pch=c(15,18, NA, 15,18), lty = c(NA,NA,1,NA,NA),
       pt.cex=2, cex=1, bty='n', lwd= 3,
       col = model_colors[c(1,2,5,3,4)])
mtext("Legend", line=-1.5, at=0.5, cex=1)
dev.off()

#### Results  ####

##best AUC increase
for (i in 1:4){
  max_metrics <- sapply(metrics.list, function(x) x[[study_species_vector[[i]]]])
  max_metrics <- max(unlist(max_metrics))
  print(max_metrics)
}

##habitat vs spatial
best.metrics$auc$rel$fits_BRM_colext / best.metrics$auc$rel$fits_habitat

##habitat-brm vs habitat
best.metrics$auc$rel$fits_BRM_habitat / best.metrics$auc$rel$fits_habitat

##habitat-brm vs brm
best.metrics$auc$rel$fits_BRM_habitat / best.metrics$auc$rel$fits_BRM_colext

##habitat-brm vs habitat vs brm
best.metrics$auc$rel$fits_BRM_habitat /
  (best.metrics$auc$rel$fits_habitat + best.metrics$auc$rel$fits_BRM_colext)

(best.metrics$auc$rel$fits_BRM_habitat - pmax(best.metrics$auc$rel$fits_habitat, best.metrics$auc$rel$fits_BRM_colext) )/
  pmin(best.metrics$auc$rel$fits_habitat, best.metrics$auc$rel$fits_BRM_colext)




best.metrics$auc_change$rel$fits_BRM_habitat / best.metrics$auc_change$rel$fits_habitat

best.metrics$auc_change$rel$fits_BRM_habitat /
  (best.metrics$auc_change$rel$fits_habitat + best.metrics$auc_change$rel$fits_BRM_colext)

(best.metrics$auc_change$rel$fits_BRM_habitat - pmax(best.metrics$auc_change$rel$fits_habitat, best.metrics$auc_change$rel$fits_BRM_colext) )/
  pmin(best.metrics$auc_change$rel$fits_habitat, best.metrics$auc_change$rel$fits_BRM_colext)