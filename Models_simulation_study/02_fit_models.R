

output.path <- file.path(Dir.Base, "colext", "Analysis", "Rdata",
                         "revision", "simulation")


library(foreach)
library(doParallel)
library(unmarked)


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
dynamics <- c("col", "ext", "colext")

coord <- landscape_cells[,c("x","y")]*1000 ## km to m
names(coord) <- c("lon", "lat")
psi0_var <- scale(psi0_var)
site_covs <- as.data.frame(cbind(psi0_var, col_var, ext_var, coord))
##Create simUMFs
simUMFs <- list()

for (d in 1:length(dynamics)){
  dumf <-list()
  for (i in 1:nscenarios){
    reps_list <- list()
    for ( rep in 1:nreps){
      detection_data <- matrix(NA, nsites, nyears*nsurveys)
      detection_data[surveyed.cells.array[,rep],]<- y.array[rep,,,,i,d]
      
      reps_list[[rep]] <- unmarkedMultFrame(
        y = detection_data,
        siteCovs = site_covs,
        numPrimary=2)
    }
    dumf[[i]] <- reps_list
  }
  simUMFs[[d]] <- dumf
}

psi_formula <- ~ psi0_var

#### Spatial colext functions ####

source(file.path(Dir.Base, "colext", "Analysis","pipeline", "revision",
                 "simulation", "2_seasons_spatial_colext-v4.R") )

#### Fit base ####

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
  d=1:length(dynamics), 
  .packages = "unmarked"
) %:%  foreach(sim_i = 1:nscenarios
               ) %:%  foreach(rep = 1:nreps)  %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  simUMF0 <- simUMFs[[d]][[sim_i]][[rep]]
  
  ##fit models
  spatial.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                              gammaformula = ~ 1, # Colonization
                              epsilonformula = ~ 1, # Extinction
                              pformula = ~ 1, # Detection
                              data = simUMF0) )
  
  model.list <- list()
  if (class(spatial.model) == "unmarkedFitColExt"){
    model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
    model.list[[2]] <- spatial.model@estimates
    model.list[[3]] <- spatial.model@AIC
    
  } 
  
  return(model.list)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_base) <- dynamics
for (d in dynamics){
  names(fits_base[[d]]) <- scenario_names
}

save(fits_base,
     file=file.path(output.path, "fits_base.Rdata"))

#### Fit habitat ####
n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)


fits_habitat <- foreach(
  d=1:length(dynamics), 
  .packages = "unmarked"
) %:%  foreach(sim_i = 1:nscenarios
               )  %:%  foreach(rep = 1:nreps) %dopar% {
  
  simUMF <- simUMFs[[d]][[sim_i]][[rep]]
  
  
  gamma <- ~col_var
  epsilon <- ~ext_var
  
  habitat.model <- try(colext(psiformula = psi_formula, # First-year occupancy
                              gammaformula = gamma, # Colonization
                              epsilonformula = epsilon, # Extinction
                              pformula = ~1, # Detection
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
             pformula = ~1, # Detection
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
  
  
  
  return(list(habitat.model, starts))
  
  
}

parallel::stopCluster(cl = my.cluster)




starts_habitat <- lapply(fits_habitat, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )

names(starts_habitat) <- dynamics
for (d in dynamics){
  names(starts_habitat[[d]]) <- scenario_names
}
save(starts_habitat, 
     file=file.path(output.path, "starts_habitat.Rdata"))


fits_habitat <- lapply(fits_habitat, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

names(fits_habitat) <- dynamics
for (d in dynamics){
  names(fits_habitat[[d]]) <- scenario_names
}
save(fits_habitat, 
     file=file.path(output.path, "fits_habitat.Rdata"))





#### Fit BRM ####
distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
names_dist <- paste0("max_dist_", c(1:3, seq(5,10,2.5), seq(15,30,5)) )

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
clusterExport(cl=my.cluster, c('spatial.colext', 'spatial.colext.fit'))

writeLines(c(""), "log.txt")

fits_BRM_colext <- foreach(
  d=1:length(dynamics), 
  .packages = c("unmarked", "raster")
) %:%  foreach(sim_i = 1:nscenarios
               )  %:%  foreach(rep = 1:nreps) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  ##Spatial.colext functions have to use unmarked environment
  environment(spatial.colext) <- asNamespace('unmarked')
  environment(spatial.colext.fit) <- asNamespace('unmarked')
  
  simUMF0 <- simUMFs[[d]][[sim_i]][[rep]]
  
  dist.models <- list()
  for (dist_i in 1:length(distances_analysed)) {
    
    
    disp.par <- distances_analysed[dist_i] + 0.1
    dist <- names_dist[dist_i]
    
    spatial.model <- try( spatial.colext(psiformula = psi_formula, # First-year occupancy
                                         data = simUMF0,
                                         spatial.col = T, 
                                         spatial.ext = T, sp.model= "BRM",
                                         disp.par= disp.par) )
    
    model.list <- list()
    if (class(spatial.model) == "unmarkedFitColExt"){
      model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
      model.list[[2]] <- spatial.model@estimates
      model.list[[3]] <- spatial.model@AIC
      dist.models[[dist]] <- model.list
    } else{
      dist.models[[dist]] <- NA
    }
    
  }
  
  return(dist.models)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_BRM_colext) <- dynamics
for (d in dynamics){
  names(fits_BRM_colext[[d]]) <- scenario_names
}

save(fits_BRM_colext,
     file=file.path(output.path, "fits_BRM_colext.Rdata"))


#### Fit BRM habitat ####
distances_analysed <- c(1:3, seq(5,10,2.5), seq(15,30,5))
names_dist <- paste0("max_dist_", c(1:3, seq(5,10,2.5), seq(15,30,5)) )

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
clusterExport(cl=my.cluster, c('spatial.colext', 'spatial.colext.fit'))

writeLines(c(""), "log.txt")

fits_BRM_habitat <- foreach(
  d=1:length(dynamics), 
  .packages = c("unmarked", "raster")
) %:%  foreach(sim_i = 1:nscenarios
               )  %:%  foreach(rep = 1:nreps) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  ##Spatial.colext functions have to use unmarked environment
  environment(spatial.colext) <- asNamespace('unmarked')
  environment(spatial.colext.fit) <- asNamespace('unmarked')
  
  simUMF0 <- simUMFs[[d]][[sim_i]][[rep]]
  
  dist.models <- list()
  for (dist_i in 1:length(distances_analysed)) {
    
    
    disp.par <- distances_analysed[dist_i] + 0.1
    dist <- names_dist[dist_i]
    
    spatial.model <- try(spatial.colext(psiformula = psi_formula, # First-year occupancy
                                        gammaformula = ~col_var,
                                        epsilonformula = ~ext_var,
                                        data = simUMF0,
                                        spatial.col = T, 
                                        spatial.ext = T, sp.model= "BRM",
                                        disp.par= disp.par) )
    
    model.list <- list()
    if (class(spatial.model) == "unmarkedFitColExt"){
      model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
      model.list[[2]] <- spatial.model@estimates
      model.list[[3]] <- spatial.model@AIC
      dist.models[[dist]] <- model.list
    } else{
      dist.models[[dist]] <- NA
    }
    
  }
  
  return(dist.models)
  
}

parallel::stopCluster(cl = my.cluster)

names(fits_BRM_habitat) <- dynamics
for (d in dynamics){
  names(fits_BRM_habitat[[d]]) <- scenario_names
}

save(fits_BRM_habitat,
     file=file.path(output.path, "fits_BRM_habitat.Rdata"))

#### Fits IFM ####

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
clusterExport(cl=my.cluster, c('spatial.colext', 'spatial.colext.fit'))

writeLines(c(""), "log.txt")

fits_IFM_colext <- foreach(
  d=1:length(dynamics), 
  .packages = c("unmarked", "raster", "IFMcpp")
) %:%  foreach(sim_i = 1:nscenarios
               )  %:%  foreach(rep = 1:nreps) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  ##Spatial.colext functions have to use unmarked environment
  environment(spatial.colext) <- asNamespace('unmarked')
  environment(spatial.colext.fit) <- asNamespace('unmarked')
  
  simUMF0 <- simUMFs[[d]][[sim_i]][[rep]]
  
  dist.models <- list()
  divs <- list()
  for (alpha in c(1,2.5,5,10,20) ){
    
    cat(paste("Analysing scenario: ", sim_i, " alpha: ", alpha,
              "\n"), file="log.txt", append=TRUE)
    
    disp.par <- alpha
    dist <- as.character(alpha)
    
    best_model <- list()
    
    for (div in c(1, 2.5, 5, 10, 25, 50, 100, 200)){
      
      spatial.model <- try( spatial.colext(psiformula = psi_formula, # First-year occupancy
                                           data = simUMF0,
                                           spatial.col = T, 
                                           spatial.ext = T, sp.model= "IFM",
                                           disp.par= disp.par, 
                                           IFM.bounding.box = 30,
                                           div = div) )
      
      model.list <- list()
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        best_model[[as.character(div)]] <- model.list
      } else{
        best_model[[as.character(div)]] <- list(NA, NA, NA)
      }
      
    }
    
    index_best <- which.min(sapply(best_model, function(x) x[[3]] ) )
    
    dist.models[[dist]] <-  best_model[[index_best]]
    divs[[dist]] <- index_best
    
  }
  return(list(dist.models, divs) )
  
}

parallel::stopCluster(cl = my.cluster)

best_div_ifm <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )
  
fits_IFM_colext <- lapply(fits_IFM_colext, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )

names(fits_IFM_colext) <- dynamics
for (d in dynamics){
  names(fits_IFM_colext[[d]]) <- scenario_names
}

names(best_div_ifm) <- dynamics
for (d in dynamics){
  names(best_div_ifm[[d]]) <- scenario_names
}

save(fits_IFM_colext,
     file=file.path(output.path, "fits_IFM_colext.Rdata"))

save(best_div_ifm,
     file=file.path(output.path, "best_div_ifm.Rdata"))

#### Fits IFM habitat ####

n.cores <- parallel::detectCores() - 4

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
clusterExport(cl=my.cluster, c('spatial.colext', 'spatial.colext.fit'))

writeLines(c(""), "log.txt")

fits_IFM_habitat <- foreach(
  d=1:length(dynamics), 
  .packages = c("unmarked", "raster", "IFMcpp")
) %:%  foreach(sim_i = 1:nscenarios
               )  %:%  foreach(rep = 1:nreps) %dopar% {
  
  # sink("log.txt", append=TRUE)
  
  ##Spatial.colext functions have to use unmarked environment
  environment(spatial.colext) <- asNamespace('unmarked')
  environment(spatial.colext.fit) <- asNamespace('unmarked')
  
  simUMF0 <- simUMFs[[d]][[sim_i]][[rep]]
  divs <- list()
  
  dist.models <- list()
  for (alpha in c(1,2.5,5,10,20) ){
    
    
    disp.par <- alpha
    dist <- as.character(alpha)
    
    best_model <- list()
    
    for (div in c(1, 2.5, 5, 10, 25, 50, 100, 200)){
      
      spatial.model <- try( spatial.colext(psiformula = psi_formula, # First-year occupancy
                                           gammaformula = ~col_var,
                                           epsilonformula = ~ext_var,
                                           data = simUMF0,
                                           spatial.col = T, 
                                           spatial.ext = T, sp.model= "IFM",
                                           disp.par= disp.par, 
                                           IFM.bounding.box = 30,
                                           div = div) )
      
      model.list <- list()
      if (class(spatial.model) == "unmarkedFitColExt"){
        model.list[[1]] <- aperm(spatial.model@projected[2,,], c(2,1)) ##matrix sites - seasons
        model.list[[2]] <- spatial.model@estimates
        model.list[[3]] <- spatial.model@AIC
        best_model[[as.character(div)]] <- model.list
      } else{
        best_model[[as.character(div)]] <- list(NA, NA, NA)
      }
      
    }
    
    index_best <- which.min(sapply(best_model, function(x) x[[3]] ) )
    
    dist.models[[dist]] <-  best_model[[index_best]]
    divs[[dist]] <- index_best
    
  }
  
  return(list(dist.models, divs) )
  
}

parallel::stopCluster(cl = my.cluster)

best_div_ifm_habitat <- lapply(fits_IFM_habitat, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[2]])) )
  
fits_IFM_habitat <- lapply(fits_IFM_habitat, function(x) 
  lapply(x, function(y) lapply(y, function(z) z[[1]])) )


names(fits_IFM_habitat) <- dynamics
for (d in dynamics){
  names(fits_IFM_habitat[[d]]) <- scenario_names
}

names(best_div_ifm_habitat) <- dynamics
for (d in dynamics){
  names(best_div_ifm_habitat[[d]]) <- scenario_names
}

save(fits_IFM_habitat,
     file=file.path(output.path, "fits_IFM_habitat.Rdata"))

save(best_div_ifm_habitat,
     file=file.path(output.path, "best_div_ifm_habitat.Rdata"))
