get_cat_site_covs <- function(){
  
  ##CORINE variables
  clc2000_1km <- raster::stack("./Data/Site_covariables/CORINE/clc2000_1km.tif")
  clc2018_1km <- raster::stack("./Data/Site_covariables/CORINE/clc2018_1km.tif")
  
  
  cat_site_covs_2000 <- raster::extract(clc2000_1km, cells_id_coord[,2:3])
  
  cells.mosaic.2000 <- get.mosaic(cat_site_covs_2000)
  
  cat_site_covs_2018 <- raster::extract(clc2018_1km, cells_id_coord[,2:3])
  
  
  ##Diff % LC between 2018 and 2000 of all study cells
  cat_LC_change <- as.data.frame(cat_site_covs_2018 - cat_site_covs_2000)
  names(cat_LC_change) <- paste0("change_", names(cat_LC_change))
  
  
  cat_lcc_type <- get.lcc.type(cat_LC_change)
  cat_lcc_type[cat_lcc_type=="broad.to.mixed"] <- "no.change"
  cat_lcc_intensity_cat <- get.change.intensity.category(cat_LC_change)
  cat_lcc_intensity <- get.change.intensity(cat_LC_change)
  
  cat_site_covs <- cat_LC_change
  
  cat_site_covs$lcc_type <- as.factor(cat_lcc_type)
  cat_site_covs$lcc_intensity <- cat_lcc_intensity
  cat_site_covs$lcc_intensity_cat <- as.factor(cat_lcc_intensity_cat)
  cat_site_covs$cell_1x1 <- cells_id_coord$cell_1x1
  
  cat_site_covs <- cbind(cat_site_covs, cat_site_covs_2000)
  cat_site_covs$lc_mosaic_2000 <- cells.mosaic.2000
  
  ##LIDAR variables
  mean.lidar.1km.2016 <- raster::stack("./Data/Site_covariables/Lidar_variables_biofisiques/mean.forest.structure.vars.1km.2016.tif")
  perc.closed.1km.2016 <- raster("./Data/Site_covariables/Lidar_variables_biofisiques/perc.closed.1km.2016.tif")
  perc.forest.1km.2016 <- raster("./Data/Site_covariables/Lidar_variables_biofisiques/perc.forest.1km.2016.tif")
  
  cells.mean.lidar.1km.2016 <- raster::extract(mean.lidar.1km.2016, cells_id_coord[,2:3])
  cells.mean.lidar.1km.2016[is.na(cells.mean.lidar.1km.2016)]<-0
  ##Percentage of forest with cc>100 within forest cells
  cells.perc.closed.1km.2016 <- raster::extract(perc.closed.1km.2016, cells_id_coord[,2:3])
  cells.perc.closed.1km.2016[is.na(cells.perc.closed.1km.2016)]<-0
  ##Percentage of forest cells
  cells.perc.forest.1km.2016 <- raster::extract(perc.forest.1km.2016, cells_id_coord[,2:3])
  
  cat_site_covs$perc.forest <- cells.perc.forest.1km.2016
  cat_site_covs <- cbind(cat_site_covs, cells.mean.lidar.1km.2016)
  cat_site_covs$perc.closed.forest <- cells.perc.closed.1km.2016
  
  save(cat_site_covs, file="./colext/cat_site_covs.Rdata")
  
  ##Add AONC3 site covs
  
  load("./colext/cat_site_covs.Rdata")
  
  # CHC <- read_sf(dsn = "./Data/Site_covariables/CHC", layer = "CHIC_v2_abril2019")
  
  
  habitat_files_dir <- "./Data/Site_covariables/predictors_atles3"
  habitat_var_aonc3 <- list.files(habitat_files_dir)[
    grepl("asc",
          list.files(habitat_files_dir))]
  
  
  
  habitat_variables <- cells_id_coord[,"cell_1x1", drop=F]
  
  
  for (var in habitat_var_aonc3){
    print(var)
    hab_raster <- raster(paste(habitat_files_dir, var, sep = "/"))
    hab_raster <- focal(hab_raster, w=matrix(1,5,5), 
                        fun=mean, na.rm = T, NAonly = T) ##Put neighbours mean values to NA cells
    
    var_name <- substring(var, 1, nchar(var)-4)
    
    habitat_variables[,var_name] <- extract(hab_raster, 
                                            as.matrix(cells_id_coord[,2:3]) )
    
  }
  
  
  cat_site_covs <- cbind(cat_site_covs, habitat_variables)
  
  save(cat_site_covs,
       file="./Data/Site_covariables/cat_site_covs.Rdata")

  return(cat_site_covs)
}



get_cat_sdm <- function(study_species){
  ##SDM AONC2
  ##Cells >10km from detected presence are NA -> 0
  
  
  sp.code <- aonc2_metadata$Species_code[aonc2_metadata$Scientific_name==study_species]
  
  
  
  cat_aonc2_sdm <- aonc2_sdm_db %>% filter(specie == sp.code) %>% 
    dplyr::select(cell_1x1 = cell, aonc2_sdm= predicted_value,
                  aonc2_sdm_cat = predicted_zone) %>%
    mutate(aonc2_sdm = as.numeric(aonc2_sdm), 
           aonc2_sdm_cat= as.factor(
             replace_na(as.numeric(aonc2_sdm_cat) ,0)) ) 
  
  
  cat_aonc2_sdm <- cells_id_coord[,"cell_1x1", drop = F] %>% left_join(cat_aonc2_sdm, by="cell_1x1")
  
  cat_aonc2_sdm$aonc2_sdm[is.na(cat_aonc2_sdm$aonc2_sdm)]<-0
  cat_aonc2_sdm$aonc2_sdm_cat[is.na(cat_aonc2_sdm$aonc2_sdm_cat)]<-0
  
  return(cat_aonc2_sdm)
  
  

}

##Not useful for predictions since it is scaled with the whole dataset instead
##Of with the study cells only
get_cat_neigh_sdm <- function(study_species, max_km){
  ##Sum neighbours
  max_dist <- max_km
  sdm_neigh <- vector(length=nrow(cells_id_coord))
  
  cat_aonc2 <- get_cat_sdm(study_species)
  cat_aonc_cat <- as.integer(as.character(cat_aonc2$aonc2_sdm_cat))
  
  ##Eliminate same cell distances
  for (i in 1:nrow(cells_id_coord)){
    points_distances <- pointDistance(cells_id_coord[i,c("lon", "lat")],
                                      cells_id_coord[,c("lon", "lat")],
                                      lonlat=F)
    points_distances <- points_distances/1000
    points_distances[i] <- NA
    
    max_dist_index <- which(points_distances<max_dist)
    sdm_neigh[i]<- sum(cat_aonc_cat[max_dist_index], na.rm=T)/length(max_dist_index)
    
  }

  sdm_neigh <- scale(sdm_neigh)[,1]
  
  return(sdm_neigh)
  
}

# cat_site_covs <- cat_site_covs %>% left_join(cat_aonc2_sdm, by="cell_1x1")