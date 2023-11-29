library(unmarked)
coordinates <- sp::coordinates


get_UMF <- function(study_species, 
                    LiDAR=F,
                    CORINE=F,
                    habitat_aonc3=F,
                    SDM=F,
                    topo=F,
                    climate=F){
  
  ####  Get survwys information of study species  ####
  
  ##Select AONC day surveys information of species of interest in study cells   
  
  spp_detection <- aonc_complete %>% filter(spp==study_species,
                                            cell_1x1 %in% cells_study$cell_1x1,
                                            Survey_timing=="d")
  
  
  #### Get detection data to calibrate colext ####
  
  ##For each year and sampling period obtain df: cell_1x1|presence
  aonc2_s1 <- spp_detection %>% filter(year<2010, survey_number==1) %>%
    dplyr::select(cell_1x1, species_counts)
  aonc2_s2 <- spp_detection %>% filter(year<2010, survey_number==2) %>%
    dplyr::select(cell_1x1, species_counts)
  aonc3_s1 <- spp_detection %>% filter(year>2010, survey_number==1) %>%
    dplyr::select(cell_1x1, species_counts)
  aonc3_s2 <- spp_detection %>% filter(year>2010, survey_number==2) %>%
    dplyr::select(cell_1x1, species_counts)
  
  ##  df with all study cells: cell_1x1| aonc2_s1 | aonc2_s2 | aonc3_s1 | aonc3_s2
  y <- cells_study[,"cell_1x1", drop=F] %>% 
    left_join(aonc2_s1, by="cell_1x1") %>%
    left_join(aonc2_s2, by="cell_1x1") %>%
    left_join(aonc3_s1, by="cell_1x1") %>%
    left_join(aonc3_s2, by="cell_1x1")
  
  y[is.na(y)] <- 0
  
  ##Transform into data ready for colext function
  yy <- as.matrix(y[,2:5])
  colnames(yy)<- c("aonc2_s1" , "aonc2_s2" , "aonc3_s1" , "aonc3_s2")
  

  #### Get obsCovs (site-year-observation covariates) ####
  
  ##For each year and sampling period obtain df: cell_1x1|julian
  aonc2_s1 <- db_aonc_survey_info %>% 
    filter(year<2010, survey_number==1, Survey_timing=="d") %>%
    dplyr::select(cell_1x1, year_day)
  aonc2_s2 <- db_aonc_survey_info %>% 
    filter(year<2010, survey_number==2, Survey_timing=="d") %>%
    dplyr::select(cell_1x1, year_day)
  aonc3_s1 <- db_aonc_survey_info %>% 
    filter(year>2010, survey_number==1, Survey_timing=="d") %>%
    dplyr::select(cell_1x1, year_day)
  aonc3_s2 <- db_aonc_survey_info %>% 
    filter(year>2010, survey_number==2, Survey_timing=="d") %>%
    dplyr::select(cell_1x1, year_day)
  
  julian_days <- cells_study[,"cell_1x1"] %>% 
    left_join(aonc2_s1, by="cell_1x1") %>%
    left_join(aonc2_s2, by="cell_1x1") %>%
    left_join(aonc3_s1, by="cell_1x1") %>%
    left_join(aonc3_s2, by="cell_1x1")
  
  
  ##Transform into data ready for colext function
  julian_days <- as.matrix(julian_days[,2:5])
  colnames(julian_days)<- c("aonc2_s1" , "aonc2_s2" , "aonc3_s1" , "aonc3_s2")
  
  julian_days[,] <- scale(c(julian_days))
  
  obsCovs <- data.frame(julian_day_std=as.vector(t(julian_days)))
  
  
  
  #### Get SiteCovs (site covariates) ####
  SiteCovs <- cells_study[,"cell_1x1", drop=F]
  
  
  if (LiDAR){
    
    vars_lidar <- c("perc.forest",              "ab",                      
                    "cc",                       "den",                     
                    "dbhm",                     "hmitjana",                
                    "perc.closed.forest")
    
    if (file.exists("./Data/Site_covariables/site.covs.Rdata")
        ){
      load("./Data/Site_covariables/site.covs.Rdata")
    }
    else{
      ##Build site.covs dataframe
      ##LIDAR variables
      mean.lidar.1km.2016 <- raster::stack("./Data/Site_covariables/Lidar_variables_biofisiques/mean.forest.structure.vars.1km.2016.tif")
      perc.closed.1km.2016 <- raster("./Data/Site_covariables/Lidar_variables_biofisiques/perc.closed.1km.2016.tif")
      perc.forest.1km.2016 <- raster("./Data/Site_covariables/Lidar_variables_biofisiques/perc.forest.1km.2016.tif")
      
      cells.mean.lidar.1km.2016 <- raster::extract(mean.lidar.1km.2016, cells_study[,2:3])
      cells.mean.lidar.1km.2016[is.na(cells.mean.lidar.1km.2016)]<-0
      ##Percentage of forest with cc>100 within forest cells
      cells.perc.closed.1km.2016 <- raster::extract(perc.closed.1km.2016, cells_study[,2:3])
      cells.perc.closed.1km.2016[is.na(cells.perc.closed.1km.2016)]<-0
      ##Percentage of forest cells
      cells.perc.forest.1km.2016 <- raster::extract(perc.forest.1km.2016, cells_study[,2:3])
      
      cells_yearly_site_covs <- cells.mean.lidar.1km.2016
      cells_yearly_site_covs$perc.forest <- cells.perc.forest.1km.2016
      cells_yearly_site_covs$perc.closed.forest <- cells.perc.closed.1km.2016
      
      site.covs <- cells_yearly_site_covs
      
      save(site.covs, file="./Data/Site_covariables/site.covs.Rdata")
    }
    
    SiteCovs <- cbind(SiteCovs, site.covs[, vars_lidar])
  }
  ##AONC2 SDM
  if (SDM){
    sp.code <- aonc2_metadata$Species_code[aonc2_metadata$Scientific_name==study_species]
    
    aonc2_sdm_db_sp <- aonc2_sdm_db %>% filter(specie == sp.code) %>% 
      dplyr::select(cell_1x1 = cell, predicted_value, predicted_zone)
    gc()
    
    aonc2_sdm <- cells_study %>% 
      left_join(aonc2_sdm_db_sp, by="cell_1x1") %>% 
      mutate(aonc2_sdm = replace_na(as.numeric(predicted_value), 0), 
             aonc2_sdm_cat= as.factor(
               replace_na(as.numeric(predicted_zone) ,0)) ) %>%
      dplyr::select(aonc2_sdm, aonc2_sdm_cat) 
    
    SiteCovs <- cbind(SiteCovs, aonc2_sdm)
  
  }
  if (habitat_aonc3){
    load("./Data/Site_covariables/habitat_vars_aonc3.Rdata")
    SiteCovs <- cbind(SiteCovs, habitat_vars_aonc3)
  }
  if (topo) {
    load("./Data/Site_covariables/topo_vars_aonc3.Rdata")
    SiteCovs <- cbind(SiteCovs, topo_vars_aonc3)
  }
  if (climate) {
    load("./Data/Site_covariables/climate_vars_aonc3.Rdata")
    SiteCovs <- cbind(SiteCovs, climate_vars_aonc3)
  }
  

  #### Get yearlySiteCovs (site-year covariates) ####
  
  atlas <- data.frame(matrix(rep(2:3, each=nrow(SiteCovs)), nrow(SiteCovs), 2))
  atlas <- data.frame(lapply(atlas, as.factor))
  ## We don't need for a model with 2 seasons since colext only take t-1 vars
  # colext_cells_yearly_site_covs <- rbind(cells_yearly_site_covs, cells_yearly_site_covs)
  # colext_cells_yearly_site_covs[seq(1,nrow(colext_cells_yearly_site_covs),2),] <- cells_yearly_site_covs
  # colext_cells_yearly_site_covs[seq(2,nrow(colext_cells_yearly_site_covs),2),] <- NA
  # 
  
  ####  Handle missing covariates  ####
  
  sites.to.remove <- apply(SiteCovs, 1, function (x)
    any(is.na(x)))
  
  yy <- yy[!sites.to.remove, , drop = FALSE]
  
  SiteCovs <- SiteCovs[!sites.to.remove, , drop = FALSE]
  atlas <- atlas[!sites.to.remove, , drop = FALSE]
  obsCovs <- obsCovs[!sites.to.remove[rep(1:nrow(yy), each = 4)], , drop = FALSE]
  
  
  simUMF <- unmarkedMultFrame(
    y = yy,
    siteCovs = SiteCovs,
    obsCovs = obsCovs,
    yearlySiteCovs=list(atlas=atlas),
    numPrimary=2)
}

