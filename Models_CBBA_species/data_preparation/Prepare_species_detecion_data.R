rm(list = ls())

##Load necessary libraries
library("readxl")
library(tidyverse)
library(sf)
library(lubridate) #for ydate
library(foreign) #read dbf

################################################################################
##Load CBBA2 and CBBA3 data
################################################################################

##AONC (CBBA) database
db_aonc_survey_info <- read_excel("./Data/AONC/Dades_1x1.xlsx", sheet = "1x1_header") #Does not read Init_time
db_aonc_species_detection <- read_excel("./Data/AONC/Dades_1x1.xlsx", sheet = "1x1_detail")
db_aonc_change <- read_excel("./Data/AONC/Dades_1x1.xlsx", sheet = "Canvi_AONC3_final")

aonc3_metadata <- read_excel("./Data/AONC/aonc3_species.xlsx")

aonc2_sdm_db <- read.csv2("./Data/AONC/aonc2.csv")
aonc2_metadata <- read_excel("./Data/AONC/aonc2_metadata.xlsx")
aonc2_metadata$Scientific_name[
  aonc2_metadata$Scientific_name=="Regulus ignicapillus"] <- "Regulus ignicapilla"

##Color palettes
aonc3_paleta_canvi <- read.dbf("./Data/AONC/aonc3_paleta_canvi.dbf")
colors_canvi <- rgb(aonc3_paleta_canvi[1:21,c(2,3,4)], maxColorValue = 255)

aonc3_paleta_model <- read.dbf("./Data/AONC/aonc3_paleta_model.dbf")
colors_sdm <- rgb(aonc3_paleta_model[1:20,c(1,2,3)], maxColorValue = 255)


################################################################################
##Data cleaning
################################################################################

## 1 Find ID survey duplicates and eliminate them
n_occur <- data.frame(table(db_aonc_survey_info$id_1x1_header))
##n_occur[n_occur$Freq > 1,] #Found 1 ID duplicate for 2 different cells
junk_cells_id <- db_aonc_survey_info$cell_1x1[db_aonc_survey_info$id_1x1_header
                                           %in% n_occur$Var1[n_occur$Freq > 1]]

db_aonc_survey_info <- db_aonc_survey_info[
  !(db_aonc_survey_info$id_1x1_header %in% n_occur$Var1[n_occur$Freq > 1]),]
db_aonc_species_detection <- db_aonc_species_detection[
  !(db_aonc_species_detection$id_1x1_header %in% n_occur$Var1[n_occur$Freq > 1]),]

db_aonc_change <- db_aonc_change[!(db_aonc_change$cell_1x1 %in% junk_cells_id),]

rm(n_occur, junk_cells_id)
## 2 Same variables same name
#Species name variable has different name for detection and change db
names(db_aonc_species_detection)[names(db_aonc_species_detection)=="sci_name"] <- "spp"

## 3 Find repeated registers of a species detection in a survey.id and only leave 1
db_aonc_species_detection <- db_aonc_species_detection[
  !duplicated(db_aonc_species_detection[,c("id_1x1_header","spp")]),]

## 4 Find surveys without date and delete them
no_date_surveys <- db_aonc_survey_info[is.na(db_aonc_survey_info$survey_date),
                                       c("id_1x1_header","cell_1x1") ]

db_aonc_survey_info <- db_aonc_survey_info[
  !(db_aonc_survey_info$id_1x1_header %in% no_date_surveys$id_1x1_header),]
db_aonc_species_detection <- db_aonc_species_detection[
  !(db_aonc_species_detection$id_1x1_header %in% no_date_surveys$id_1x1_header),]

##It could be that a cell has presence (in change db) but no surveys detect it, 
##but should not change analysis
##Should check it!
rm(no_date_surveys)

################################################################################
##Join data
################################################################################
survey_info <- c("id_1x1_header","cell_1x1", "survey_number",
                 "Survey_timing", "survey_date" )
change_info <- c("cell_1x1", "spp", "aonc2", "aonc3", "latitude", "longitude")
detection_info <- c("id_1x1_header", "spp", "species_counts")     

aonc_complete <- db_aonc_survey_info[,survey_info] %>%
  left_join(db_aonc_change[,change_info], by= "cell_1x1" ) %>%
  left_join(db_aonc_species_detection[,detection_info],
            by=c("id_1x1_header","spp"))  %>%
  mutate(species_counts=ifelse(is.na(species_counts),0,species_counts),
         year=as.numeric(substr(survey_date,1,4)),
         year_day=yday(survey_date)) 

rm(survey_info, change_info, detection_info)

################################################################################
##Select study cells (cells surveyed at least 2 times each atlas)
################################################################################

cells_info <- db_aonc_change %>% dplyr::select(cell_1x1, longitude, latitude) %>%
  distinct()
##Change coordinates reference system
cell_coord<-st_as_sf(cells_info[,c("longitude","latitude")],
                     coords = c("longitude","latitude"),
                     crs = st_crs(4326) ) ## Old: World Geodetic System 1984, used in GPS
cell_coord <-st_transform(cell_coord, crs = 25831) ##New: ETRS89 / UTM zone 31N
cells_info[,c("longitude","latitude")] <- st_coordinates(cell_coord)



##  Add year and yday to survey info db
db_aonc_survey_info <- db_aonc_survey_info %>%
  mutate(year = as.numeric(substr(survey_date,1,4)),
         year_day=yday(survey_date)) 

#### Get number of night and day surveys for each cell and aonc ####
cells_sampling_periods <- db_aonc_survey_info %>% group_by(cell_1x1) %>%
  summarise(num_d_2= length(which(year<2010 & Survey_timing == "d")),
            num_n_2= length(which(year<2010 & Survey_timing == "n")),
            num_d_3= length(which(year>2010 & Survey_timing == "d")),
            num_n_3= length(which(year>2010 & Survey_timing == "n")))

#### Define cells used for colext calibration ####

##Cells ID (cell_1x1) with 2 sampling periods in both aonc2 and aonc3
cells_study <- cells_sampling_periods %>%
  filter(num_d_2>=2, num_d_3>=2) %>% dplyr::select(cell_1x1) %>%
  filter(cell_1x1 %in% cells_info$cell_1x1)

cells_study <- left_join(cells_study, cells_info, by="cell_1x1")


##Only keep cells_study
rm(cells_sampling_periods, db_aonc_survey_info, cells_info, cell_coord)

################################################################################
##Landscape grid cells
################################################################################
##Grid cells_1x1  UTM coordinates
cells_1x1 <- read_sf("./Data/AONC/grid_1x1", layer = "grid_1x1")
cells_id_coord <- data.frame(cbind(cells_1x1$CELL, 
                                   st_coordinates(st_centroid(cells_1x1)) )) #ETRS89 / UTM zone 31N
names(cells_id_coord) <- c("cell_1x1", "lon", "lat")

##Grid cells_10x10  UTM coordinates
grid_10x10 <- read_sf("./Data/AONC/quadricules-utm-v1r0-2021",
                      layer = "quadricules-utm-10km-v1r0-2021") 
cells_10x10_id_coord <- data.frame(cbind(grid_10x10$COORD_10K, 
                                         st_coordinates(st_centroid(grid_10x10)) )) #ETRS89 / UTM zone 31N
names(cells_10x10_id_coord) <- c("cell_10x10", "lon", "lat")
cells_10x10_id_coord$cell_10x10 <- sapply(cells_10x10_id_coord$cell_10x10, 
                                          function(x)substring(x,4)) 


################################################################################
##Species info
################################################################################

##DeC?ceres_2013 openland and forest species
species.groups <- read_excel("./Data/AONC/Bird_species_groups.xlsx")
species.groups <- data.frame(spp=as.vector(as.matrix(species.groups)),
                             group=rep(names(species.groups), each=nrow(species.groups) ))
species.groups <- species.groups[!is.na(species.groups$spp),]

##Spp abundance
spp_abundance <- db_aonc_change %>% group_by(spp) %>%
  summarise(n_occu_aonc2 = sum(aonc2, na.rm=T), 
            n_occu_aonc3= sum(aonc3, na.rm=T)/2)

species.groups <- left_join(species.groups, spp_abundance, by="spp" )
species.groups <- species.groups[!is.na(species.groups$n_occu_aonc2),]

rm(spp_abundance)