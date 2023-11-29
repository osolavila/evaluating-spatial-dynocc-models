library(raster)
# library(NLMR)
library(foreach)
library(doParallel)
library(IFMcpp)
library(ggplot2)



output.path <- file.path(Dir.Base, "colext", "Analysis", "Rdata",
                         "revision", "simulation")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}

####Simulate occurrence scenarios

####  1. Landscape design  ####
if (file.exists(file.path(output.path, "landscape.Rdata"))){
  load(file.path(output.path, "landscape.Rdata"))
} else{
  side <- 180
  nsites <- side^2
  size_quadrant <- 100
  num_quadrants <- nsites/size_quadrant
  side_quadrant <- sqrt(nsites/num_quadrants)
  
  
  xcoord <- 1:side
  ycoord <- 1:side
  grid <- as.matrix(expand.grid(x=xcoord, y=ycoord))
  
  
  ##Raster landscape grids
  landscape_raster <- rasterFromXYZ(cbind(grid, 1:nrow(grid)))
  landscape_raster[] <- 1:length(landscape_raster)
  ##Raster quadrants
  quadrants_raster <- aggregate(landscape_raster, side_quadrant)
  quadrants_raster[] <- 1:length(quadrants_raster)
  ##Get quadrants id in landscape grids
  landscape_raster_quadrants_id <- disaggregate(quadrants_raster, side_quadrant)
  
  ##Dataframe with landcape cell info
  landscape_cells <- data.frame(Cell_ID=landscape_raster[],
                                Quadrant_ID=landscape_raster_quadrants_id[])
  landscape_cells <- cbind(landscape_cells, raster::coordinates(landscape_raster))
  
  
  
  ####  Simulate Gaussian field  ####
  rfs <- raster::stack()
  for ( i in 1:10){ 
    
    gaussian_field <- nlm_gaussianfield(ncol = side, nrow = side, autocorr_range = 30,
                                        mag_var = 8, nug = 5, user_seed = i)
    
    rfs <- raster::stack(rfs, gaussian_field)
    # plot(gaussian_field)
  }
  
  plot(rfs) ##visually select variables
  
  gf <- extract(rfs, landscape_cells[,c("x","y")])
  
  psi0_var <- gf[, 4]
  col_var <- gf[, 8]
  ext_var <- gf[, 7]
  
  save(side, nsites, size_quadrant, num_quadrants, side_quadrant,
       landscape_cells, psi0_var, col_var, ext_var, 
       file=file.path(output.path, "landscape.Rdata"))
}

##Autocovariate
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



#### Simulate occurrence  ####
set.seed(10)

nyears <- 2
nscenarios <- 4
ndynamics <- 3
nreps <- 10
zs_array <- array(dim=c(nsites, nyears, nscenarios, ndynamics, nreps))

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

for ( rep in 1:nreps){
  ##Occu t0 
  
  
  ##Probability occupancy t0
  psi0 <- plogis(mean.psi1+beta.psi1*scale(psi0_var))
  
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],psi0)),
  #              col=topo.colors(20), box=FALSE) 
  
  z0 <- rbinom(nsites, 1, psi0)
  zs_array[,1,,,rep] <- z0
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  
  sp_autocovariate <- sapply(index_neighs,
                             function(x) sum(z0[x], na.rm=T)/length(x) ) 
  # hist(sp_autocovariate)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  sp_autocovariate)))
  
  #### Fit habitat  ####
  ##col
  gamma.habitat <- plogis(mean.gamma+beta.gamma*col_var)
  
  
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  gamma.habitat))) 
  
  z.habitat <- rbinom(nsites, 1, z0+(1-z0)*gamma.habitat)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat))) 
  zs_array[,2,1,1,rep] <- z.habitat
  
  ##ext
  epsilon.habitat <- plogis(mean.epsilon+beta.epsilon*ext_var)
  
  
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  epsilon.habitat))) 
  
  z.habitat <- rbinom(nsites, 1, z0*(1-epsilon.habitat))
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat))) 
  zs_array[,2,1,2,rep] <- z.habitat
  ##colext
  z.habitat <- rbinom(nsites, 1, z0*(1-epsilon.habitat)+(1-z0)*gamma.habitat)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat))) 
  zs_array[,2,1,3,rep] <- z.habitat
  
  #### Fit BRM  ####
  
  ##Probability colonization
  gamma.brm <-  plogis(cbind(rep(1, nsites), sp_autocovariate)
                       %*% gamma.params.brm)
  
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  gamma.brm))) 
  # 
  z.brm <- rbinom(nsites, 1, z0+(1-z0)*gamma.brm)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.brm))) 
  # 
  zs_array[,2,2,1,rep] <- z.brm
  
  ##Ext
 
  
  
  epsilon.brm <-  plogis(cbind(rep(1, nsites), sp_autocovariate)
                         %*% epsilon.params.brm)
  
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  epsilon.brm))) 
  
  z.brm <- rbinom(nsites, 1, z0*(1-epsilon.brm) )
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.brm))) 
  
  zs_array[,2,2,2,rep] <- z.brm
  
  ##colext
  z.brm <- rbinom(nsites, 1, z0*(1-epsilon.brm)+(1-z0)*gamma.brm)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.brm))) 
  zs_array[,2,2,3,rep] <- z.brm
  
  #### Fit habitat+BRM  ####
  ##col
  gamma.habitat.brm.1 <- plogis(cbind(rep(1, nsites), col_var, sp_autocovariate)
                                %*% gamma.params.habitat.brm.1)
  gamma.habitat.brm.2 <- plogis(cbind(rep(1, nsites), col_var, sp_autocovariate)
                                %*% gamma.params.habitat.brm.2)
  
  z.habitat.brm.1 <- rbinom(nsites, 1, z0+(1-z0)*gamma.habitat.brm.1)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.1))) 
  
  z.habitat.brm.2 <- rbinom(nsites, 1, z0+(1-z0)*gamma.habitat.brm.2)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.2))) 
  
  zs_array[,2,3,1,rep] <- z.habitat.brm.1
  zs_array[,2,4,1,rep] <- z.habitat.brm.2
  
  ##ext
  epsilon.habitat.brm.1 <- plogis(cbind(rep(1, nsites), ext_var, sp_autocovariate)
                                  %*% epsilon.params.habitat.brm.1)
  epsilon.habitat.brm.2 <- plogis(cbind(rep(1, nsites), ext_var, sp_autocovariate)
                                  %*% epsilon.params.habitat.brm.2)
  
  z.habitat.brm.1 <- rbinom(nsites, 1, z0*(1-epsilon.habitat.brm.1) )
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.1))) 
  
  z.habitat.brm.2 <- rbinom(nsites, 1,  z0*(1-epsilon.habitat.brm.1))
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.2))) 
  
  zs_array[,2,3,2,rep] <- z.habitat.brm.1
  zs_array[,2,4,2,rep] <- z.habitat.brm.2
  
  
  ##colext
  z.habitat.brm.1 <- rbinom(nsites, 1, z0*(1-epsilon.habitat.brm.1)+
                              (1-z0)*gamma.habitat.brm.1)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.1))) 
  
  z.habitat.brm.2 <- rbinom(nsites, 1, z0*(1-epsilon.habitat.brm.2)+
                              (1-z0)*gamma.habitat.brm.2)
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z0))) 
  # raster::plot(rasterFromXYZ(cbind(landscape_cells[,c("x","y")],
  #                                  z.habitat.brm.2))) 
  
  zs_array[,2,3,3,rep] <- z.habitat.brm.1
  zs_array[,2,4,3,rep] <- z.habitat.brm.2


}

save(mean.psi1, beta.psi1,
     mean.gamma, beta.gamma,
     gamma.params.brm,
     gamma.params.habitat.brm.1,
     gamma.params.habitat.brm.2,
     nyears, nscenarios, ndynamics, nreps, 
     zs_array,
     file=file.path(output.path, "occu_sim.Rdata"))

#### Plot probability of colonisation ####

xyz.df <- expand.grid(conn=seq(0,1,0.1),hab=seq(0,1,0.1))

xyz.df$habitat <- plogis(mean.gamma + xyz.df$hab*beta.gamma) 
xyz.df$brm <- plogis(gamma.params.brm[1] + xyz.df$conn*gamma.params.brm[2])
xyz.df$habitat_brm_1 <- plogis(gamma.params.habitat.brm.1[1] + 
                                 xyz.df$hab*gamma.params.habitat.brm.1[2] +
                                 xyz.df$conn*gamma.params.habitat.brm.1[3])
xyz.df$habitat_brm_2 <- plogis(gamma.params.habitat.brm.2[1] + 
                                 xyz.df$hab*gamma.params.habitat.brm.2[2] +
                                 xyz.df$conn*gamma.params.habitat.brm.2[3])

# plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +  # Adjust colors as needed
#   labs(x = "connectivity", y = "colonization habitat covariable", fill = "prob col",
#        title = "habitat") +
#   theme_minimal() +
#   theme( plot.title = element_text(hjust = 0.5) )
# 
# # Display the plot
# print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +  # Adjust colors as needed
  labs(x = "connectivity", y = "colonization covariable", fill = "prob col") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = brm)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +  # Adjust colors as needed
  labs(x = "connectivity", y = "colonization covariable", fill = "prob col") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat_brm_1)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +  # Adjust colors as needed
  labs(x = "connectivity", y = "colonization covariable", fill = "prob col") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat_brm_2)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "green", limits = c(0, 1)) +  # Adjust colors as needed
  labs(x = "connectivity", y = "colonization covariable", fill = "prob col") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)


#### Plot probability of extinction ####

xyz.df <- expand.grid(conn=seq(0,1,0.1),hab=seq(0,1,0.1))

xyz.df$habitat <- plogis(mean.epsilon + xyz.df$hab*beta.epsilon) 
xyz.df$brm <- plogis(epsilon.params.brm[1] + xyz.df$conn*epsilon.params.brm[2])
xyz.df$habitat_brm_1 <- plogis(epsilon.params.habitat.brm.1[1] + 
                                 xyz.df$hab*epsilon.params.habitat.brm.1[2] +
                                 xyz.df$conn*epsilon.params.habitat.brm.1[3])
xyz.df$habitat_brm_2 <- plogis(epsilon.params.habitat.brm.2[1] + 
                                 xyz.df$hab*epsilon.params.habitat.brm.2[2] +
                                 xyz.df$conn*epsilon.params.habitat.brm.2[3])

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +  # Adjust extors as needed
  labs(x = "connectivity", y = "extinction covariable", fill = "prob ext") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = brm)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +  # Adjust extors as needed
  labs(x = "connectivity", y = "extinction covariable", fill = "prob ext") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat_brm_1)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +  # Adjust extors as needed
  labs(x = "connectivity", y = "extinction covariable", fill = "prob ext") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )

# Display the plot
print(plot)

plot <- ggplot(xyz.df, aes(x = conn, y = hab, fill = habitat_brm_2)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +  # Adjust extors as needed
  labs(x = "connectivity", y = "extinction covariable", fill = "prob ext") +
  theme_minimal() + guides(fill = "none") +
  theme(
    axis.title.x = element_text(size = 26),  # Adjust the size as needed
    axis.title.y = element_text(size = 26),   # Adjust the size as needed
    axis.text.x = element_text(size = 25),  # Adjust the size as needed
    axis.text.y = element_text(size = 25)   # Adjust the size as needed
  )


# Display the plot
print(plot)

#### Legends prob col and prob ext ####
# Create a separate horizontal legend without the main plot
legend_plot <- ggplot() +
  geom_tile(aes(x = 1, y = 1, fill = c(0, 1)), width = 1, height = 1) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_void() +
  theme(
    legend.position = "bottom",  # Place the legend at the bottom
    legend.direction = "horizontal"  # Make the legend horizontal
  ) +
  guides(fill = guide_colorbar(title = "prob ext"))

# Display the legend plot
print(legend_plot)

# Create a separate horizontal legend without the main plot
legend_plot <- ggplot() +
  geom_tile(aes(x = 1, y = 1, fill = c(0, 1)), width = 1, height = 1) +
  scale_fill_gradient(low = "white", high = "green") +
  theme_void() +
  theme(
    legend.position = "bottom",  # Place the legend at the bottom
    legend.direction = "horizontal"  # Make the legend horizontal
  ) +
  guides(fill = guide_colorbar(title = "prob col"))

# Display the legend plot
print(legend_plot)


####  Sampling scenarios (spatial coverage)  ####
set.seed(14)
prop.surveyed <- 0.05

##array: sites,prop.surveyed.scenarios, replicates
surveyed.cells.array <- array(NA, 
                              dim=c(nsites, nreps) )

for (replicate in 1:nreps){
  
  # sample_cells <- function(prop.surveyed)
  surveyed_cells <- unlist(lapply(1:num_quadrants, function(x)
    landscape_cells$Cell_ID[landscape_cells$Quadrant_ID == x][
      sample(1:size_quadrant, floor(size_quadrant*prop.surveyed))]
  ))
  surveyed_cells <- 1:nsites %in% surveyed_cells
  
  surveyed.cells.array[,replicate] <- surveyed_cells
  
  
}

save(surveyed.cells.array, prop.surveyed,
     file=file.path(output.path, paste0("surveyed.cells.array.Rdata") ) )

#### Detection scenarios  ####

nsurveys <- 2
nsurveyed <- nsites*prop.surveyed
p.det <- 0.5


y.array <- array(NA, dim=c(nreps,nsurveyed, nsurveys, nyears,
                           nscenarios,ndynamics))

set.seed(29)

##
for ( rep in 1:nreps){
  z <- zs_array[surveyed.cells.array[,rep],1,1,1,rep]
  for (s in 1:nsurveys){
    y.array[rep,,s,1,,] <- rbinom(length(z), 1, z*p.det)
  }
}


for ( rep in 1:nreps){
  for (i in 1:ndynamics){ ## col, ext or colext
    for (scn in 1:nscenarios){
      z <- zs_array[surveyed.cells.array[,rep],2,scn,i,rep]
      for (s in 1:nsurveys){
        y.array[rep,,s,2,scn,i] <- rbinom(length(z), 1, z*p.det)
      }
      
    }
  }
}
save(y.array, nsurveys, nsurveyed, p.det,
     file=file.path(output.path, paste0("y.array.Rdata") ) )

