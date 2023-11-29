###############################################################################
##NAME: Two-season IFMs occupancy model
##      with prefixed neighbours occurrence at t0
##AUTHOR: Oriol Solà - o.sola@creaf.uab.cat

##CONTENTS: 
##1. Modifications of the unmarked colext function for the spatial dynocc model
##  as defined by Chandler et al. 2015. Only to be used for a two seasons model
##  with known probability of occurrence at t0 (i.e. SDM for the first season).
##2. Worked example extracted from Solà et al. 2023 (not published yet)

###############################################################################
####  Load required libraries ####
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "Rcpp", "unmarked", "raster"
)
sapply(package_vec, install.load.package)


####  1. Spatial IFM modified colext functions  ####

##Cumulative probability of colonisation or survival (eq. 5) Solà et al. 2023
##Cpp code
src <-
  "
  NumericVector rcpp_y(List disp_kernel_times_psi0, NumericVector base_delta){
    NumericVector gammas (disp_kernel_times_psi0.length());
    for(int i=0; i<disp_kernel_times_psi0.length(); ++i){
      
      double gamma_i = 1;
      NumericVector disp = as<NumericVector>(disp_kernel_times_psi0[i]);
      
      for(int j=0; j<disp.length(); ++j){
        gamma_i *= 1-base_delta[i]*disp[j];
      }
      gammas[i]= 1-gamma_i;
  
    }
  
    return gammas;
  }
  "
Rcpp::cppFunction(src)

##Modified colext function

##New arguments:

##1. Spatial col: whether we use the spatial model for colonisation or not

##2. Spatial ext: whether we use the spatial model for extinction or not

##3. col_disp_kernel_times_psi0 and  ext_disp_kernel_times_psi0: 
##List of vectors with length nº study cells
##Vectors of length nº neighbours study cell i
##For each neighbour j calculates: exp(- distance(i,j) / alpha) * psi(j)
##(The two arguments can be the same list or a different one depending on alpha)

colext2 <- function (psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, 
                     pformula = ~1, data, starts, method = "BFGS", se = TRUE,
                     pre.calcul = T, spatial.col = F, spatial.ext = F,
                     col_disp_kernel_times_psi0, ext_disp_kernel_times_psi0,
                     ...) 
{
  K <- 1
  data@y[data@y > K] <- K
  y <- getY(data)
  J <- numY(data)/data@numPrimary
  M <- nrow(y)
  nY <- ncol(y)/J
  n.det <- sum(apply(y > 0, 1, any, na.rm = TRUE))
  fc <- match.call()
  fc[[1]] <- as.name("colext.fit")
  formula <- list(psiformula = psiformula, gammaformula = gammaformula, 
                  epsilonformula = epsilonformula, pformula = pformula)
  check_no_support(formula)
  fc$formula <- as.name("formula")
  fc$bootstrap.se <- fc$covdata.site <- fc$covdata.obs <- fc$data <- fc$B <- fc$psiformula <- fc$gammaformula <- fc$epsilonformula <- fc$pformula <- NULL
  fc$data <- as.name("data")
  fc$J <- as.name("J")
  fc$method <- as.name("method")
  fc$getHessian <- as.name("se")
  fc$nSites <- as.name("M") ##Track total number of surveyed cells
  fc$pre.calcul <- as.name("pre.calcul")
  fc$spatial.col <- as.name("spatial.col")
  fc$spatial.ext <- as.name("spatial.ext")
  fc$se <- NULL
  if (missing(starts)) {
    fc$starts <- NULL
  }
  else {
    fc$starts <- eval(starts)
  }
  if (missing(col_disp_kernel_times_psi0)) {
    fc$col_disp_kernel_times_psi0 <- NULL
  }
  else {
    fc$col_disp_kernel_times_psi0 <- eval(col_disp_kernel_times_psi0)
  }
  if (missing(ext_disp_kernel_times_psi0)) {
    fc$ext_disp_kernel_times_psi0 <- NULL
  }
  else {
    fc$ext_disp_kernel_times_psi0 <- eval(ext_disp_kernel_times_psi0)
  }
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(fc)))
    for (a in names(extras)[existing]) fc[[a]] <- extras[[a]]
    if (any(!existing)) {
      fc <- as.call(c(as.list(fc), extras[!existing]))
    }
  }
  fm <- eval(fc)
  fm$n.det <- n.det
  opt <- fm$opt
  nP <- fm$nP
  M <- fm$M
  nDP <- fm$nDP
  nGP <- fm$nGP
  nEP <- fm$nEP
  nSP <- fm$nSP
  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  names(ests) <- fm$mle$names
  fmAIC <- 2 * opt$value + 2 * nP
  psiParams <- ests[1:nSP]
  colParams <- ests[(nSP + 1):(nSP + nGP)]
  extParams <- ests[(nSP + nGP + 1):(nSP + nGP + nEP)]
  detParams <- ests[(nSP + nGP + nEP + 1):nP]
  psi <- unmarkedEstimate(name = "Initial", short.name = "psi", 
                          estimates = psiParams, covMat = as.matrix(covMat[1:nSP, 
                                                                           1:nSP]), invlink = "logistic", invlinkGrad = "logistic.grad")
  col <- unmarkedEstimate(name = "Colonization", short.name = "col", 
                          estimates = colParams, covMat = as.matrix(covMat[(nSP + 
                                                                              1):(nSP + nGP), (nSP + 1):(nSP + nGP)]), invlink = "logistic", 
                          invlinkGrad = "logistic.grad")
  ext <- unmarkedEstimate(name = "Extinction", short.name = "ext", 
                          estimates = extParams, covMat = as.matrix(covMat[(nSP + 
                                                                              nGP + 1):(nSP + nGP + nEP), (nSP + nGP + 1):(nSP + 
                                                                                                                             nGP + nEP)]), invlink = "logistic", invlinkGrad = "logistic.grad")
  det <- unmarkedEstimate(name = "Detection", short.name = "p", 
                          estimates = detParams, covMat = as.matrix(covMat[(nSP + 
                                                                              nGP + nEP + 1):nP, (nSP + nGP + nEP + 1):nP]), invlink = "logistic", 
                          invlinkGrad = "logistic.grad")
  estimateList <- unmarkedEstimateList(list(psi = psi, col = col, 
                                            ext = ext, det = det))
  umfit <- new("unmarkedFitColExt", fitType = "colext", 
               call = match.call(), formula = as.formula(paste(unlist(formula), 
                                                               collapse = " ")), psiformula = psiformula, 
               gamformula = gammaformula, epsformula = epsilonformula, 
               detformula = pformula, data = data, sitesRemoved = fm$designMats$removed.sites, 
               estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = opt$value, 
               nllFun = fm$nll, projected = fm$projected, projected.mean = fm$projected.mean, 
               smoothed = fm$smoothed, smoothed.mean = fm$smoothed.mean)
  return(umfit)
}



colext.fit2 <- function (formula, data, J, starts = NULL, method, getHessian = TRUE, 
                         wts, nSites,
                         pre.calcul, spatial.col, spatial.ext,
                         col_disp_kernel_times_psi0= NULL, 
                         ext_disp_kernel_times_psi0= NULL,
                         ...) 
{
  ##Spatial IFM probability of colonisation
  get.gamma <- function(colParams){
    
    ##“baseline” colonization probability
    base_delta <- plogis(X.it.gam %*% colParams)
    
    ##Only study sites kept by the colext function
    sites.index <- setdiff(1:nSites, designMats$removed.sites)
    
    if (nY <= 2){
      
      if (pre.calcul){
        
        gammas <- rcpp_y(col_disp_kernel_times_psi0[sites.index], base_delta)
        
      } else{
        stop("It is necessary to precalculate exp*psi", 
             call. = FALSE)
      }
      
    } else {
      stop("Spatial colonization function not implemented for nY>2 yet", 
           call. = FALSE)
    }
    
    
    return(gammas)
    
  }
  
  ##Spatial IFM probability of extinction
  get.epsilon <- function(extParams){
    
    ##“baseline” colonization probability
    base_delta <- plogis(X.it.eps %*% extParams)
    
    
    sites.index <- setdiff(1:nSites, designMats$removed.sites)
    
    if (nY <= 2){
      
      if (pre.calcul){
        
        epsilons <- 1-rcpp_y(ext_disp_kernel_times_psi0[sites.index], base_delta)
        
      } else{
        
        stop("It is necessary to precalculate exp*psi", 
             call. = FALSE)
      }
      
    } else {
      stop("Spatial colonization function not implemented for nY>2 yet", 
           call. = FALSE)
    }
    
    
    return(epsilons)
    
  }
  
  K <- 1
  designMats <- getDesign(data, formula = as.formula(paste(unlist(formula), 
                                                           collapse = " ")))
  V.itjk <- designMats$V
  X.it.gam <- designMats$X.gam
  X.it.eps <- designMats$X.eps
  W.i <- designMats$W
  detParms <- colnames(V.itjk)
  gamParms <- colnames(X.it.gam)
  epsParms <- colnames(X.it.eps)
  psiParms <- colnames(W.i)
  y <- designMats$y
  M <- nrow(y)
  nY <- ncol(y)/J
  if (missing(wts)) 
    wts <- rep(1, M)
  X.it.gam <- as.matrix(X.it.gam[-seq(nY, M * nY, by = nY), 
  ])
  X.it.eps <- as.matrix(X.it.eps[-seq(nY, M * nY, by = nY), 
  ])
  nDP <- length(detParms)
  nGP <- length(gamParms)
  nEP <- length(epsParms)
  nSP <- length(psiParms)
  nDMP <- 1
  nP <- nDP + nSP + nGP + nEP
  y.itj <- as.numeric(t(y))
  y.itj[is.na(y.itj)] <- 99
  V.itjk[is.na(V.itjk)] <- 9999
  y.it <- matrix(t(y), nY * M, J, byrow = TRUE)
  J.it <- rowSums(!is.na(y.it))
  V.arr <- array(t(V.itjk), c(nDP, nDMP, J, nY, M))
  V.arr <- aperm(V.arr, c(2, 1, 5, 4, 3))
  y.arr <- array(y.itj, c(J, nY, M))
  y.arr <- aperm(y.arr, c(3:1))
  storage.mode(J.it) <- storage.mode(y.arr) <- storage.mode(K) <- "integer"
  alpha <- array(NA, c(K + 1, nY, M))
  forward <- function(detParms, phis, psis, storeAlpha = FALSE) {
    negloglike <- 0
    psiSite <- matrix(c(1 - psis, psis), K + 1, M, byrow = TRUE)
    mp <- array(V.itjk %*% detParms, c(nDMP, J, nY, M))
    for (t in 1:nY) {
      storage.mode(t) <- "integer"
      detVecs <- .Call("getDetVecs", y.arr, mp, J.it[seq(from = t, 
                                                         to = length(J.it) - nY + t, by = nY)], t, K, 
                       PACKAGE = "unmarked")
      psiSite <- psiSite * detVecs
      if (storeAlpha) 
        alpha[, t, ] <<- psiSite[, ]
      if (t < nY) {
        for (i in 1:M) {
          psiSite[, i] <- phis[, , t, i] %*% psiSite[, 
                                                     i]
        }
      } 
      else {
        negloglike <- negloglike - sum(wts * log(colSums(psiSite)))
      }
    }
    negloglike
  }
  backward <- function(detParams, phis) {
    beta <- array(NA, c(K + 1, nY, M))
    for (i in 1:M) {
      backP <- rep(1, K + 1)
      for (t in nY:1) {
        beta[, t, i] <- backP
        detVec <- rep(1, K + 1)
        for (j in 1:J) {
          if (y.arr[i, t, j] != 99) {
            mp <- V.arr[, , i, t, j] %*% detParams
            detVecObs <- .Call("getSingleDetVec", 
                               y.arr[i, t, j], mp, K, PACKAGE = "unmarked")
            detVec <- detVec * detVecObs
          }
        }
        if (t > 1) 
          backP <- t(phis[, , t - 1, i]) %*% (detVec * 
                                                backP)
      }
    }
    return(beta)
  }
  X.gam <- X.it.gam %x% c(-1, 1)
  X.eps <- X.it.eps %x% c(-1, 1)
  phis <- array(NA, c(2, 2, nY - 1, M))
  
  
  nll <- function(params) {
    
    psis <- plogis(W.i %*% params[1:nSP])
    
    colParams <- params[(nSP + 1):(nSP + nGP)]
    extParams <- params[(nSP + nGP + 1):(nSP + nGP + nEP)]
    detParams <- params[(nSP + nGP + nEP + 1):nP]
    if (spatial.col == T){
      gammas <- get.gamma(colParams)
      ##Eliminate boundary probabilities of colonisation (i.e., 0 or 1)
      gammas[gammas==0] <- 10^-16
      gammas[gammas==1] <- 1 - 10^-16
      phis[1, 1, , ] <- 1-gammas
      phis[2, 1, , ] <- gammas
    } else{
      phis[, 1, , ] <- plogis(X.gam %*% colParams)
    }
    if (spatial.ext == T){
      epsilons <- get.epsilon(extParams)
      ##Eliminate boundary probabilities of extinction (i.e., 0 or 1)
      epsilons[epsilons==0] <- 10^-16
      epsilons[epsilons==1] <- 1 - 10^-16
      phis[1, 2, , ] <- epsilons
      phis[2, 2, , ] <- 1-epsilons
    } else {
      phis[, 2, , ] <- plogis(X.eps %*% -extParams)
    }
    forward(detParams, phis, psis) + 0.001 * sqrt(sum(params^2))
  }
  
  
  if (is.null(starts)) 
    starts <- rep(0, nP)
  set.seed(5)
  fm <- optim(starts, nll, method = method, hessian = getHessian, 
              ...)
  mle <- fm$par
  
  psis <- plogis(W.i %*% mle[1:nSP])
  
  colParams <- mle[(nSP + 1):(nSP + nGP)]
  extParams <- mle[(nSP + nGP + 1):(nSP + nGP + nEP)]
  detParams <- mle[(nSP + nGP + nEP + 1):nP]
  
  if (spatial.col == T){
    gammas <- get.gamma(colParams)
    gammas[gammas==0] <- 10^-16
    gammas[gammas==1] <- 1 - 10^-16
    phis[1, 1, , ] <- 1-gammas
    phis[2, 1, , ] <- gammas
  } else{
    phis[, 1, , ] <- plogis(X.gam %*% colParams)
  }
  
  if (spatial.ext == T){
    epsilons <- get.epsilon(extParams)
    epsilons[epsilons==0] <- 10^-16
    epsilons[epsilons==1] <- 1 - 10^-16
    phis[1, 2, , ] <- epsilons
    phis[2, 2, , ] <- 1-epsilons
  } else {
    phis[, 2, , ] <- plogis(X.eps %*% -extParams)
  }
  
  
  projected <- array(NA, c(2, nY, M))
  projected[1, 1, ] <- 1 - psis
  projected[2, 1, ] <- psis
  for (i in 1:M) {
    for (t in 2:nY) {
      projected[, t, i] <- phis[, , t - 1, i] %*% projected[, 
                                                            t - 1, i]
    }
  }
  projected.mean <- apply(projected, 1:2, mean)
  rownames(projected.mean) <- c("unoccupied", "occupied")
  colnames(projected.mean) <- 1:nY
  forward(detParams, phis, psis, storeAlpha = TRUE)
  beta <- backward(detParams, phis)
  gamma <- array(NA, c(K + 1, nY, M))
  for (i in 1:M) {
    for (t in 1:nY) {
      gamma[, t, i] <- alpha[, t, i] * beta[, t, i]/sum(alpha[, 
                                                              t, i] * beta[, t, i])
    }
  }
  smoothed.mean <- apply(gamma, 1:2, mean)
  rownames(smoothed.mean) <- c("unoccupied", "occupied")
  colnames(smoothed.mean) <- 1:nY
  parm.names <- c(psiParms, gamParms, epsParms, detParms)
  mle.df <- data.frame(names = parm.names, value = mle)
  rownames(mle.df) <- paste(c(rep("psi", nSP), rep("col", 
                                                   nGP), rep("ext", nEP), rep("det", nDP)), 
                            c(1:nSP, 1:nGP, 1:nEP, 1:nDP))
  list(mle = mle.df, opt = fm, nP = nP, M = M, nDP = nDP, nGP = nGP, 
       nEP = nEP, nSP = nSP, nllFun = nll, designMats = designMats, 
       projected = projected, projected.mean = projected.mean, 
       smoothed = gamma, smoothed.mean = smoothed.mean)
}

####  2. Tutorial  ####

##0. Load example data

Dir.Base <- getwd()

##Study cells (i.e. surveyed cells)
load(file.path(Dir.Base, "study_cells.Rdata")) ##df: cell_ID|lon|lat

##Landscape cells (i.e. surveyed an unsurveyed cells)
load(file.path(Dir.Base, "landscape_cells.Rdata")) ##df: cell_ID|lon|lat|sdm_cbba2

##Detection data Luscinia megarhynchos
load(file.path(Dir.Base, "detection_data.Rdata"))##df: season-observation (major minor order)



##1. Substitute default colext functions for modified spatial colext functions
environment(colext.fit2) <- asNamespace('unmarked')
assignInNamespace("colext.fit", colext.fit2, ns = 'unmarked')

environment(colext2) <- asNamespace('unmarked')
assignInNamespace("colext", colext2, ns = 'unmarked')



##2. Adjacency lists with neighbour distances and indexes to each study cell

##Index study cells in landscape cells df
index_study_in_landscape <- match(study_cells$cell_id, landscape_cells$cell_id)

##Distances matrix -> rows: study cells, columns: landscape cells
points_distances <- pointDistance(study_cells[,c("lon", "lat")],
                                  landscape_cells[,c("lon", "lat")],
                                  lonlat=F) ## 397 MB
##Distances in km
points_distances <- points_distances/1000

##Eliminate same cell distances
for (i in 1:nrow(points_distances)){
  points_distances[i, index_study_in_landscape[i]] <- NA
}

##Adjacency list with neighbours (at 30km max) indexes
index_cells_30km <- list()
##Adjacency list with neighbours distances
distance_cells_30km <- list()
for (i in 1:nrow(study_cells)){
  neigh <- which(points_distances[i, ] < 30)
  index_cells_30km[[i]] <- neigh
  distance_cells_30km[[i]] <- points_distances[i, neigh]
}



##3. Calculate exp dispersal kernel multiplied by psi of neighbours

alpha <- 2 ##Same alpha for col and ext dispersal kernel, but it could be different
div <- 1 ##Omega in (Eq.3) for likelihood optimisation fine tuning
##List of vectors with length nº study cells
##Vectors of length nº neighbours study cell i
##For each neighbour j calculates: exp(- distance(i,j) / alpha) * psi(j)
disp_kernel_times_psi0 <- list()
for (i in 1:length(distance_cells_30km)){
  product <-
    exp(-distance_cells_30km[[i]]/alpha)*
    landscape_cells$sdm_cbba2[index_cells_30km[[i]]]/div
  disp_kernel_times_psi0[[i]] <- product[product>0.001] ##bounding box to speed up likelihood optimisation
}



##4. Fit IFM spatial model


## Create unmarked multiframe (umf)
SiteCovs <- data.frame(psi0 = 
                         landscape_cells$sdm_cbba2[index_study_in_landscape])
simUMF <- unmarkedMultFrame(
  y = detection_data,
  siteCovs = SiteCovs,
  numPrimary=2)

##Fit model
spatial.model <- unmarked::colext(
  psiformula = ~ psi0, # First-year occupancy
  gammaformula = ~ 1, # Colonization
  epsilonformula = ~ 1, # Extinction
  pformula = ~ 1, # Detection
  data = simUMF,
  spatial.col = T,
  spatial.ext = T,
  col_disp_kernel_times_psi0= disp_kernel_times_psi0,
  ext_disp_kernel_times_psi0 = disp_kernel_times_psi0) 


# fixed.model <- unmarked::colext(
#   psiformula = ~ psi0, # First-year occupancy
#   gammaformula = ~ 1, # Colonization
#   epsilonformula = ~ 1, # Extinction
#   pformula = ~ 1, # Detection
#   data = simUMF,
#   spatial.col = F,
#   spatial.ext = F)



##5. Predict occupancy dynamics study cells

##Probability occupancy t0 (psi0), colonisation (col), extinction (ext)
psi0 <- predict(spatial.model, type="psi")[,1] ##or spatial.model@projected[2,1,]
base_col <- predict(spatial.model, type="col")[,1]
base_ext <- predict(spatial.model, type="ext")[,1]
col <- rcpp_y(disp_kernel_times_psi0, base_col)
ext <- 1-rcpp_y(disp_kernel_times_psi0, base_ext)

##Occupancy time 1
psi1 <- psi0*(1-ext) + (1-psi0)*col
