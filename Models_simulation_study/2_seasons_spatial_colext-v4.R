####  Load required libraries ####
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "Rcpp", "unmarked", "raster", "tictoc"
)
sapply(package_vec, install.load.package)


##Cumulative probability of colonisation or survival (eq. 5) Solà et al. 2023
##Cpp code
# src <-
#   "
#   NumericVector rcpp_y(List disp_kernel_times_psi0, NumericVector base_delta){
#     NumericVector gammas (disp_kernel_times_psi0.length());
#     for(int i=0; i<disp_kernel_times_psi0.length(); ++i){
#       
#       double gamma_i = 1;
#       NumericVector disp = as<NumericVector>(disp_kernel_times_psi0[i]);
#       
#       for(int j=0; j<disp.length(); ++j){
#         gamma_i *= 1-base_delta[i]*disp[j];
#       }
#       gammas[i]= 1-gamma_i;
#   
#     }
#   
#     return gammas;
#   }
#   "
# Rcpp::cppFunction(src)
# 
# 
# src_2 <-
#   "
#   NumericVector rcpp_phi(
#                         List distance_neighs, 
#                         List index_neighs,
#                         NumericVector base_delta,
#                         NumericVector psis_landscape,
#                         double alpha,
#                         double div){
#     NumericVector gammas (distance_neighs.length());
#     for(int i=0; i<distance_neighs.length(); ++i){
#       
#       double gamma_i = 1;
#       NumericVector dist = as<NumericVector>(distance_neighs[i]);
#       NumericVector index = as<NumericVector>(index_neighs[i]);
#       
#       for(int j=0; j<dist.length(); ++j){
#         gamma_i *= 1 - base_delta[i] * std::exp(-dist[j]/alpha) *
#         psis_landscape[ index[j] - 1 ] / div;
#       }
#       gammas[i]= 1-gamma_i;
#   
#     }
#   
#     return gammas;
#   }
#   "
# 
# Rcpp::cppFunction(src_2)

##Modified colext function

##New arguments:

##1. Spatial col: whether we use the spatial model for colonisation or not

##2. Spatial ext: whether we use the spatial model for extinction or not

##3. ...

spatial.colext <- function (psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, 
                     pformula = ~1, data, starts, method = "BFGS", se = TRUE,
                     spatial.col = F, spatial.ext = F, sp.model, 
                     disp.par, div=1, IFM.bounding.box,
                     fixed.landscape.psi0 = F, psi0.var.name,
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
  fc[[1]] <- as.name("spatial.colext.fit")
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
  fc$fixed.landscape.psi0 <- as.name("fixed.landscape.psi0")
  fc$spatial.col <- as.name("spatial.col")
  fc$spatial.ext <- as.name("spatial.ext")
  fc$sp.model <- as.name("sp.model")
  fc$disp.par <- as.name("disp.par")
  fc$div <- as.name("div")
  fc$se <- NULL
  if (missing(starts)) {
    fc$starts <- NULL
  }
  else {
    fc$starts <- eval(starts)
  }
  if (missing(psi0.var.name)) {
    fc$psi0.var.name <- NULL
  }
  else {
    fc$psi0.var.name <- eval(psi0.var.name)
  }
  if (missing(IFM.bounding.box)) {
    fc$IFM.bounding.box <- NULL
  }
  else {
    fc$IFM.bounding.box <- eval(IFM.bounding.box)
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



spatial.colext.fit <- function (formula, data, J, starts = NULL, method, getHessian = TRUE, 
                         wts, nSites,
                         spatial.col = F, spatial.ext = F, sp.model, 
                         disp.par, div, IFM.bounding.box = NULL,
                         fixed.landscape.psi0 = F, psi0.var.name,
                         ...) 
{
  
  ##Spatial probability of colonisation or extinction
  get.phis <- function(phiParams, X.phi,
                       psis.landscape){
    
    if (sp.model == "IFM"){
      ##“baseline” colonization probability
      base_delta <- plogis(X.phi %*% phiParams)
      
      alpha <- disp.par
      
      if(is.null(disp_kernel_times_psi0) ){
        phis <- rcpp_phi(distance_neighs, index_neighs,
                         base_delta, psis.landscape, alpha, div)
      } else{
        phis <- rcpp_y(disp_kernel_times_psi0, base_delta)
      }
      
      ##Eliminate boundary probabilities (i.e., 0 or 1)
      phis[phis==0] <- 10^-16
      phis[phis==1] <- 1 - 10^-16
      
    } else if(sp.model == "BRM"){
      if(is.null(sp_autocovariate) ){
        sp_autocovariate <- sapply(index_neighs,
                                   function(x) sum(psis.landscape[x], na.rm=T)/length(x) ) 
      }
      X.phi <- cbind(X.phi, sp_autocovariate)
      phis <- plogis(X.phi %*% phiParams)
    }
    
    return(phis)
  }
  
  conditional.prob <- function(psis, detParms, y.t.arr, t){
    
    
    prob.no.obs <- 1 - array(V.itjk %*% detParms, c(J, nY, M))
    prob.no.obs <- aperm(prob.no.obs, c(3:1))
    
    cum.not.obs <- apply( prob.no.obs[,t,], 1, prod)
    prob.present.not.obs <- psis*cum.not.obs
    
    conditional.prob <- prob.present.not.obs/ ((1-psis) + prob.present.not.obs) #Mackenzy 2006 pg 97
    
    detected <- apply(y.t.arr, 1, any) #vector length(M) species detected, time t
    
    conditional.prob[detected] <- 1
    
    conditional.prob
    
  }
  
  
  K <- 1
  
  
  ##Design mats for all landscape cells
  landscape.designMats <- unmarked:::getDesign(
    data, 
    formula = as.formula(paste(unlist(formula), collapse = " ")),
    na.rm=F)
  W.landscape <- landscape.designMats$W
  X.gam.landscape <- landscape.designMats$X.gam
  X.eps.landscape <- landscape.designMats$X.eps
  
  
  ##Design mats for study cells
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
  if (sp.model == "BRM"){
    if (spatial.col == T){
      nGP <- nGP + 1
      gamParms <- c(gamParms, "col_autocovariate")
    }
    if (spatial.ext ==T){
      nEP <- nEP + 1
      epsParms <- c(epsParms, "ext_autocovariate")
    }
  }
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
  
  ####  Spatial module  ####
  ##Index study cells in landscape cells df
  index_study_in_landscape <- setdiff(1:nSites, designMats$removed.sites)
  coord <- data@siteCovs[,c("lon", "lat")]
  #' ####################################################################### #
  #'  !To do: Raise error message if "lon" "lan" no present
  #'  !Check memory: memory: nrow(coord)*nrow(coord_study)*8*2/1000000 MB
  #'  !Calculate max_dist
  #' ####################################################################### #
  coord_study <- coord[index_study_in_landscape,]
  ##Adjacency list with neighbours indexes
  index_neighs <- list()
  ##Adjacency list with neighbours distances
  distance_neighs <- list()  
  for (i in 1:M){
    
    ##Neighbours distances
    points_distances <- pointDistance(coord_study[i ,c("lon", "lat")],
                                      coord[,c("lon", "lat")],
                                      lonlat=F)
    points_distances <- points_distances/1000
    ##Eliminate same cell distances
    points_distances[index_study_in_landscape[i]] <- NA
    
    ##Neighbours to consider
    if (sp.model=="BRM"){
      neigh <- which(points_distances<disp.par)
    } else if (!is.null(IFM.bounding.box) ){
      neigh <- which(points_distances<IFM.bounding.box)
    }else{
      neigh <- (1:length(points_distances) )[-index_study_in_landscape[i]]
    }
    
    index_neighs[[i]] <- neigh
    distance_neighs[[i]] <- points_distances[neigh]
  }
  
  ##Precaculate spatial autocovariate (BRM) or dispersal kernel times psi0 (IFM)
  if (fixed.landscape.psi0){
    
    psis.landscape <- data@siteCovs[,psi0.var.name]
    
    if (sp.model == "IFM"){
      
      disp_kernel_times_psi0 <- list()
      for (i in 1:M ){
        product <-
          exp(-distance_neighs[[i]]/disp.par)*
          psis.landscape[index_neighs[[i]]]/div
        disp_kernel_times_psi0[[i]] <- product 
      }
    } else if (sp.model == "BRM") {
      sp_autocovariate <- sapply(index_neighs,
                                 function(x) sum(psis.landscape[x], na.rm=T)/length(x)  )
    }
    
  } else{
    disp_kernel_times_psi0 <- NULL
    sp_autocovariate <- NULL
  }
  
  #####
  
  nll <- function(params) {
    
    colParams <- params[(nSP + 1):(nSP + nGP)]
    extParams <- params[(nSP + nGP + 1):(nSP + nGP + nEP)]
    detParams <- params[(nSP + nGP + nEP + 1):nP]
    
    psis <- plogis(W.i %*% params[1:nSP])
    
    cond.prob <- conditional.prob(psis, detParams, y.arr[,1,],t = 1)
    
    psis.landscape <-  as.vector( plogis(W.landscape %*% params[1:nSP]) )
    psis.landscape[index_study_in_landscape] <- cond.prob
    
    if (spatial.col == T){
      gammas <- get.phis(colParams, 
                         X.it.gam[seq(1, M * (nY-1), by = nY-1),
                                  , drop = FALSE],
                         psis.landscape)
      #Remember phis: array(NA, c(2, 2, nY - 1, M))
      phis[1, 1, 1 , ] <- 1-gammas 
      phis[2, 1, 1, ] <- gammas
    } else{
      phis[, 1, , ] <- plogis(X.gam %*% colParams)
    }
    if (spatial.ext == T){
      epsilons <- 1-get.phis(extParams, 
                             X.it.eps[seq(1, M * (nY-1), by = nY-1),
                                      , drop = FALSE],
                             psis.landscape)
      phis[1, 2, 1, ] <- epsilons
      phis[2, 2, 1, ] <- 1-epsilons
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
  if (!fixed.landscape.psi0){
    psis.landscape <-  plogis(W.landscape %*% mle[1:nSP])
  }
    
  
  colParams <- mle[(nSP + 1):(nSP + nGP)]
  extParams <- mle[(nSP + nGP + 1):(nSP + nGP + nEP)]
  detParams <- mle[(nSP + nGP + nEP + 1):nP]
  
  if (spatial.col == T){
    gammas <- get.phis(colParams, 
                       X.it.gam[seq(1, M * (nY-1), by = nY-1),
                                , drop = FALSE],
                       psis.landscape)
    #Remember phis: array(NA, c(2, 2, nY - 1, M))
    phis[1, 1, 1 , ] <- 1-gammas 
    phis[2, 1, 1, ] <- gammas
  } else{
    phis[, 1, , ] <- plogis(X.gam %*% colParams)
  }
  
  if (spatial.ext == T){
    epsilons <- 1-get.phis(extParams, 
                           X.it.eps[seq(1, M * (nY-1), by = nY-1),
                                    , drop = FALSE],
                           psis.landscape)
    phis[1, 2, 1, ] <- epsilons
    phis[2, 2, 1, ] <- 1-epsilons
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


#### Test ####

##1. Substitute default colext functions for modified spatial colext functions
environment(spatial.colext) <- asNamespace('unmarked')
# assignInNamespace("spatial.colext", spatial.colext, ns = 'unmarked')

environment(spatial.colext.fit) <- asNamespace('unmarked')
# assignInNamespace("colext", colext2, ns = 'unmarked')

####  Load data ####

# Dir.Base <- "C:/Users/uriso/OneDrive - CREAF/Paper_colext/S4"
# 
# ##Study cells (i.e. surveyed cells)
# load(file.path(Dir.Base, "study_cells.Rdata")) ##df: cell_ID|lon|lat
# 
# ##Landscape cells (i.e. surveyed an unsurveyed cells)
# load(file.path(Dir.Base, "landscape_cells.Rdata")) ##df: cell_ID|lon|lat|sdm_cbba2
# 
# ##Detection data Luscinia megarhynchos
# load(file.path(Dir.Base, "detection_data.Rdata"))##df: season-observation (major minor order)
# 
# 
# ####  Create unmarked multiframe with landscape info  ####
# 
# ##Detection data
# landscape_detection_data <- matrix(nrow=nrow(landscape_cells), 
#                                    ncol=4)
# 
# index_study_in_landscape <- match(study_cells$cell_id, landscape_cells$cell_id)
# 
# landscape_detection_data[index_study_in_landscape,] <- detection_data
# 
# 
# ## Create unmarked multiframe (umf)
# 
# simUMF <- unmarkedMultFrame(
#   y = landscape_detection_data,
#   siteCovs = landscape_cells,
#   numPrimary=2)
# 
# 
# ####Fit model ####
# a<-spatial.colext(psiformula = ~sdm_cbba2, # First-year occupancy
#                data = simUMF, fixed.landscape.psi0 = T, 
#                psi0.var.name = "sdm_cbba2",
#                spatial.col = T, spatial.ext = T, sp.model= "BRM", 
#                disp.par= 2)
# 
# a@negLogLike
# 
# b<-spatial.colext(psiformula = ~sdm_cbba2, # First-year occupancy
#                   data = simUMF, 
#                   spatial.col = T, spatial.ext = T, sp.model= "BRM", 
#                   disp.par= 2)
# 
# b@negLogLike
# 
# c<-spatial.colext(psiformula = ~sdm_cbba2, # First-year occupancy
#                   data = simUMF, fixed.landscape.psi0 = T, 
#                   psi0.var.name = "sdm_cbba2",
#                   spatial.col = T, spatial.ext = T, sp.model= "IFM", 
#                   disp.par= 2, IFM.bounding.box = 5)
# 
# c@negLogLike
# 
# tic()
# d<-spatial.colext(psiformula = ~sdm_cbba2, # First-year occupancy
#                   data = simUMF, 
#                   spatial.col = T, spatial.ext = T, sp.model= "IFM", 
#                   disp.par= 2, IFM.bounding.box = 5)
# toc()
# d@negLogLike
# 
# tic()
# e<-spatial.colext(psiformula = ~sdm_cbba2, # First-year occupancy
#                   data = simUMF, 
#                   spatial.col = T, spatial.ext = T, sp.model= "IFM", 
#                   disp.par= 2, IFM.bounding.box = 20)
# toc()
# e@negLogLike
# 
# 
