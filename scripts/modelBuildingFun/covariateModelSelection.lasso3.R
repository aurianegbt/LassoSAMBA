covariateModelSelection.lasso <- function(nfolds = 5,
                                          alpha = 1,
                                          stabilitySelection = TRUE,
                                          nSS=1000,
                                          thresholdsSS=0.80,
                                          covFix = NULL,
                                          pen.coef=NULL,
                                          weight=1,
                                          paramToUse="all",
                                          eta=NULL,
                                          p.max=1,
                                          sp0=NULL,
                                          covariate.model=NULL,
                                          criterion="BIC",
                                          ncrit=20,
                                          lambda.grid=NULL,
                                          printFrequencySS = FALSE,
                                          iter=1){
  # Simulate Individual Parameters and setup parameters
  
  sp.df <- Rsmlx:::mlx.getSimulatedIndividualParameters()
  if (is.null(sp.df$rep))
    sp.df$rep <- 1
  if (!is.null(sp0)) {
    nrep0 <- max(sp0$rep)
    sp.df$rep <- sp.df$rep + nrep0
    dn <- setdiff(names(sp.df), names(sp0))
    sp0[dn] <- sp.df[dn]
    sp.df <- rbind(sp0, sp.df)
  }
  nrep <- max(sp.df$rep)
  ind.dist <- Rsmlx:::mlx.getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  n.param <- length(param.names)
  cov.info <- Rsmlx:::mlx.getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  tcov.names <- NULL
  covariates <- cov.info$covariate
  cov.cat <- cov.names[cov.types == "categorical"]
  covariates[cov.cat] <- lapply(covariates[cov.cat], as.factor)
  indvar <- Rsmlx:::mlx.getIndividualParameterModel()$variability$id
  indvar[setdiff(param.names, paramToUse)] <- FALSE
  
  cov.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
  r <- res <- r.cov0 <- list()
  eps <- 1e-15

  # Y transformation
  Y=sp.df[,c("rep","id",names(indvar)[which(indvar==TRUE)])]

  # Update eta cov0 and pweight
  eta.list =list()
  cov0.list = list()
  pweight.list = list()
  for(j in 1:n.param){
    dj <- ind.dist[j]
    nj <- names(dj)
    if (indvar[j]) {
      if (tolower(dj) == "lognormal") {
        Y[,nj] <- log(Y[,nj] + eps)
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("log.", nj)
      }else if (tolower(dj) == "logitnormal") {
        Y[,nj] <- log((Y[,nj] + eps)/(1 - Y[,nj] + eps))
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("logit.", nj)
      }else if (tolower(dj) == "probitnormal") {
        Y[,nj] <- qnorm(Y[,nj])
        colnames(Y)[which(colnames(Y)==nj)]  <- paste0("probit.", nj)
      }
      if (!is.null(eta)) {
        eta.list <- append(eta.list,list(eta[paste0("eta_", nj)]))
      }else {eta.list <- append(eta.list,list(NULL))}
      names(eta.list)[length(eta.list)] <- nj
      if (length(covFix) > 0){
        cmj <- cov.model[[nj]][covFix]
        cov0.list <- append(cov0.list,list(names(which(!cmj))))
      }else {
        cov0.list <-  append(cov0.list,list(NULL))
      }
      names(cov0.list)[length(cov0.list)] <- nj
      pweight.list <- append(pweight.list,list(weight[nj, ]))
      names(pweight.list)[length(pweight.list)] <- nj
    }
  }
  cov0.list <- updateCov0(Y, eta.list, covariates, p.max, covFix,
                            pen.coef, pweight.list, cov0.list)
  
  Sigma=diag(Rsmlx:::mlx.getEstimatedPopulationParameters()[paste0("omega_",param.names[which(indvar)])]**2)
  colnames(Sigma) <- param.names[which(indvar)]
  rownames(Sigma) <- param.names[which(indvar)]
  # # param.names for order in r :
  N = length(unique(Y$id))
  Y.mat = sapply(Y[,-c(1,2)],function(x){rowMeans(matrix(x,nrow=N))}) #1 : rep 2 : id
  X.mat  = covariates[,setdiff(colnames(covariates),"id")]
  
  
  if(iter==1){
    pathToSavePlot = paste0(dirname(lixoftConnectors::getProjectSettings()$directory),"/CalibrationPlot")
    if(!dir.exists(pathToSavePlot)){dir.create(pathToSavePlot)}
    pathToSavePlot <- paste0(pathToSavePlot,"/calibrationPlot_iter",iter,"_")
    
    outputs = setNames(foreach(p = names(indvar)[which(indvar)],.packages=c("sharp","ggplot2","grDevices"),.export = "sharpCalibration") %dopar% { 
      res = sharpCalibration(Y.mat[,stringr::str_detect(colnames(Y.mat),p),drop=F],
                             X.mat,
                             Sigma[p,p],
                             nfolds,alpha,nSS,p.name=p,iter=iter,
                             n_cores = floor(parallel::detectCores()/sum(which(indvar))))
      
      ggsave(plot=res$plot, filename =paste0(pathToSavePlot,p,".jpeg"),
             height=1500,width=3000,units = "px",device=grDevices::jpeg)
      
      return(res)
    },names(indvar)[which(indvar)])
    
    thresholdsSS = setNames(lapply(outputs,FUN=function(x){unname(x$thresholdsSS)}),names(indvar)[which(indvar)])
    
    r.var <- lapply(param.names[which(indvar)],FUN=function(x){
      applyMethodLasso(Y.mat[,stringr::str_detect(colnames(Y.mat),x),drop=F],
                       X.mat,Sigma[x,x],alpha,cov0.list[[x]],
                       stabilitySelection,nfolds,nSS,outputs[[x]]$thresholdsSS,
                       criterion,ncrit,covariate.model[[x]],
                       lambda.grid=outputs[[x]]$lambda.grid,
                       printFrequencySS,p.name=x)})
    
  }else{
    r.var <- lapply(param.names[which(indvar)],FUN=function(x){
      applyMethodLasso(Y.mat[,stringr::str_detect(colnames(Y.mat),x),drop=F],
                       X.mat,Sigma[x,x],alpha,cov0.list[[x]],
                       stabilitySelection,nfolds,nSS,thresholdsSS[[x]],
                       criterion,ncrit,covariate.model[[x]],
                       lambda.grid=NULL,
                       printFrequencySS,p.name=x)}) 
  }
  
  r <- res <- r.cov0 <- list()
  for(j in 1:n.param){
    nj=param.names[j]
    if(indvar[j]){
      ind = which(unlist(lapply(r.var,FUN=function(x){x$p.name}))==nj)
      r[[j]] <- r.var[[ind]]

      res[[j]] <- r[[j]]$res
      names(res)[j] <- param.names[j]
      r.cov0 <- append(r.cov0,list(r[[j]]$cov0))
      names(r.cov0)[length(r.cov0)] <- param.names[j]

    }else{
      r[[j]] <- list(model="fixed")
      res[[j]] <- "none"
      names(res)[j] <- param.names[j]
      r[[j]]$p.name <- nj
    }
  }
  
  e <- as.data.frame(lapply(r[indvar], function(x) {
    x$model$residuals
  }))
  e.names <- unlist(lapply(r[indvar], function(x) {
    x$p.name
  }))
  names(e) <- paste0("eta_", e.names)
  if (!is.null(sp.df["id"]))
    e <- cbind(sp.df["id"], e)
  if (!is.null(sp.df["rep"]))
    e <- cbind(sp.df["rep"], e)

  covariate.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
  covariate <- Rsmlx:::mlx.getCovariateInformation()$covariate
  js <- 0
  trs <- list()
  tr0 <- NULL
  for (k in (1:n.param)) {
    if (!identical(res[[k]], "none")) {
      covariate.model[[k]][1:length(covariate.model[[k]])] <- FALSE
      if (indvar[k]) {
        ck <- attr(r[[k]]$model$terms, "term.labels")
        if (length(ck) > 0) {
          for (j in (1:length(ck))) {
            ckj <- ck[j]
            if (identical(substr(ckj, 1, 4), "log.")) {
              js <- js + 1
              ckj.name <- sub("log.", "", ckj)
              covkj <- covariate[[ckj.name]]
              lckj <- paste0("l", ckj.name)
              tr.str <- paste0(lckj, " = \"log(", ckj.name,
                               "/", signif(mean(covkj), digits = 2),
                               ")\"")
              trs[[js]] <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",
                                  tr.str, ")")
              tr0 <- unique(c(tr0, ckj.name))
              covariate.model[[k]][lckj] <- TRUE
            }
            else {
              covariate.model[[k]][ckj] <- TRUE
            }
          }
        }
      }
    }
  }
  res <- Rsmlx:::formatCovariateModel(res)

  return(list(model = covariate.model, residuals = e, res = res,
              add.covariate = trs, sp = sp.df, tr0 = tr0, r.cov0 = r.cov0,thresholdsSS=thresholdsSS))
}

