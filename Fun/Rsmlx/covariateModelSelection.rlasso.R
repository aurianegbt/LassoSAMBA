covariateModelSelection.rlasso <- function(nfolds = 5,
                                          alpha = 1,
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
                                          ncrit=20){
  # Simulate Individual Parameters and setup parameters
  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC 
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }
  
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
  
  Sigma=Rsmlx::getEstimatedCovarianceMatrix()$cov.matrix
  colnames(Sigma) <- colnames(Rsmlx::getEstimatedCovarianceMatrix()$cor.matrix)
  rownames(Sigma) <- rownames(Rsmlx::getEstimatedCovarianceMatrix()$cor.matrix)
  # # param.names for order in r :
  N = length(unique(Y$id))
  Y.mat = lapply(split(Y,f = Y$rep),FUN=function(x){x[,-which(colnames(x) %in% c("rep","id"))]})
  
  X.mat  = covariates[,setdiff(colnames(covariates),"id")]
  runREP = foreach(k = 1:nrep) %dopar% {
    source("Fun/Rsmlx/applyMethodLasso.R")
    source("Fun/Rsmlx/lassoSelection.R")
    source("Fun/Rsmlx/modelFromSelection.R")
    Yk.mat = as.matrix(Y.mat[[k]])
    (sapply(param.names[which(indvar)],FUN=function(x){
      applyMethodLasso(Yk.mat[,stringr::str_detect(colnames(Yk.mat),x),drop=F],
                       X.mat,Sigma[x,x],alpha,cov0.list[[x]],
                       stabilitySelection=FALSE,nfolds=nfolds,
                       thresholdsSS=NULL,
                       criterion=criterion,ncrit=ncrit,
                       covariate.model=NULL,
                       printFrequencySS=FALSE,p.name=x)$res}))
  }
  
  # cov0 compute prior 
  Y.mat =  sapply(Y[,-c(1,2)],function(x){rowMeans(matrix(x,nrow=N))})
  
  r.AUX = Reduce("+",runREP)
  resSelection = t(r.AUX/nrep)
  
  if(length(thresholdsSS)!=1){
    sel=t(sapply(param.names[which(indvar)],FUN=function(p){
      selection.p <- resSelection[p,]
      Yp.mat = Y.mat[,stringr::str_detect(colnames(Y.mat),p),drop=F]
      
      critThresholds = setNames(sapply(thresholdsSS,FUN=function(tSS){
        res = selection.p>=tSS
        
        if(all(!res)){
          crit.tSS = critFUN(lm(Yp.mat ~ NULL))
        }else{
          Xkeep = as.matrix(X.mat[,res])
          crit.tSS = critFUN(lm(Yp.mat~Xkeep))
        }
        return(crit.tSS)
      }),thresholdsSS)
      
      thresholdsFinal = min(thresholdsSS[which(critThresholds==min(critThresholds))])
      
      res.p = selection.p >=thresholdsFinal
      mode(res.p) <- "integer"
      return(res.p)
    }))
  }else{
    sel = apply(resSelection>=thresholdsSS,MARGIN=2,FUN=as.numeric)
    rownames(sel) <- rownames(resSelection)
  }
  
  
  r.var = lapply(param.names[which(indvar)],
                 FUN=function(x){
                   mod = modelFromSelection(Y.mat[,stringr::str_detect(colnames(Y.mat),x),drop=F],X.mat,setNames(sel[x,],colnames(sel)))
                   list(model=mod,res=setNames(as.logical(sel[x,]),colnames(sel)),cov0=cov0.list[[x]],p.name=x)})
    
  
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
              add.covariate = trs, sp = sp.df, tr0 = tr0, r.cov0 = r.cov0))
}

