lassoSelection <- function(Y,X,Sigma=NULL,alpha=1,cov0.list=NULL, # cov0.list dans le même ordre que Y !!!!
                           stabilitySelection=TRUE,nfolds=5,nSS=1000,thresholdsSS=0.90,
                           printFrequencySS = FALSE){
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }

  cov.names = colnames(X)
  param.names = names(cov0.list)

  Ysc <- scale(apply(Y,2,FUN=as.numeric))
  Xsc <- scale(apply(X,2,FUN=as.numeric))

  if(!is.null(Sigma)){ rootInvSigma = solve(chol(Sigma)) }else{ rootInvSigma = diag(ncol(Y)) }
  Ywh <- as.numeric( Ysc %*% rootInvSigma)
  Xwh <- kronecker(t(rootInvSigma),Xsc)


  if(is.null(cov0.list)){
    exclude = NULL
  }else{
    if(ncol(Y)==1){
      if(is.list(cov0.list)){
        exclude = which(cov.names %in% cov0.list[[1]])
      }else{
        exclude = which(cov.names %in% cov0.list)
      }
    }else{
      exclude = numeric(0)
      ncov.tot = length(cov.names)
      for(j in 1:length(cov0.list)){
        exclude <- c(exclude,ncov.tot*(j-1)+which(cov.names %in% cov0.list[[j]]))
      }
    }
  }

  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ # sinon glmnet ne fonctionne pas
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1 # je décide qu'on sélectionne la seule covariable présente s'il y en a une (au pire on la retirera avec le test de p.valeur)
  }else{
    resCV = glmnet::cv.glmnet(Xwh, Ywh, alpha = alpha , nfolds=nfolds , exclude=exclude)

    lambda=resCV$lambda.min

    if(!stabilitySelection){
      res = glmnet::glmnet(Xwh, Ywh, alpha=alpha, lambda=lambda, exclude=exclude)
      selection = tabulate(Matrix:::which(res$beta != 0),ncol(Xwh))
    }else{
      selection = data.frame()
      for(i in 1:nSS){
        indSampled = sample(1:length(Ywh), floor(length(Ywh)/2))
        res = glmnet::glmnet(Xwh[indSampled,], Ywh[indSampled], alpha=alpha, lambda=lambda, exclude=exclude)

        selection <- rbind(selection,tabulate(Matrix::which(res$beta != 0),ncol(Xwh)))
      }
      resSelection = colSums(selection)/nSS

      if(printFrequencySS){
        cat(paste0("Selection frequency for ",nSS," runs :\n"))
        for.print <- matrix(resSelection,nrow=ncol(Y),byrow=TRUE)
        colnames(for.print) <- cov.names
        rownames(for.print)<- colnames(Y)
        print(for.print)
      }

      resSelection[resSelection>=thresholdsSS] <- 1
      resSelection[resSelection<thresholdsSS] <- 0

      selection <- resSelection
    }
  }

  if(printFrequencySS){
    cat("\nCovariate selection results :\n")
  }
  selection <- matrix(selection,nrow=ncol(Y),byrow=TRUE)
  colnames(selection) <- cov.names
  rownames(selection)<- colnames(Y)
  return(selection)
}
