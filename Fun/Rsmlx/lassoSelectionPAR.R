lassoSelection <- function(Y,X,
                           omega=NULL,
                           alpha=1,
                           cov0=NULL, # cov0.list dans le même ordre que Y !!!!
                           stabilitySelection=TRUE,
                           nfolds=5,
                           nSS=1000,
                           thresholdsSS=seq(0.5,0.95,0.05), # can be a vector of thresholdsSS
                           criterion="BIC",
                           printFrequencySS = FALSE,
                           MSE=FALSE){
  # depends on for each 
  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC 
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }
  
  if(!is.matrix(Y)){
    Y <- as.matrix(Y)
  }
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }

  cov.names = colnames(X)

  Xsc <- scale(apply(X,2,FUN=as.numeric))

  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  Ywh <- Y %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)


  if(is.null(cov0)){
    exclude = NULL
  }else{
    exclude = which(cov.names %in% cov0)
  }

  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))
    
    lambda=NULL
    thresholdsFinal=NULL
    stop=TRUE
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ # sinon glmnet ne fonctionne pas
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1 # je décide qu'on sélectionne la seule covariable présente s'il y en a une (au pire on la retirera avec le test de p.valeur)
    
    lambda=NULL
    thresholdsFinal = NULL
    stop=TRUE
  }else{
    resCV = glmnet::cv.glmnet(Xwh, Ywh, alpha = alpha , nfolds=nfolds , exclude=exclude)

    lambda=resCV$lambda.min

    if(!stabilitySelection){
      res = glmnet::glmnet(Xwh, Ywh, alpha=alpha, lambda=lambda, exclude=exclude)
      selection = tabulate(Matrix:::which(res$beta != 0),ncol(Xwh))
    }else{
      selection =
        foreach(i = 1:nSS,.combine = "rbind") %dopar% {
        indSampled = sample(1:length(Ywh), floor(length(Ywh)/2))
        res = glmnet::glmnet(Xwh[indSampled,], Ywh[indSampled], alpha=alpha, lambda=lambda, exclude=exclude)

       tabulate(Matrix::which(res$beta != 0),ncol(Xwh))
      }
      resSelection = colSums(selection)/nSS

      if(printFrequencySS){
        cat(paste0("Selection frequency for ",nSS," runs :\n"))
        for.print <- matrix(resSelection,nrow=ncol(Y),byrow=TRUE)
        colnames(for.print) <- cov.names
        rownames(for.print)<- colnames(Y)
        print(for.print)
      }
      
      resThresholds = lapply(thresholdsSS,function(t){
        resSelection[resSelection>=t] <- 1
        resSelection[resSelection<t] <- 0
        return(resSelection)
      })
      
      if(MSE){
        critThresholds = setNames(sapply(resThresholds,function(res){
          if(all(!as.logical(res))){
            model=lm(Ywh ~ NULL)
          }else{
            Xkeep = Xwh[,as.logical(res)]
            model = lm(Ywh~Xkeep)
          }
          
          criterion = mean(summary(model)$residuals^2)
          return(criterion)
        }),thresholdsSS)
      }else{
        critThresholds = setNames(sapply(resThresholds,function(res){
          if(all(!as.logical(res))){
            criterion = critFUN(lm(Ywh ~ NULL))
          }else{
            Xkeep = Xwh[,as.logical(res)]
            criterion = critFUN(lm(Ywh~Xkeep))
          }
          return(criterion)
        }),thresholdsSS)
      }
      thresholdsFinal = min(thresholdsSS[which(critThresholds==min(critThresholds))])
      
      resSelection[resSelection>=thresholdsFinal] <- 1
      resSelection[resSelection<thresholdsFinal] <- 0
      
      selection <- resSelection
    }
    stop=FALSE
  }
  
  if(all(!as.logical(selection))){
    criterion = critFUN(lm(Ywh ~ NULL))
  }else{
    Xkeep = Xwh[,as.logical(selection)]
    criterion = critFUN(lm(Ywh~Xkeep))
  }

  selection <- setNames(as.logical(selection),cov.names)
  if(stabilitySelection){
    return(list(selection=selection,criterion=criterion,param=c(lambda=lambda,thresholds=thresholdsFinal),stop=stop))
  }else{
    return(list(selection=selection,criterion=criterion,param=c(lambda=lambda),stop=stop))
  }
}
