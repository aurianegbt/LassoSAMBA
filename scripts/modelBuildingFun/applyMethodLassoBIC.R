applyMethodLassoBIC <- function(Y,X,omega,cov0,
                             alpha=1,
                             criterion="BIC",
                             covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1,
                             iter=1){
  to.cat = c()
  
  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC 
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }
  
  cov.names = colnames(X)
  tparam.names = colnames(Y)
  
  if(!is.matrix(Y)){
    Yaux <- as.matrix(Y)
  }else{Yaux=Y}
  if(!is.matrix(X)){
    Xaux <- as.matrix(X)
  }else{Xaux=X}
  
  Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
  
  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  Ywh <- Yaux %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)
  colnames(Xwh) <- cov.names
  
  if(!is.null(covariate.model)){
    if(any(!(cov.names %in% names(covariate.model)))){
      savedSelection <- setNames(rep(0,length(cov.names)),cov.names)
      savedSelection[names(covariate.model)] <- as.numeric(covariate.model)
    }else{
      savedSelection = setNames(as.numeric(covariate.model),names(covariate.model))
    }
    prevSelection = covariate.model
    
    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(prevSelection)[which(prevSelection)]]
      oldCriterion = critFUN(lm(Ywh ~ Xkeep))
    }
    to.cat <- c(to.cat,paste0("\n Lasso selection, calibrated using BIC method, improving the ",criterion," criterion for ",p.name," :\n "))
    to.cat <- c(to.cat,paste0("       -> Old Criterion : ",round(oldCriterion,digits=2)),"\n")
  }
  
  if(is.null(cov0)){
    exclude = NULL
  }else{
    exclude = which(cov.names %in% cov0)
  }
  
  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))
    
    to.cat.here = ""
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ 
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1
    if(all(!as.logical(selection))){
      newcriterion = critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
      newcriterion = critFUN(lm(Ywh~Xkeep))
    }
    
    to.cat.here = ""
  }else{
    
    
    
    fit <- glmnet::glmnet(Xwh,Ywh,alpha=alpha,exclude = exclude)
    
    tLL <- fit$nulldev - deviance(fit) 
    # 2*(loglike_sat -loglike(Null) - 2*(loglike_sat - loglike) =  2loglik - 2loglike(Null) = 2LL - cst 
    k <- fit$df
    n <- fit$nobs
    
    fit.BIC <- log(n)*k - tLL
    argmax_id = which.min(fit.BIC)
    
    coef.final = coef(fit,s=fit$lambda[which.min(fit.BIC)])
    
    selection = setNames(as.numeric(coef.final[-1,1]!=0),names(coef.final[-1,1]))
    if(all(!as.logical(selection))){
      newcriterion = critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
      newcriterion = critFUN(lm(Ywh~Xkeep))
    }
    
    
    
    to.cat.here =  paste0("\n              > parameter values : ",
                          paste0("lambda=",round(fit$lambda[which.min(fit.BIC)],digits=3)))
  }
  
  if(newcriterion >= oldCriterion){
    to.cat <- c(to.cat,paste0("        No model improving the criterion as been find, the previous covariate model is kept."))
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,to.cat.here)
  }
  
  model.list = modelFromSelection(Y,X,selection)
  
  
  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat))
}