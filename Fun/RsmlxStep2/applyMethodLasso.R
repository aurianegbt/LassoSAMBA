applyMethodLasso <- function(Y,X,Sigma=NULL,alpha=1,cov0.list=NULL,
                             stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.9,
                             critMV="BIC",ncrit=20,covariate.model=NULL){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre

  cov.names = colnames(X)
  param.names = names(cov0.list)
  n.param = ncol(Y)
  tparam.names = colnames(Y)

  if(!is.null(covariate.model)){
    covariate.model = lapply(covariate.model[param.names],FUN=function(x){x[cov.names]})
    savedSelection = t(as.data.frame(lapply(covariate.model,FUN=as.numeric)))
    colnames(savedSelection) <- cov.names
    prevSelection = unlist(covariate.model,use.names = FALSE)
    if(!is.matrix(Y)){
      Yaux <- as.matrix(Y)
    }else{Yaux=Y}
    if(!is.matrix(X)){
      Xaux <- as.matrix(X)
    }else{Xaux=X}
    Ysc <- scale(apply(Yaux,2,FUN=as.numeric))
    Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
    if(!is.null(Sigma)){ rootInvSigma = solve(chol(Sigma)) }else{ rootInvSigma = diag(ncol(Y)) }
    Ywh <- as.numeric( Ysc %*% rootInvSigma)
    Xwh <- kronecker(t(rootInvSigma),Xsc)
    if(all(!prevSelection)){
      eval(parse(text=(paste0("oldCriterion =",critMV,"(lm(Ywh ~ NULL))"))))
    }else{
      Xkeep = Xwh[,prevSelection]
      eval(parse(text=(paste0("oldCriterion =",critMV,"(lm(Ywh ~ Xkeep))"))))
    }
    cat("Search of a lasso selection improving the ",critMV," criterion (max ",ncrit," searchs) :\n ")
    cat(paste0("  ▶ old Criterion ",critMV," : ",round(oldCriterion,digits=2)))
  }

  resSelection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS,critMV)
  selection = resSelection$selection
  newCriterion = resSelection$criterion
  if(!is.null(covariate.model)){
    tested = 1
    cat(paste0("\n    - new Criterion ",tested, " : ",round(newCriterion,digits=2)))
    flag = newCriterion < oldCriterion
    while(tested < ncrit & !flag){
      resSelection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS,critMV)
      selection = resSelection$selection
      newCriterion = resSelection$criterion

      tested = tested + 1
      cat(paste0("\n    - new Criterion ",tested, " : ",round(newCriterion,digits=2)))
      flag = newCriterion < oldCriterion
    }

    if(!flag){
      selection=savedSelection
    }
    cat("\n")
  }

  model.list = modelFromLassoSelection(Y,X,selection,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS)

  r <- list()
  if(n.param==1){
    return(list(model=model.list,res=selection,cov0=cov0.list,p.name=param.names))
  }else{
    for(j in 1:n.param){
      r[[j]] <- list(model = model.list[[j]], res = selection[j,],cov0=cov0.list[[j]],p.name=param.names[j])
    }
    return(r)
  }
}
