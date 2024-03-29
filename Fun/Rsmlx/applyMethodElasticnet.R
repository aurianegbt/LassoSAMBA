applyMethodElasticnet <- function(Y,X,omega,cov0,
                             stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.9,
                             criterion="BIC",ncrit=20,covariate.model=NULL,
                             printFrequencySS = TRUE,p.name=NULL){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètr
  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC 
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }

  cov.names = colnames(X)
  tparam.names = colnames(Y)

  if(!is.null(covariate.model)){
    savedSelection = setNames(as.numeric(covariate.model),cov.names)
    prevSelection = covariate.model
    if(!is.matrix(Y)){
      Yaux <- as.matrix(Y)
    }else{Yaux=Y}
    if(!is.matrix(X)){
      Xaux <- as.matrix(X)
    }else{Xaux=X}
    
    Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
    
    if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
    Ywh <- Y %*% rootInvOmega
    Xwh <- kronecker(t(rootInvOmega),Xsc)
    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,prevSelection]
      oldCriterion =critFUN(lm(Ywh ~ Xkeep))
    }
    cat("\nSearch of a elastic net selection improving the ",criterion," criterion (max ",ncrit," searchs) for ",p.name," :\n ")
    cat(paste0("  ▶ old Criterion ",criterion," : ",round(oldCriterion,digits=2)))
  }

  resSelection = elasticnetSelection(Y,X,omega,cov0,stabilitySelection,nfolds,nSS,thresholdsSS,criterion,printFrequencySS = printFrequencySS)
  
  selection = resSelection$selection
  newCriterion = resSelection$criterion
  stop=resSelection$stop
  if(stop){
    cat( " - not enough covariates in search scope.\n")
  }
  if(!stop & !is.null(covariate.model)){
    tested = 1
    cat("\n    - new Criterion ",tested, " : ",round(newCriterion,digits=2))
    cat("\n           > parameter values : ",
        paste0(names(resSelection$param),"=",
               sapply(resSelection$param,function(x){
                 10**floor(log10(x))*round(x/(10**floor(log10(x))),digits=2)
                 }),collapse=","))
    flag = newCriterion < oldCriterion
    while(tested < ncrit & !flag){
      resSelection = elasticnetSelection(Y,X,omega,cov0,stabilitySelection,nfolds,nSS,thresholdsSS,criterion,printFrequencySS = printFrequencySS)
      selection = resSelection$selection
      newCriterion = resSelection$criterion

      tested = tested + 1
      cat("\n    - new Criterion ",tested, " : ",round(newCriterion,digits=2))
      cat("\n           > parameter values : ",
          paste0(names(resSelection$param),"=",
                 sapply(resSelection$param,function(x){
                   10**floor(log10(x))*round(x/(10**floor(log10(x))),digits=2)
                   }),collapse=","))
      flag = newCriterion < oldCriterion
    }

    if(!flag){
      selection=savedSelection
    }
    cat("\n")
  }

  model.list = modelFromSelection(Y,X,selection)

  r <- list()
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name))
}
