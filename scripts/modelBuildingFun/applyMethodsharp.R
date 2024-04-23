applyMethodsharp <- function(Y,X,omega,cov0,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             criterion="BIC",ncrit=20,covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre
  
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
    savedSelection = setNames(as.numeric(covariate.model),names(covariate.model))
    prevSelection = covariate.model
    
    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,prevSelection]
      oldCriterion =critFUN(lm(Ywh ~ Xkeep))
    }
    to.cat <- c(to.cat,paste0("\n Lasso selection, calibrated using sharp method, improving the ",criterion," criterion for ",p.name," :\n "))
    to.cat <- c(to.cat,paste0("       -> Old Criterion : ",round(oldCriterion,digits=2)),"\n")
  }
  
  if(is.null(cov0)){
    exclude = NULL
  }else{
    exclude = which(cov.names %in% cov0)
  }
  
  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ 
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1
  }else{
    VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,exclude=exclude,nfolds=nfolds,alpha=alpha,K=nSS,n_cores=n_cores)
    selection = sharp::SelectedVariables(VariableSelection.outputs)
  }
  
  if(all(!as.logical(selection))){
    newcriterion = critFUN(lm(Ywh ~ NULL))
  }else{
    Xkeep = Xwh[,as.logical(selection)]
    newcriterion = critFUN(lm(Ywh~Xkeep))
  }
  
  if(newcriterion > oldCriterion){
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)),", previous model kept.")
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,paste0("\n              > parameter values : ",
        paste0(c("lambda","thresholds"),"=",c("","0."),
               sapply(sharp::ArgmaxId(VariableSelection.outputs),function(x){
                 10**floor(log10(x))*round(x/(10**floor(log10(x))),digits=2)
               }),collapse=", ")))
  }
  
  model.list = modelFromSelection(Y,X,selection)
  
  
  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat))
}
