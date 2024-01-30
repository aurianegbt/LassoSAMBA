applyMethodsharp <- function(Y,X,omega,cov0,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             p.name=NULL){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètr
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
    selection = sharp::SelectedVariables(sharp::VariableSelection(Xwh,Ywh,exclude=exclude,nfolds=nfolds,alpha=alpha,K=nSS))
  }
  
  
  model.list = modelFromSelection(Ywh,Xwh,selection)
  
  
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name))
}
