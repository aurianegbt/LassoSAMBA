applyMethodLasso <- function(Y,X,Sigma=NULL,alpha=1,cov0.list=NULL,
                              stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.90,
                              parallel=FALSE,ncores=1){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre 
  
  selection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS,parallel,ncores)
  model.list = modelFromLassoSelection(Y,X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)
  
  n.param = ncol(Y)
  param.names = names(cov0.list)
  tparam.names = colnames(Y)
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