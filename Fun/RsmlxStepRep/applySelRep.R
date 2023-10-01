applySelRep <- function(Y,X,Sigma=NULL,alpha=1,cov0.list=NULL,
                              stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.90,
                              parallel=FALSE,ncores=1,ncrit=20,fitList=NULL,critMV="BIC"){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre 
  if(!is.null(fitList)){
    oldCriterion = mvIC(fitList = fitList,criterion = critMV)
    print(paste0("old Criterion ",critMV," : ",round(oldCriterion,digits=2)))
    
    prevSel = #TODO 
  }
  selection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS,parallel,ncores)
  model.list = modelFromLassoSelection(Y,X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)
  
  if(!is.null(fitList)){
    newCriterion = mvIC(model.list,criterion=critMV)
    tested = 1
    print(paste0("new Criterion ",tested, " : ",round(newCriterion,digits=2)))
    flag = newCriterion < oldCriterion
    
    while(tested < ncrit & !flag){
      selection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS,parallel,ncores)
      model.list = modelFromLassoSelection(Y,X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)
      tested = tested + 1
      newCriterion = mvIC(model.list,criterion=critMV)
      print(paste0("new Criterion ",critMV," ",tested, " : ",round(newCriterion,digits=2)))
      flag = newCriterion < oldCriterion
    }
    if(!flag){
      selection = prevSel
    }
  }
  return(selection) # prevSel
}