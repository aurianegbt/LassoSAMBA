applyMethodLasso <- function(Y,X,Sigma=NULL,alpha=1,cov0.list=NULL,
                              stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.90,
                              ncrit=20,covariate.model=NULL,critMV="BIC"){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre
  param.names = names(cov0.list)
  tparam.names = colnames(Y)
  n.param = ncol(Y)


  fitData = cbind(Y,X)
  fit.list = list()
  for(j in 1:n.param){
    nj = param.names[j]
    tnj = tparam.names[j]
    covaAux = names(covariate.model[[nj]][covariate.model[[nj]]])
    if(length(covaAux)!=0){
      formula = paste0(tnj," ~ 1 + ",paste0(covaAux,collapse=" + "))
    }else{
      formula = paste0(tnj," ~ 1")
    }

    fitModel = eval(parse(text=paste0("lm(",formula,",data=fitData)")))

    fit.list <- append(fit.list,list(fitModel))
  }


  oldCriterion = mvIC(fitList = fit.list,criterion = critMV)
  print(paste0("old Criterion ",critMV," : ",round(oldCriterion,digits=2)))

  selection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS)
  model.list = modelFromLassoSelection(Y,X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS)

  newCriterion = mvIC(model.list,criterion="BIC")
  tested = 1
  print(paste0("new Criterion ",tested, " : ",round(newCriterion,digits=2)))
  flag = newCriterion < oldCriterion

  while(tested < ncrit & !flag){
    selection = lassoSelection(Y,X,Sigma,alpha,cov0.list,stabilitySelection,nfolds,nSS,thresholdsSS)
    model.list = modelFromLassoSelection(Y,X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS)
    tested = tested + 1
    newCriterion = mvIC(model.list,criterion="BIC")
    print(paste0("new Criterion ",tested, " : ",round(newCriterion,digits=2)))
    flag = newCriterion < oldCriterion
  }


  if(!flag){
    model.list=fit.list
  }

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
