applyMethodLasso <- function(Y.list,X,Sigma=NULL,alpha=1,cov0.list=NULL,
                        stabilitySelection=FALSE,nfolds=5,nSS=1000,thresholdsSS=0.90,
                        parallel=FALSE,ncores=1,ncrit=20,covariate.model=NULL,critMV="BIC"){

  nrep = length(Y.list)
  param.names = names(cov0.list)
  tparam.names = colnames(Y.list[[1]])
  n.param = ncol(Y.list[[1]])

  selectionSummary = matrix(data = 0,nrow=n.param,ncol=ncol(X))
  colnames(selectionSummary) <- colnames(X)
  rownames(selectionSummary) <- tparam.names

  for(k in 1:nrep){
    fitData = cbind(Y.list[[k]],X)
    fit.list = list()
    oldSelection = data.frame()
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

      oldSelection <- rbind(oldSelection,as.numeric(covariate.model[[nj]][colnames(X)]))
      rownames(oldSelection)[j] <- tnj
    }
    colnames(oldSelection) <- colnames(X)

    oldCriterion = mvIC(fitList = fit.list,criterion = critMV)

    selection = lassoSelection(Y.list[[k]],X,Sigma,alpha,cov0.list,stabilitySelection=FALSE,nfolds,nSS,thresholdsSS,parallel,ncores)
    model.list = modelFromLassoSelection(Y.list[[k]],X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = FALSE,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)


    newCriterion = mvIC(model.list,criterion="BIC")
    tested = 1
    flag = newCriterion < oldCriterion

    while(tested < ncrit & !flag){
      selection = lassoSelection(Y.list[[k]],X,Sigma,alpha,cov0.list,stabilitySelection=FALSE,nfolds,nSS,thresholdsSS,parallel,ncores)
      model.list = modelFromLassoSelection(Y.list[[k]],X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = FALSE,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)
      tested = tested + 1
      newCriterion = mvIC(model.list,criterion="BIC")
      flag = newCriterion < oldCriterion
    }

    if(!flag){
      model.list=fit.list
      selection = oldSelection
    }

    selectionSummary =  selectionSummary + selection

  }
  selectionSummary <- selectionSummary/nrep

  selection = selectionSummary
  selection[selection>=thresholdsSS] = 1
  selection[selection<thresholdsSS] = 0

  model.list = modelFromLassoSelection(do.call(rbind, Y.list),X,selection=selection,Sigma = Sigma,alpha = alpha,cov0.list = cov0.list,stabilitySelection = stabilitySelection,nfolds = nfolds,nSS = nSS,thresholdsSS = thresholdsSS,parallel = parallel,ncores = ncores)


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
