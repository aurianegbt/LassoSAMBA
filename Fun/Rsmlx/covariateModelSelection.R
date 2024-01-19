covariateModelSelection <- function(buildMethod,
                                    nfolds = 5,
                                    alpha = 1,
                                    stabilitySelection = TRUE,
                                    nSS=1000,
                                    thresholdsSS=0.90,
                                    thresholdsRep=0.75,
                                    covFix = NULL,
                                    pen.coef=NULL,
                                    weight=1,
                                    n.full = 10,
                                    nb.model = 1,
                                    direction="both",
                                    paramToUse="all",
                                    eta=NULL,
                                    p.max=1,
                                    steps=1000,
                                    sp0=NULL,
                                    iter=1,
                                    printFrequencySS=TRUE,
                                    correlation.model=NULL,
                                    covariate.model=NULL,
                                    criterion = "BIC",
                                    ncrit=20,
                                    replicatesSS=FALSE){
  if(buildMethod=="reg"){
    covariateModelSelection.reg(covFix,pen.coef,weight,n.full,nb.model,direction,paramToUse,eta,p.max,steps,sp0,iter,correlation.model)
  }else if(buildMethod=="StepAIC"){
    covariateModelSelection.StepAIC(covFix,pen.coef,weight,n.full,nb.model,direction,paramToUse,eta,p.max,steps,sp0,iter)
  }else if(buildMethod=="lasso"){
    covariateModelSelection.lasso(nfolds,alpha,stabilitySelection,nSS,thresholdsSS,thresholdsRep,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS,replicatesSS)
  }else if(buildMethod=="rlasso"){
    covariateModelSelection.rlasso(nfolds,alpha,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS)
  }else if(buildMethod=="relasticnet"){
    covariateModelSelection.rlasso(nfolds,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS)
  }else if(buildMethod=="elasticnet"){
    covariateModelSelection.elasticnet(nfolds,stabilitySelection,nSS,thresholdsSS,thresholdsRep,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS,replicatesSS)
  }
}


