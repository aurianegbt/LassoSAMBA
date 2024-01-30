covariateModelSelection <- function(buildMethod,
                                    nfolds = 5,
                                    alpha = 1,
                                    nSS=1000,
                                    thresholdsSS=0.90,
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
                                    ncrit=20){
  if(buildMethod %in% c("reg","stepAIC")){
    covariateModelSelection.reg(covFix,pen.coef,weight,n.full,nb.model,direction,paramToUse,eta,p.max,steps,sp0,iter,correlation.model)
  }else if(buildMethod=="lasso"){
    covariateModelSelection.lasso(nfolds,alpha,stabilitySelection=FALSE,nSS,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS)
  }else if(buildMethod=="lassoSS"){
    covariateModelSelection.lasso(nfolds,alpha,stabilitySelection=TRUE,nSS,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS=FALSE)
  }else if(buildMethod=="elasticnet"){
    covariateModelSelection.elasticnet(nfolds,stabilitySelection=FALSE,nSS,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS)
  }else if(buildMethod=="elasticnetSS"){
    covariateModelSelection.elasticnet(nfolds,stabilitySelection=TRUE,nSS,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit,printFrequencySS=FALSE)
  }else if(buildMethod=="rlasso"){
    covariateModelSelection.rlasso(nfolds,alpha,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit)
  }else if(buildMethod=="relasticnet"){
    covariateModelSelection.relasticnet(nfolds,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,covariate.model,criterion,ncrit)
  }else if(buildMethod=="rsharp"){
    covariateModelSelection.rsharp(nfolds,alpha,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0)
  }else if(buildMethod=="sharp"){
    covariateModelSelection.sharp(nfolds,alpha,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,nSS)
  }
}


