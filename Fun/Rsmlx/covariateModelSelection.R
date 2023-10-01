covariateModelSelection <- function(buildMethod,
                                    nfolds = 5,
                                    alpha = 1,
                                    parallel = FALSE,
                                    ncores=1,
                                    stabilitySelection = TRUE, 
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
                                    correlation.model=NULL){
  if(buildMethod=="reg"){
    covariateModelSelection.reg(covFix,pen.coef,weight,n.full,nb.model,direction,paramToUse,eta,p.max,steps,sp0,iter,correlation.model)
  }else if(buildMethod=="lasso"){
    covariateModelSelection.lasso(nfolds,alpha,parallel,ncores,stabilitySelection,nSS,thresholdsSS,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0)
  }
}


