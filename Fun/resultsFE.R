resultsFE <- function(covariateSize,covariateType,Latex=FALSE){
  load("filesComputed.RData")
  EstimatedStandardDeviation = data.frame()
  EstimatedPopulationParameters = data.frame()
  
  
  for(f in files){
    pathToRes = paste0("Results_Estim/",covariateSize,covariateType,"/estimResults_",f,".RData")
    load(pathToRes)
    
    EstimatedStandardDeviation = rbind(EstimatedStandardDeviation, t(res[,"StandardError",drop=FALSE]))
    EstimatedPopulationParameters = rbind(EstimatedPopulationParameters,t(res[,"Parameters",drop=FALSE]))
    
    rownames(EstimatedPopulationParameters)[length(EstimatedPopulationParameters)] <- rownames(EstimatedStandardDeviation)[length(EstimatedStandardDeviation)] <- paste0("estim_",f)
  }
  
  load("PopulationParameters_theo.RData")
  
  stats_estim = data.frame(target_values = t(PopulationParameters))
  #NaN : 
  rmSD = which(EstimatedStandardDeviation=="NaN",arr.ind = TRUE)[,"col"]
  for(c in rmSD){
    EstimatedStandardDeviation[,c] <- as.numeric(EstimatedStandardDeviation[,c])
  }
  
  for(c in which(EstimatedPopulationParameters=="NaN",arr.ind = TRUE)[,"col"]){
    EstimatedPopulationParameters[,c] <- as.numeric(EstimatedPopulationParameters[,c])
  }
  if(length(rmSD)!=0){
    stats_estim <- cbind(stats_estim,mean=rowMeans(EstimatedPopulationParameters),
                         bias=((rowMeans(EstimatedPopulationParameters)-stats_estim[,'target_values']))/t(PopulationParameters),
                         sd_emp = apply(EstimatedPopulationParameters,1,sd),
                         sd_est = rowMeans(EstimatedStandardDeviation[,-rmSD]))
  }else{
    stats_estim <- cbind(stats_estim,mean=rowMeans(EstimatedPopulationParameters),
                         bias=((rowMeans(EstimatedPopulationParameters)-stats_estim[,'target_values']))/t(PopulationParameters),
                         sd_emp = apply(EstimatedPopulationParameters,1,sd),
                         sd_est = rowMeans(EstimatedStandardDeviation))
  }
  cover = rep(0,10)
  n = ncol(EstimatedPopulationParameters)
  for (j in 1:n){
    sdhat = EstimatedStandardDeviation[,j] 
    lw = EstimatedPopulationParameters[,j] - 1.96*sdhat
    up = EstimatedPopulationParameters[,j] + 1.96*sdhat
    bool = as.numeric(PopulationParameters > lw & PopulationParameters < up)
    cover = cover + bool
  }
  cover <- cover/n 
  stats_estim <- cbind(stats_estim,cover)
  
  
  Power = c(NA,0,NA,0,NA,0,NA,NA,NA,NA)
  ind = c(2,4,6)
  for(j in 1:n){
    sdhat = EstimatedStandardDeviation[ind,j] 
    lw = EstimatedPopulationParameters[ind,j] - 1.96*sdhat
    up = EstimatedPopulationParameters[ind,j] + 1.96*sdhat
    bool = as.numeric(0 > lw & 0 < up)
    for(k in 1:3){
      Power[[ind[k]]] = Power[[ind[k]]] + (1-bool[k])/n
    }
  }
  
  stats_estim <- cbind(stats_estim,Power)
  
  stats_estim[is.na(stats_estim)[,"Power"],"Power"] <- "."
  colnames(stats_estim) <- c("Target Value","Mean of estimation","Relative Bias","Empirical Standard Deviation","Estimated Standard Deviation","Cover","Puissance")
  
  if(Latex){
    stats_estimLTX=stats_estim
    tex2_names = c("$\\varphi_{S,pop}$",
                   "$\\beta_{\\varphi_{S},cAGE}$",
                   "$\\varphi_{L,pop}$",
                   "$\\beta_{\\varphi_{L},RACE_{Eur}}$",
                   "$\\delta_{Ab,pop}$",
                   "$\\beta_{\\delta_{Ab},SEX_M}$",
                   "$\\omega_{\\varphi_{S}}$",
                   "$\\omega_{\\varphi_L}$",
                   "$\\omega_{\\delta_{Ab}}$",
                   "$\\sigma_{Ab}$")
    cat("\n \n ----------- LATEX CODE FOR RESULTS -----------\n \n")
    rownames(stats_estimLTX) <- tex2_names
    colnames(stats_estimLTX) <- c("Target Value","Mean of estimation","Relative Bias","Empirical Standard Deviation","Estimated Standard Deviation","Cover","Puissance")
    print(xtable::xtable(stats_estimLTX,digits=3),sanitize.text.function=function(x){x})
  }
  
  
  return(invisible(stats_estim))
}