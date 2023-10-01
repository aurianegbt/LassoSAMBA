estimateFS <- function(covariateSize,covariateType,fileNumber,print=TRUE){
  pathToSim=paste0("Files/",covariateSize,covariateType,"/simulation/simulation_",fileNumber,".txt")
  load(paste0("Files/",covariateSize,covariateType,"/headerTypes.RData"))
  
  ## Function core : 
  newProject(data = list(dataFile = pathToSim,
                         headerTypes = headerTypes),
             modelFile = "model.txt")
  
  
  setCovariateModel(phi_S=list(cAGE=TRUE),phi_L=list(RACE=TRUE),delta_AB=list(SEX=TRUE))
  
  setIndividualParameterVariability(delta_S=FALSE,delta_L=FALSE)
  
  setPopulationParameterInformation(
    delta_S_pop=list(initialValue=0.231,method="FIXED"),
    delta_L_pop=list(initialValue=0.000316,method="FIXED"))
  
  contcov <- names(getCovariateInformation()$type[getCovariateInformation()$type=="continuous"])
  
  setErrorModel(yAB_="constant")
  
  scenario <- getScenario()
  scenario$tasks['standardErrorEstimation'] <- TRUE
  setScenario(scenario)
  
  runScenario()
  
  EstimatedPopulationParameters = getEstimatedPopulationParameters()[-c(1,2)] # delta_S, delta_L fixed
  
  EstimatedStandardErrors = getEstimatedStandardErrors()[[1]]
  
  res = data.frame(Parameters = EstimatedPopulationParameters,StandardError = EstimatedStandardErrors$se)
  
  # pour la lisibilité
  aux <- rownames(res)
  aux[aux=="a"] <- "sigma_AB"
  rownames(res) <- aux
  
  cat("______________________________________________________\n")
  return(invisible(res))
}