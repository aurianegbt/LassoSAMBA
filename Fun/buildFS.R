buildFS <- function(pathToSim,covariateSize,covariateType,
                    project="Pasin",
                    temporaryDirectory = NULL,
                    buildMethod = "reg",
                    thresholdsSS = 0.90,
                    weight=NULL,
                    p.max=0.1){

  load(paste0("Files/Files",project,"/",covariateSize,covariateType,"/headerTypes.RData"))


  # Project initialization
  suppressWarnings({
    newProject(data = list(dataFile = pathToSim,
                           headerTypes = headerTypes),
               modelFile = paste0(project,".txt"))
  })

  if(project=="Pasin" | project=="PasinWeird"){
    setIndividualParameterVariability(delta_S=FALSE,delta_L=FALSE)
    
    setPopulationParameterInformation(
      delta_S_pop=list(initialValue=0.23,method="FIXED"),
      delta_L_pop=list(initialValue=0.000316,method="FIXED"))
    
    model="covariate"
    
    
    setErrorModel(yAB_="constant")
  }else if(project=="PK"){
    model="covariate"
  }

  # Scenario
  scenario <- getScenario()
  scenario$tasks['standardErrorEstimation'] <- TRUE
  setScenario(scenario)
  
  saveProject(projectFile=paste0(temporaryDirectory,"/Build.mlxtran"))
  
  
  # Model Building
  res = buildmlx(project = paste0(temporaryDirectory,"/Build.mlxtran"),
                 buildMethod = buildMethod,
                 model=model,
                 weight=if(is.null(weight)){NULL}else{list(covariate=weight)},
                 test=FALSE,
                 thresholdsSS=thresholdsSS,
                 p.max=p.max)

  Model <- Rsmlx:::mlx.getIndividualParameterModel()
  return(list(Model=Model,time=res$time,iter=res$iter))
}
