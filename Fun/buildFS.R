buildFS <- function(pathToSim,covariateSize,covariateType,
                    project="Pasin",
                    temporaryDirectory = NULL,
                    buildMethod = "reg",
                    stabilitySelection = TRUE,
                    thresholdsSS = 0.90,
                    thresholdsRep= 0.75,
                    cluster=FALSE,
                    weight=NULL,
                    ncrit=20,replicatesSS=FALSE){

  load(paste0("Files/Files",project,"/",covariateSize,covariateType,"/headerTypes.RData"))


  # Project initialization
  newProject(data = list(dataFile = pathToSim,
                         headerTypes = headerTypes),
             modelFile = paste0(project,".txt"))

  if(project=="Pasin" | project=="PasinMore"){
    setIndividualParameterVariability(delta_S=FALSE,delta_L=FALSE)
    
    setPopulationParameterInformation(
      delta_S_pop=list(initialValue=0.23,method="FIXED"),
      delta_L_pop=list(initialValue=0.000316,method="FIXED"))
    
    model="covariate"
  }else if(project=="Warfarine"){
    model=c("covariate","correlation")
  }

  # Scenario
  scenario <- getScenario()
  scenario$tasks['standardErrorEstimation'] <- TRUE
  setScenario(scenario)

  saveProject(projectFile=paste0(temporaryDirectory,"/Build.mlxtran"))
  
  # Model Building
  res = buildmlx(project =paste0(temporaryDirectory,"/Build.mlxtran"),
                 buildMethod = buildMethod,
                 stabilitySelection = stabilitySelection,
                 model=model,
                 test=FALSE,
                 thresholdsSS=thresholdsSS,
                 thresholdsRep=thresholdsRep,
                 ncrit=ncrit,
                 replicatesSS=replicatesSS)

  Model <- Rsmlx:::mlx.getIndividualParameterModel()
  return(list(Model=Model,time=res$time,iter=res$iter))
}
