buildFS <- function(pathToSim,
                    project,
                    temporaryDirectory = NULL,
                    buildMethod = "stepAIC",
                    FDP_thr=0.1){
  
  load(paste0("data/simulationFiles/Files",project,"/headerTypes.RData"))
  
  # Project initialization
  suppressWarnings({
    newProject(data = list(dataFile = pathToSim,
                           headerTypes = headerTypes),
               modelFile = paste0("data/modelFiles/",if(stringr::str_detect(project,"Pasin")){"Pasin"}else{project},".txt"))
  })
  
  if(stringr::str_detect(project,"Pasin")){
    model = "covariate"
    
    obs.name=getMapping()$mapping[[1]]$model
    eval(parse(text=paste0("setErrorModel(",obs.name,"='constant')")))
    
    setIndividualParameterVariability(delta_S=FALSE,delta_L=FALSE)
    setPopulationParameterInformation(
      delta_S_pop=list(initialValue=0.23,method="FIXED"),
      delta_L_pop=list(initialValue=0.000316,method="FIXED"))
  }else if(project=="Naveau"){
    model = c("covariate","correlation")
    
    obs.name=getMapping()$mapping[[1]]$model
    eval(parse(text=paste0("setErrorModel(",obs.name,"='constant')")))
    
    setIndividualParameterVariability(D=FALSE,V=FALSE)
    setPopulationParameterInformation(D_pop=list(initialValue=100,method="FIXED"),
                                      V_pop=list(initialValue=30,method="FIXED"))
    
    setIndividualParameterDistribution(phi1="normal")
    setIndividualParameterDistribution(phi2="normal")
    
  }
  
  saveProject(projectFile=paste0(temporaryDirectory,"/Build.mlxtran"))
  
  # Model Building
  res = buildmlx(project = paste0(temporaryDirectory,"/Build.mlxtran"),
                 buildMethod = if(buildMethod=="reg"){"stepAIC"}else{"lasso"},
                 model=model,
                 test=FALSE,
                 FDP_thr=FDP_thr)
  
  Model <- Rsmlx:::mlx.getIndividualParameterModel()
  LL <- Rsmlx:::mlx.getEstimatedLogLikelihood()
  return(list(Model=Model,time=res$time,iter=res$iter,LL=LL))
}
