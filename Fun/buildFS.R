buildFS <- function(pathToSim,covariateSize,covariateType,
                    temporaryDirectory = NULL,
                    buildMethod = "reg",
                    stabilitySelection = TRUE,
                    cluster=FALSE,
                    p.max = 0.1,
                    weight=NULL,
                    nrep = 8){

  # Simulation and variables
  if(cluster){
    sim = read.csv(pathToSim)

    data <- sim[,c("id","time","yAB")]
    covariate <- sim[sim$time==0,-c(2,3)]

    nCovariate = colnames(covariate[,-1])

    rownames(covariate) <- covariate$id
    covariate <- covariate[,-1]

    covSD0 = names(which(apply(covariate,FUN=sd,MARGIN=2)==0))
    covariate <- covariate[,setdiff(nCovariate,covSD0)]

    treeCovariate <- hclustvar(covariate)
    agregLevel = treeCovariate$height
    agregSmooth = lm(agregLevel ~ poly(cluster,3), data=data.frame(agregLevel=rev(agregLevel),cluster=1:length(agregLevel)))
    cl = as.numeric(names(which(diff(agregSmooth$fitted.values)==max(diff(agregSmooth$fitted.values)))))


    Group = cutreevar(treeCovariate,cl,matsim=TRUE)

    variableCluster=Group$cluster
    variableSynth=cbind(id=rownames(covariate),Group$scores)

    newSim = merge(data,variableSynth,by="id")

    headerTypes = c("id","time","observation",rep("contcov",cl))

    pathToSim = paste0(temporaryDirectory,"/simulation.txt")
    write.csv(newSim,file=pathToSim,quote = FALSE,row.names = FALSE)
  }else{
    load(paste0("Files/",covariateSize,covariateType,"/headerTypes.RData"))
  }


  # Project initialization
  newProject(data = list(dataFile = pathToSim,
                         headerTypes = headerTypes),
             modelFile = "model.txt")

  setIndividualParameterVariability(delta_S=FALSE,delta_L=FALSE)

  setPopulationParameterInformation(
    delta_S_pop=list(initialValue=0.231,method="FIXED"),
    delta_L_pop=list(initialValue=0.000316,method="FIXED"))

  # setCovariateModel(phi_S=list(cAGE=TRUE),phi_L=list(RACE=TRUE),delta_AB=list(SEX=TRUE))

  condSettings = getConditionalDistributionSamplingSettings()

  condSettings$nbsimulatedparameters <- nrep

  setConditionalDistributionSamplingSettings(condSettings)

  # Scenario
  scenario <- getScenario()
  scenario$tasks['standardErrorEstimation'] <- TRUE
  setScenario(scenario)

  saveProject(projectFile=paste0(temporaryDirectory,"/Build.mlxtran"))

  # Model Building
  res = buildmlx(project =paste0(temporaryDirectory,"/Build.mlxtran"),
                 buildMethod = buildMethod,
                 stabilitySelection = stabilitySelection,
                 model="covariate",
                 test=FALSE,
                 p.max = p.max,
                 weight = weight)

  Model <- Rsmlx:::mlx.getIndividualParameterModel()

  if(cluster){
    covariateModel = Model$covariateModel

    empty = rep(FALSE,length(nCovariate))
    names(empty) <- nCovariate
    for(par in 1:length(covariateModel)){
      toSet = empty

      clSelected = stringr::str_remove(names(which(covariateModel[[par]])),"cluster")
      covSelected = names(variableCluster[which(variableCluster %in% clSelected)])
      toSet[covSelected] <- TRUE

      covariateModel[[par]] <- toSet
    }

    Model$covariateModel <- covariateModel
  }
  return(list(Model=Model,time=res$time,iter=res$iter))
}
