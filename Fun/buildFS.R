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


    pathToSimOriginal = pathToSim
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
    cat("\nCurrent clustered model : \n")
    print(covariateModel)
    cat("Selection whithin each cluster...")


    pos = c(1:length(nCovariate))
    names(pos) <- nCovariate

    cov0.list = list()


    Rsmlx:::mlx.runConditionalDistributionSampling()
    sp.df <- Rsmlx:::mlx.getSimulatedIndividualParameters()
    nrep <- max(sp.df$rep)
    ind.dist <- Rsmlx:::mlx.getIndividualParameterModel()$distribution
    param.names <- names(ind.dist)
    n.param <- length(param.names)
    cov.info <- Rsmlx:::mlx.getCovariateInformation()
    cov.names <- cov.info$name
    cov.types <- cov.info$type
    tcov.names <- NULL
    covariates <- cov.info$covariate
    cov.cat <- cov.names[cov.types == "categorical"]
    covariates[cov.cat] <- lapply(covariates[cov.cat], as.factor)
    indvar <- Rsmlx:::mlx.getIndividualParameterModel()$variability$id
    cov.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
    r <- res <- r.cov0 <- list()
    eps <- 1e-15

    # Y transformation
    Y=sp.df[,c("rep","id",names(indvar)[which(indvar==TRUE)])]
    for(j in 1:n.param){
      dj <- ind.dist[j]
      nj <- names(dj)
      if (indvar[j]) {
        if (tolower(dj) == "lognormal") {
          Y[,nj] <- log(Y[,nj] + eps)
          colnames(Y)[which(colnames(Y)==nj)] <- paste0("log.", nj)
        }else if (tolower(dj) == "logitnormal") {
          Y[,nj] <- log((Y[,nj] + eps)/(1 - Y[,nj] + eps))
          colnames(Y)[which(colnames(Y)==nj)] <- paste0("logit.", nj)
        }else if (tolower(dj) == "probitnormal") {
          Y[,nj] <- qnorm(Y[,nj])
          colnames(Y)[which(colnames(Y)==nj)]  <- paste0("probit.", nj)
        }
      }
    }
    # param.names for order in r :

    Sigma=Rsmlx::getEstimatedCovarianceMatrix()$cov.matrix
    rownames(Sigma) <- colnames(Sigma) <- rownames(Rsmlx::getEstimatedCovarianceMatrix()$cor.matrix)


    for(par in names(indvar)[which(indvar==TRUE)]){
      clSelected = stringr::str_remove(names(which(covariateModel[[par]])),"cluster")
      covSelected = names(variableCluster[which(variableCluster %in% clSelected)])

      cov0.list <- append(cov0.list,list(setdiff(nCovariate,covSelected)))
      names(cov0.list)[length(cov0.list)] <- par
    }

    N = length(unique(Y$id))
    Y.mat = sapply(Y[,-c(1,2)],function(x){rowMeans(matrix(x,nrow=N))}) #1 : rep 2 : id
    covariates = read.csv(pathToSimOriginal) %>% filter(time==0)
    covariates <- covariates[,nCovariate]
    X.mat  = covariates

    selection = lassoSelection(Y.mat,X.mat,Sigma,cov0.list=cov0.list)

    empty = rep(FALSE,length(nCovariate))
    names(empty) <- nCovariate
    for(par in names(indvar)[which(indvar==TRUE)]){
      toSet = empty
      covSelected = colnames(selection)[which(as.logical(selection[par,]))]
      toSet[covSelected] <- TRUE

      covariateModel[[par]] <- toSet
    }

    Model$covariateModel <- covariateModel

    cat("\nFinal Covariate model:\n")
    print(covariateModel)

  }
  return(list(Model=Model,time=res$time,iter=res$iter))
}
