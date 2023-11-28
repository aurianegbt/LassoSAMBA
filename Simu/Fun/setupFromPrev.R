setupFromPrev <- function(previousProjectNameNumber,previousProjectType,nbcovToAdd,randomCov = TRUE,covSettingsToAdd=NULL,corcovSettingsToAdd=NULL,covariateChar="X",print=TRUE,individualsNumber=100,NORM=TRUE){
  ## Variable needed
  source("Fun/genCov.R")
  source("Fun/addCov.R")
  source("Fun/randomCovariate.R")
  catDist = c("binom","sample","binary","binomial")
  contDist = setdiff(c("beta","binom","cauchy","chisq","exp","f","gamma","geom","norm","pois","sample","unif","beta","exponential","gamma","normal","noZeroPoisson","poisson","uniform"),catDist)

  ## definition and files initialization
  newProjectNameNumber = previousProjectNameNumber + nbcovToAdd
  ProjectName = paste0(newProjectNameNumber,previousProjectType)

  SubProjectName = setdiff(stringr::str_sub(list.dirs(paste0("Files/",previousProjectNameNumber,previousProjectType)),start=(7+nchar(paste0(previousProjectNameNumber,previousProjectType))+1)),c("","covTable"))

  Folder = paste0("Files/",ProjectName)
  if(dir.exists(Folder)){unlink(Folder, recursive=TRUE)}
  dir.create(Folder)
  dir.create(paste0(Folder,"/covTable"))
  for(p in paste0(Folder,"/",SubProjectName)){dir.create(p)}
  rm(p)

  ## Cov Settings
  if(!randomCov){
    if(is.null(covSettingsToAdd) && is.null(corcovSettingsToAdd)){
      stop("Please provide a settings for covariate if random definition is disabled.")
    }else if(length(covSettingsToAdd)+length(corcovSettingsToAdd)!=nbcovToAdd){
      warning(paste0("The length of covSettings to add, ",length(covSettingsToAdd),", and the number of covariates to add, ",nbcovToAdd," does not match."))
  }}else{
    covSettingsToAdd = randomCovariate(n=nbcovToAdd,NORM=NORM)
    names(covSettingsToAdd) <- paste0(covariateChar,1:nbcovToAdd)
  }
  covSettings <- covSettingsToAdd
  names(covSettings) <- names(covSettingsToAdd)
  corcovSettings <- corcovSettingsToAdd
  names(corcovSettings) <- names(corcovSettingsToAdd)
  if(!is.null(corcovSettingsToAdd)){
    save(covSettings,corcovSettings,file=paste0(Folder,"/covSettingsAdded.RData"))
  }else{save(covSettings,file=paste0(Folder,"/covSettingsAdded.RData"))}

  distribution = unlist(lapply(c(covSettingsToAdd,corcovSettingsToAdd),function(x){x$distribution}),use.names = FALSE)


  # for sim.setup update
  headerAdded = distribution %in% catDist
  headerAdded[headerAdded] <- "catcov"
  headerAdded[headerAdded==FALSE] <- "contcov"

  load(paste0("Files/",previousProjectNameNumber,previousProjectType,"/headerTypes.RData"))
  headerTypes <- c(headerTypes,headerAdded)
  save(headerTypes,file=paste0("Files/",ProjectName,"/headerTypes.RData"))

  ## Files creation
  for(i in 1:100){
    covNew = genCorCov(covSettings,corcovSettings,individualsNumber,oldCov=read.csv(paste0("Files/",previousProjectNameNumber,previousProjectType,"/covTable/covTable_",i,".txt")))
    pathToCov =  paste0("Files/",previousProjectNameNumber,previousProjectType,"/covTable/covTable_",i,".txt")
    pathToSaveCov = paste0("Files/",ProjectName,"/covTable/covTable_",i,".txt")

    pathToSimS = paste0("Files/",previousProjectNameNumber,previousProjectType,"/",SubProjectName,"/simulation_",i,".txt")
    pathToSaveSimS = paste0("Files/",ProjectName,"/",SubProjectName,"/simulation_",i,".txt")

    addCov(pathToCov,pathToSaveCov,covNew,pathToSimS,pathToSaveSimS,print=print)
  }
}
