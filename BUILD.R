# empty the memory
rm(list=ls())

# Required Packages :
suppressMessages(library(Rsmlx))
suppressMessages(library(lixoftConnectors))
suppressMessages(library(stats))
suppressMessages(library(MASS))
suppressMessages(library(glmnet))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(slurmR))

# Load Functions and softwares
invisible(purrr::quietly(lixoftConnectors::initializeLixoftConnectors)(software='monolix', force=T))

# Arguments of batch
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)

project = args[1]

covariateType <- args[2]

buildMethod <- args[3]
if(buildMethod == "lassoSS"){
  buildMethod="lasso"
  stabilitySelection = TRUE
  ncrit=20
}else if(buildMethod=="elasticnetSS"){
  buildMethod="elasticnet"
  stabilitySelection = TRUE
  ncrit=20
}else{
  stabilitySelection = FALSE
  ncrit=100
}    

covariateSize <- as.numeric(args[4])

cluster <- as.logical(args[5])

weight <- args[6]
    if(weight=="NULL"){weight <- NULL}else{
      Rsmlx <- paste0(Rsmlx,"PEN",weight);weight = list(covariate=as.numeric(weight))}

eval(parse(text=paste0("thresholdsSS = ",args[7])))


source("Fun/Rsmlx/fun.R")
source("Fun/buildFS.R")

## Presentation of work :
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat(" • Method used : ")
if(buildMethod=="reg"){cat("stepAIC")}else if(buildMethod=="lasso"){cat("lasso")}else if(buildMethod=="elasticnet"){cat("elastic net")}else{cat(buildMethod)}
if(stabilitySelection){cat(" with stability selection")}
if(length(thresholdsSS)!=1){cat(" using multiple thresholds")}
cat(";\n")
cat(" • Clusterisation : ")
cat(cluster)
cat(";\n")
cat(" • Penalization : ")
if(is.null(weight)){cat(" None")}else{cat(paste0(" ",weight))}
cat(";\n")
cat(paste0("Working on file n°",arr," from ",project," simulations with ",covariateSize))
if(covariateType=="corcov"){cat(" correlated")}
cat(" covariates.\n")
cat("===========================================================================\n")


# Folder and paths
pathToResults = paste0("Results/Results",project)
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToResults = paste0(pathToResults,"/")
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToResults = paste0(pathToResults,"Results_",buildMethod)
if(stabilitySelection){pathToResults=paste0(pathToResults,"SS")}
if(buildMethod!="reg" & length(thresholdsSS)!=1){pathToResults=paste0(pathToResults,"Crit")}
if(cluster){pathToResults=paste0(pathToResults,"Cl")}
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToResults <- paste0(pathToResults,"/",covariateSize,covariateType)
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToFolderResults = pathToResults
pathToCompResults=paste0(pathToResults,"/compResults_",arr,".RData")
pathToResults=paste0(pathToResults,"/buildResults_",arr,".RData")

# Temporary Files :
temporaryDirectory = paste0("tmp",covariateSize,arr,covariateType)
if(dir.exists(temporaryDirectory)){
  warning("A folder already exists for this job. It has been deleted.")
  unlink(temporaryDirectory,recursive=T,force=T)
}
dir.create(temporaryDirectory)

pathToSim = paste0("Files/Files",project,"/",covariateSize,covariateType,"/simulation/simulation_",arr,".txt")



if(file.exists(pathToResults) && file.exists(pathToCompResults)){cat("Model already built.\n")}else{


  resBuild = buildFS(pathToSim,covariateSize,covariateType,
                     project,
                     temporaryDirectory = temporaryDirectory,
                     buildMethod = buildMethod,
                     stabilitySelection = stabilitySelection,
                     thresholdsSS = thresholdsSS,
                     cluster=cluster,
                     weight= weight,
                     ncrit=ncrit)

  model = resBuild$Model
  time = resBuild$time
  iter = resBuild$iter

  save(model,file=pathToResults)
  save(time,iter,file=pathToCompResults)
}
unlink(temporaryDirectory,recursive = TRUE)
cat("\n===========================================================================\n")
