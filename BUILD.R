# empty the memory
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Required Packages :
# suppressMessages(library(Rsmlx))
# suppressMessages(library(stats))
# suppressMessages(library(MASS))
# suppressMessages(library(glmnet))
# suppressMessages(library(stringr))
# suppressMessages(library(dplyr))
# suppressMessages(library(slurmR))
suppressMessages(library(lixoftConnectors))

# Load Functions and softwares
invisible(purrr::quietly(lixoftConnectors::initializeLixoftConnectors)(software='monolix', force=T))

# Arguments of batch
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)
# args=c("Warfarine","reg","cov","10","FALSE","NULL","0.7")
# project,  method,type, size, rep, weight, thresholds 

project = args[1]
buildMethod <- args[2]
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
  ncrit=20
}    

covariateType <- args[3]
covariateSize <- as.numeric(args[4])


REP = as.logical(args[5])
weight <- if(args[6]=="NULL"){NULL}else{list(covariate=as.numeric(args[6]))}
PEN <- !is.null(weight)
eval(parse(text=paste0("thresholdsSS = ",args[7])))

source("Fun/Rsmlx/fun.R")
# source("Fun/buildFS.R")

## Presentation of work :
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat(" • Method used : ",
    if(buildMethod=="reg"){"stepAIC"}else if(buildMethod=="lasso"){"lasso"}else
      if(buildMethod=="elasticnet"){"elastic net"}," ",
    if(stabilitySelection){"with stability selection"}," ",
    if(length(thresholdsSS)!=1){"using multiple thresholds"},";\n")
cat(" • Penalization : ",if(is.null(weight)){"None"}else{weight},";\n")
cat(" • StabSel on replicates : ",REP,";\n")
cat("Working on file n°",arr," from ",project," simulations with ",covariateSize,if(covariateType=="corcov"){" correlated"}," covariates.\n")
cat("===========================================================================\n")


# Folder and paths
pathToResults = paste0("Results/Results",project)
  dir(pathToResults)
pathToResults = paste0(pathToResults,"/Results_",buildMethod,
                       if(stabilitySelection){"SS"},
                       if(buildMethod!="reg" & length(thresholdsSS)!=1){"Crit"},
                       if(PEN){paste0("PEN",paste0(weight,collapse="-"))},
                       if(REP){"REP"})
  dir(pathToResults)
  
pathToResults <- paste0(pathToResults,"/",covariateSize,covariateType)
  dir(pathToResults)
  
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
replicatesSS=REP


if(file.exists(pathToResults) && file.exists(pathToCompResults)){cat("Model already built.\n")}else{


  resBuild = buildFS(pathToSim,covariateSize,covariateType,
                     project,
                     temporaryDirectory = temporaryDirectory,
                     buildMethod = buildMethod,
                     stabilitySelection = stabilitySelection,
                     thresholdsSS = thresholdsSS,
                     cluster=cluster,
                     weight= weight,
                     ncrit=ncrit,
                     replicatesSS=REP)

  model = resBuild$Model
  time = resBuild$time
  iter = resBuild$iter

  save(model,file=pathToResults)
  save(time,iter,file=pathToCompResults)
}
unlink(temporaryDirectory,recursive = TRUE)
cat("\n===========================================================================\n")
