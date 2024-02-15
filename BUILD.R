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
suppressWarnings(library(foreach))
suppressWarnings(suppressMessages(library(lixoftConnectors)))

# Load Functions and softwares
invisible(purrr::quietly(lixoftConnectors::initializeLixoftConnectors)(software='monolix', force=T))

# Arguments of batch
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)
# args = c("Pasin","sharp","cov","50","NULL","NULL","FALSE")
# project,  method,type, size, weight, thresholdsSS , nocov0
# seq(0.6,0.95,length.out=50)

project = args[1]
buildMethod <- args[2]
covariateType <- args[3]
covariateSize <- as.numeric(args[4])
weight <- if(args[5]=="NULL"){NULL}else{as.numeric(args[5])}
  PEN <- !is.null(weight)
eval(parse(text=paste0("thresholdsSS = ",args[6])))
noCov0 = as.logical(args[7])
source("Fun/Rsmlx/fun.R")

## Presentation of work :
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat(" • Method used : ",
    if(buildMethod=="reg"){
      "stepAIC"
      }else if(buildMethod=="lasso"){
        "lasso"
      }else if(buildMethod=="lassoSS"){
        "lasso with stability selection"
      }else if(buildMethod=="elasticnet"){
        "elastic net"
      }else if(buildMethod=="elasticnetSS"){
        "elastic net with stability selection"
      }else if(buildMethod=="rlasso"){
        "lasso with stability selection on replicates"
      }else if(buildMethod=="relasticnet"){
        "elastic net with stability selection on replicates"
      }else if(buildMethod=="sharp"){
        "sharp"
      }else if(buildMethod=="rsharp"){
        "sharp on replicates"
      }else{buildMethod}," ",
    if((buildMethod %in% c("lassoSS","elasticnetSS","rlasso","relasticnet"))
       & length(thresholdsSS)!=1){
      "using multiple thresholds"
      },";\n")
cat(" • Penalization : ",if(is.null(weight)){"None"}else{weight},";\n")
cat("Working on file n°",arr," from ",project," simulations with ",covariateSize,if(covariateType=="corcov"){" correlated"}," covariates.\n")
if(noCov0){cat("\n Model building done without statistical tests.")}
cat("===========================================================================\n")


# Folder and paths
pathToResults = paste0("Results/Results",project)
  dir(pathToResults)
pathToResults = paste0(pathToResults,"/Results_",buildMethod,
                       if((buildMethod %in% c("lassoSS","elasticnetSS","rlasso","relasticnet"))
                          & length(thresholdsSS)!=1){"Crit"},
                       if(PEN){paste0("PEN",paste0(weight,collapse="-"))},
                       if(noCov0){"noCov0"})
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


if(file.exists(pathToResults) && file.exists(pathToCompResults)){cat("Model already built.\n")}else{
  
  resBuild = buildFS(pathToSim,covariateSize,covariateType,
                     project,
                     temporaryDirectory = temporaryDirectory,
                     buildMethod = buildMethod,
                     thresholdsSS = thresholdsSS,
                     weight = weight,
                     p.max=if(noCov0){1}else{0.1})

  model = resBuild$Model
  time = resBuild$time
  iter = resBuild$iter

  save(model,file=pathToResults)
  save(time,iter,file=pathToCompResults)
}
unlink(temporaryDirectory,recursive = TRUE)
cat("\n===========================================================================\n")
