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
suppressMessages(library(ClustOfVar))

# Load Functions and softwares
invisible(purrr::quietly(initializeLixoftConnectors)(software='monolix', force=T))

# Arguments of batch
arr <- as.numeric(Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)

covariateType <- args[1]
buildMethod <- args[2]
covariateSize <- as.numeric(args[3])
Rsmlx <- args[4]
cluster <- as.logical(args[5])
buildOptions <- args[6]
weight <- args[7]
nrep <- as.integer(args[8])
if(weight=="NULL"){
  weight <- NULL
}

if(Rsmlx=="RsmlxStep" || Rsmlx=="RsmlxStepRep"){
 suppressWarnings(library(mvIC))
}

source(paste0("Fun/",Rsmlx,"/fun.R"))


# Arguments handling
if(buildMethod == "lassoSS"){
  buildMethod="lasso"
  stabilitySelection = TRUE
}else{
  stabilitySelection = FALSE
}
if(buildOptions=="Mstar"){
  source("Fun/buildFStar.R")
}else{
  source("Fun/buildFS.R")
}

if(buildOptions=="0"){
  p.max = 1
}else{
  p.max = 0.1
}
if(!is.null(weight)){
  Rsmlx <- paste0(Rsmlx,"PEN",weight)
  weight = list(covariate=as.numeric(weight))
}
if(nrep != 8){
  Rsmlx <- paste0(Rsmlx,"REP",nrep)
}



for(arr in 1:100){
  ## Presentation of work :
  cat("===========================================================================\n
    - - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
  cat(" • Method used : ")
  if(buildMethod=="reg"){cat("SAMBA")}else if(buildMethod=="lassoSS"){cat("lasso with stability selection within SAMBA")}else(cat(buildMethod))
  cat("\n")
  cat(" • Version : ")
  if(Rsmlx=="Rsmlx"){cat("v.1 SAMBA standard")}else if(Rsmlx=="Rsmlx2"){cat("v.2 SAMBA with stop test on overall model")}else if(Rsmlx=="RsmlxStep"){cat("v.3 SAMBA with IC test on lasso model")}else(cat(Rsmlx))
  cat("\n")
  cat(" • Clusterisation done : ")
  cat(cluster)
  cat("\n")
  cat(" • Options : ")
  if(buildOptions==""){cat("None")}else if(buildOptions=="Mstar"){cat("starting from M*")}else if(buildOptions=="0"){cat(" without statistical test to exclude covariates from search")}else(cat(buildOptions))
  cat("\n")
  cat(" • Penalization : ")
  if(is.null(weight)){cat(" None")}else{cat(paste0(" ",weight))}
  cat("\n")
  cat(" • Replicas : ")
  cat(nrep)
  cat("\n\n")
  cat(paste0("Working on file n°",arr," with ",covariateSize))
  if(covariateType=="corcov"){cat(" correlated")}
  cat(" covariates.\n")
  cat("===========================================================================\n")
  # Folder and paths
  pathToResults = paste0("Results",buildOptions,"/Results",stringr::str_remove(Rsmlx,"Rsmlx"),"/")
  if(!dir.exists(pathToResults)){dir.create(pathToResults)}
  pathToResults = paste0(pathToResults,"Results_",buildMethod)
  if(stabilitySelection){pathToResults=paste0(pathToResults,"SS")}
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

  pathToSim = paste0("Files/",covariateSize,covariateType,"/simulation/simulation_",arr,".txt")


  # Work !
  if(file.exists(pathToResults) && file.exists(pathToCompResults)){cat("Model already built.\n")}else{


    resBuild = buildFS(pathToSim,covariateSize,covariateType,
                       temporaryDirectory = temporaryDirectory,
                       buildMethod,
                       stabilitySelection = stabilitySelection,
                       cluster=cluster,
                       p.max = p.max,
                       weight= weight,
                       nrep = nrep)

    model = resBuild$Model
    time = resBuild$time
    iter = resBuild$iter

    save(model,file=pathToResults)
    save(time,iter,file=pathToCompResults)
  }
  unlink(temporaryDirectory,recursive = TRUE)
  cat("\n===========================================================================\n")
}
