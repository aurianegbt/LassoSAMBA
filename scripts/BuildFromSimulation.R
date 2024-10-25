# empty the memory
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Load Functions and softwares
ibrary(lixoftConnectors)
library(foreach)
initializeLixoftConnectors(software = "monolix")

## Fun summarize and call
source("scripts/modelBuildingFun/buildFS.R")

# Function just modified
source("scripts/modelBuildingFun/buildmlx.R")
source("scripts/modelBuildingFun/covariateModelSelection.R")
source("scripts/modelBuildingFun/covariateModelSelection.reg.R")

# Function implemented
source("scripts/modelBuildingFun/modelFromSelection.R")
source("scripts/modelBuildingFun/updateCov0.R")
source("scripts/modelBuildingFun/covariateModelSelection.lasso.R")
source("scripts/modelBuildingFun/applyMethodLasso.R")



# Arguments of batch
# For sbatch use
arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID')) # number of replicates 

args = commandArgs(trailingOnly=TRUE)
# project,  method 
# args=c("GaussianPasin","lasso")
project = args[1]     # "Pasin", "GaussianPasin" or "Naveau" 
buildMethod = args[2] # "lasso" or "reg"
FDP_thr = as.numeric(args[3])

## Presentation of work :
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat(" • Method used : ",buildMethod,";\n")
cat("Working on file n°",arr," with",if(project=="Naveau"){"500"}else{"200 correlated"}," covariates covariates, on project ",project,".")
if(buildMethod=="lasso"){cat("\n Model building done without statistical tests.")}
cat("\n===========================================================================\n")


# Folder and paths
dir("outputs/buildingResults/")
dir("outputs/buildingResults/simulation")
pathToResults = paste0("outputs/buildingResults/simulation/Results",project)
dir(pathToResults)
pathToResults = paste0(pathToResults,"/Results_",buildMethod,if(buildMethod=="lasso"){paste0("FDP",FDP_thr*100)})
dir(pathToResults)

pathToFolderResults = pathToResults
pathToCompResults=paste0(pathToResults,"/compResults_",arr,".RData")
pathToResults=paste0(pathToResults,"/buildResults_",arr,".RData")

# Temporary Files :
temporaryDirectory = paste0("tmp",arr,project,buildMethod,paste0("FDP",FDP_thr*100))
if(dir.exists(temporaryDirectory)){
  warning("A folder already exists for this job. It has been deleted.")
  unlink(temporaryDirectory,recursive=T,force=T)
}
dir.create(temporaryDirectory)


pathToSim = paste0("data/simulationFiles/Files",project,"/simulation/simulation_",arr,".txt")

if(file.exists(pathToResults) && file.exists(pathToCompResults)){cat("Model already built.\n")}else{
  
  resBuild = buildFS(pathToSim,
                     project,
                     temporaryDirectory = temporaryDirectory,
                     buildMethod = buildMethod,
                     FDP_thr=FDP_thr)
  
  model = resBuild$Model
  time = resBuild$time
  iter = resBuild$iter
  LL = resBuild$LL
  
  save(model,LL,file=pathToResults)
  save(time,iter,file=pathToCompResults)
}
unlink(temporaryDirectory,recursive = TRUE)
cat("\n===========================================================================\n")
