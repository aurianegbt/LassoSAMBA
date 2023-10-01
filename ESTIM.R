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
invisible(purrr::quietly(initializeLixoftConnectors)(software='monolix', force=T))
source("Fun/estimateFS.R")


# Arguments of batch 
arr <- as.numeric(Slurm_env(x='SLURM_ARRAY_TASK_ID'))

## Presentation of work : 
cat("===========================================================================\n       - - - - - - - - - - - Work Information - - - - - - - - - - -       ")
cat("Method used : SAMBA")
cat(paste0("Working on file n°",arr))
cat("===========================================================================\n")



# Folder and paths 
pathToSim = paste0("Files/3cov//simulation_",arr,".txt")
pathToCov = paste0("Files/3cov/covTable/covTable_",arr,".txt")
pathToResults = paste0("Results/Results_Estim")
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToResults = paste0(pathToResults,"/3cov")
if(!dir.exists(pathToResults)){dir.create(pathToResults)}
pathToResults=paste0(pathToResults,"/estimResults_",arr,".RData")

# Work ! 
if(file.exists(pathToResults)){print("Estimate already computed.\n")}else{
  res = estimateFS(3,"cov",arr)
  
  rownames(res)=c("phi_S_pop","beta_cAGE","phi_L_pop","beta_RACE_E","delta_AB_pop","beta_SEX_M","omega_S","omega_L","omega_AB","sigma_AB")
  
  save(res,file=pathToResults)
}
cat("===========================================================================\n")


