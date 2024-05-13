# Empty memory ------------------------------------------------------------
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Load packages, arguments, functions and  data  --------------------------
suppressWarnings({
  suppressMessages({
    library(lixoftConnectors)
    library(foreach)
    library(dplyr)
    initializeLixoftConnectors(software = "monolix",
                               path = "/cm/shared/dev/modules/generic/apps/tools/monolix/2023R1bis/",
                               force = TRUE)
  })
})

nbBatch <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)

buildMethod = args[1]
eval(parse(text=paste0("thresholdsSS = ",args[2])))
noCov0 = as.logical(args[3])
eval(parse(text=paste0("model = ",args[4])))
nbCov = as.numeric(args[5])

source("scripts/MBFun.R")


# Presentation of Job -----------------------------------------------------
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat("Application to Prevac-UP data, batched version\n")
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
    }else{buildMethod}," ",
    if((buildMethod %in% c("lassoSS","elasticnetSS","rlasso","relasticnet"))
       & length(thresholdsSS)!=1){
      "using multiple thresholds"
    },"\n")
cat(" • Model(s) built : ",paste0(model,collapse=", "))
if(!is.null(thresholdsSS))
  cat("\n • thresholds used :",paste0(thresholdsSS,collapse=", "))
cat("\n • Number of genes left : ",nbCov)
if(noCov0){cat("\n Model building done without statistical tests.")}
cat("\n===========================================================================\n")


# Folder and paths  -------------------------------------------------------
pathToResults = paste0("outputs/buildingResults/application/ResultsPREVAC_batch")
dir(pathToResults)
pathToResults = paste0(pathToResults,"/Results_",buildMethod,
                       if((buildMethod %in% c("lassoSS","elasticnetSS"))
                          & length(thresholdsSS)!=1){"Crit"},
                       if(noCov0){"noCov0"})
dir(pathToResults)
pathToResults = paste0(pathToResults,"/",paste0(model,collapse=", "),"_",nbBatch,".RData")


# Temporary Files ---------------------------------------------------------
temporaryDirectory = paste0("/beegfs/agabaut/tmpPREVAC",nbCov,nbBatch)
if(dir.exists(temporaryDirectory)){
  warning("A folder already exists for this job. It has been deleted.")
  unlink(temporaryDirectory,recursive=T,force=T)
}
dir.create(temporaryDirectory)

# Job ---------------------------------------------------------------------
data = read.csv("data/applicationFiles/arm1/dataTrans_1000.txt")
genesNames = colnames(data)[5:1004]
sampledGenes = sample(genesNames,nbCov,replace = F)
rankGenes = tabulate(which(genesNames %in% sampledGenes),1000)
data = data %>% select(id,time,yAb,init,all_of(sampledGenes))
write.csv(data,file=paste0(temporaryDirectory,"/data",nbBatch,".txt"),quote = F,row.names = F)

newProject(modelFile = "data/modelFiles/PasinReg.txt",
           data = list(dataFile=paste0(temporaryDirectory,"/data",nbBatch,".txt"),
                       headerTypes = c("id","time","observation","regressor",rep("contcov",nbCov))))

setIndividualParameterVariability(delta_L=FALSE)
setErrorModel(yAb="constant")

load("data/applicationFiles/PopulationParameterInit.RData")

setPopulationParameterInformation(PopulationParameterInformation)

saveProject(paste0(temporaryDirectory,"/Build.mlxtran"))

res = buildmlx(project = paste0(temporaryDirectory,"/Build.mlxtran"),
               buildMethod = buildMethod,
               model=model,
               test=FALSE,
               thresholdsSS=thresholdsSS,
               p.max=if(noCov0){1}else{0.1})

Model <- Rsmlx:::mlx.getIndividualParameterModel()
covModel <- lapply(Model$covariateModel,FUN=function(cov){names(cov[which(cov)])})
time=res$time
iter=res$iter

save(Model,covModel,time,iter,sampledGenes,rankGenes,file=pathToResults)

unlink(temporaryDirectory,force = T,recursive = T)
