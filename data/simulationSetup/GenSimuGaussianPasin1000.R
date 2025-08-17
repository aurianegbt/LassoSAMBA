###############################################################################
#
#     Simulation with gaussian distributed covariates 
#
###############################################################################
set.seed(1710)
dir <- function(path){if(!dir.exists(path)){dir.create(path)}}

## Load the library required 
suppressMessages({
  library(dplyr)
  library(lixoftConnectors)
  library(simstudy)
  library(data.table)
  library(ggcorrplot)
  library(ggpubr,quietly=TRUE)
  
  initializeLixoftConnectors("simulx")
})
source("~/Travail/00_Theme.R")

load("data/simulationSetup/distribPasin1000.RData")
## Generate 1000 correlated covariates and then create 100 replicates
loadProject("data/simulationSetup/Pasin.smlx")
setNbReplicates(1)

covTableALL = as.data.frame(mvtnorm::rmvnorm(n=100*100,mean=mu,sigma = genCovMat))


colnames(covTableALL) <- c("AGE","G1","G2",paste0("Gen",1:997))

for(i in 1:100){
  covTable = covTableALL[(1+(i-1)*100):(i*100),]
  covTable <- cbind(id=1:100,covTable)
  
  write.csv(covTable[,1:4],paste0("tmpfile",i,".txt"),quote = F,row.names = F)
  
  defineCovariateElement(name=paste0("covTable",i),
                         element = paste0("tmpfile",i,".txt"))
  
  if(i==1){
    setGroupElement(group=paste0("simulationGroup",i), elements = c(paste0("covTable",i)))
  }else{
    addGroup(paste0("simulationGroup",i))
    setGroupElement(group=paste0("simulationGroup",i), elements = c(paste0("covTable",i)))
  }
}

runSimulation()
sim <- getSimulationResults()


for(i in 1:100){
  dataset = sim$res$yAB[sim$res$yAB$group==paste0("simulationGroup",i),c("original_id","time","yAB")] %>%
    rename(id=original_id)
  
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% 
    mutate(AGE = AGE - mean(covTableALL[(1+(i-1)*100):(i*100),]$AGE)) %>%
    mutate(id=1:100,.before = AGE) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  dir("data/simulationFiles/FilesGaussianPasin1000")
  dir("data/simulationFiles/FilesGaussianPasin1000/covTable")
  dir("data/simulationFiles/FilesGaussianPasin1000/simulation")
  
  if(i==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesGaussianPasin1000/covTable/covTable_",i,".txt"),quote = F,row.names = F) 
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesGaussianPasin1000/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}


headerTypes = c("id","time","observation",rep("contcov",1000))
save(headerTypes,file="data/simulationFiles/FilesGaussianPasin1000/headerTypes.RData")
