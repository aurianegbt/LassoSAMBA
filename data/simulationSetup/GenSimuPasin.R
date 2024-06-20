###############################################################################
#
#     Simulation with randomly distributed covariates 
#
###############################################################################
set.seed(81511807)
dir <- function(path){if(!dir.exists(path)){dir.create(path)}}


## Load the library required 
suppressMessages({
  library(dplyr)
  library(lixoftConnectors)
  initializeLixoftConnectors("simulx")
  loadProject("data/simulationSetup/Pasin.smlx")
  library(simstudy)
  library(data.table)
  library(ggcorrplot)
  library(foreach)
  library(ggpubr,quietly=TRUE)
  source("~/Travail/00_Theme.R")
  source("scripts/simulationFun/genCov.R")
  source("scripts/simulationFun/randomCovariate.R")
})

## Generate distribution 
distribution = randomCovariate(n=197)

## Correlation matrix
load("data/simulationSetup/distribPasin.RData")


## Generate 200 correlated covariates and then create 100 replicates

def <- defData(varname="AGE",formula =35,variance=16, dist = "normal")
def <- defData(def,varname="G1",formula=0,variance=1,dist="normal")
def <- defData(def,varname="G2",formula=0,variance=1,dist="normal")

for(i in 1:197){
  d= distribution[[i]]
  if(d$distribution=="gamma"){
    theta = d$elements$scale
    k = d$elements$shape
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=k*theta,variance = theta)
  }else if(d$distribution=="normal"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=d$elements$mean,variance=d$elements$sd)
  }else if(d$distribution=="poisson"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=d$elements$lambda)
  }else if(d$distribution=="uniform"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=paste0(d$elements,collapse = ";"))
  }
}
setNbReplicates(1)
genCorMat <- genCorMat + 1e-10 * diag(ncol(genCorMat))
covTableALL = genCorFlex(100*100,def,corMatrix=genCorMat)

for(i in 1:100){
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% mutate(id = (id-1)%%100+1)

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
  
  dir("data/simulationFiles/FilesPasin")
  dir("data/simulationFiles/FilesPasin/covTable")
  dir("data/simulationFiles/FilesPasin/simulation")
  
  if(i ==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesPasin/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesPasin/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}

headerTypes = c("id","time","observation",rep("contcov",3),sapply(sapply(distribution,FUN=function(x){x$distribution})=="poisson",FUN=function(x){if(x){"catcov"}else{"contcov"}}))
save(headerTypes,file="data/simulationFiles/FilesPasin/headerTypes.RData")
