###############################################################################
#
#     Simulation with randomly distributed covariates and numerous covariates 
#
###############################################################################
set.seed(61651537)
dir <- function(path){if(!dir.exists(path)){dir.create(path)}}


## Load the library required 
suppressMessages({
  library(dplyr)
  library(lixoftConnectors)
  initializeLixoftConnectors("simulx")
  loadProject("data/simulationSetup/CplxPasin.smlx")
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
load("data/applicationFiles/DataTransD63.RData")
load("data/simulationSetup/ComplexPasin/genesKept.RData")

aux = genesD63[,-1] 
genesKept = c(genesKept,setdiff(names(sort(apply(aux,2,sd),decreasing = TRUE)),genesKept))[1:200]

aux <- aux %>% select(all_of(genesKept))
# colnames(aux) <- c(1:ncol(aux))

genCorMat <- cor(aux,method="spearman")
epsilon <- 1e-10
genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))


# diag(genCorMat)[1:19] <- sd_genes**2
sd_genes = apply(aux[,1:19],2,sd)
mu <- apply(aux[1:19],2,mean)
save(genCorMat,file="data/simulationSetup/corrMatrixCplxPasin.RData")

# Draw distribution
distribution = randomCovariate(n=(200-19))

## Generate 200 correlated covariates and then create 100 replicates

def <- defData(varname="G1",formula = mu[1],variance=sd_genes[1]**2, dist = "normal")
for(i in 2:19){
  def <- defData(def,varname=paste0("G",i),formula=mu[i],variance=sd_genes[i]**2,dist="normal")
}
for(i in 20:200){
  d= distribution[[i-19]]
  if(d$distribution=="gamma"){
    theta = d$elements$scale
    k = d$elements$shape
    def <- defData(def,varname=paste0("G",i),dist=d$distribution,
                   formula=k*theta,variance = theta)
  }else if(d$distribution=="normal"){
    def <- defData(def,varname=paste0("G",i),dist=d$distribution,
                   formula=d$elements$mean,variance=d$elements$sd)
  }else if(d$distribution=="poisson"){
    def <- defData(def,varname=paste0("G",i),dist=d$distribution,
                   formula=d$elements$lambda)
  }else if(d$distribution=="uniform"){
    def <- defData(def,varname=paste0("G",i),dist=d$distribution,
                   formula=paste0(d$elements,collapse = ";"))
  }
}
setNbReplicates(1)
covTableALL = genCorFlex(100*100,def,corMatrix=genCorMat)

for(i in 1:100){
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% mutate(id = (id-1)%%100+1)

  write.csv(covTable[,1:17],paste0("tmpfile",i,".txt"),quote = F,row.names = F)
  
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
  dataset = sim$res$yAb[sim$res$yAb$group==paste0("simulationGroup",i),c("original_id","time","yAb")] %>%
    rename(id=original_id,yAB=yAb)
  
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% 
    mutate(id=1:100)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  dir("data/simulationFiles/FilesCplxPasin")
  dir("data/simulationFiles/FilesCplxPasin/covTable")
  dir("data/simulationFiles/FilesCplxPasin/simulation")
  
  write.csv(covTable,file=paste0("data/simulationFiles/FilesCplxPasin/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("data/simulationFiles/FilesCplxPasin/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}
headerTypes = c("id","time","observation",rep("contcov",19),sapply(sapply(distribution,FUN=function(x){x$distribution})=="poisson",FUN=function(x){if(x){"catcov"}else{"contcov"}}))

save(headerTypes,file="data/simulationFiles/FilesCplxPasin/headerTypes.RData")
