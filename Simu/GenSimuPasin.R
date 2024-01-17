## On va générer l'ensemble des covariables puis créer les simulations selon celle choisi comme significative 
set.seed(1710)


## Load the library required 
library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors("simulx")
library(simstudy)
library(data.table)
library(ggcorrplot)
library(foreach)
library(ggpubr,quietly=TRUE)
source("~/Travail/00_Theme.R")
source("Simu/Fun/genCov.R")
source("Simu/Fun/randomCovariate.R")

doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))


## Generate 500 cov (correlated or not) and then create 500-200-50-10 simulation framework  
# cov
distribution = randomCovariate(n=997,NORM = TRUE)

loadProject("Simu/Pasin.smlx")

sim <- getSimulationResults()

dir <- function(x){if(!dir.exists(x)){dir.create(x)}}



foreach(i = 1:100) %dopar% {
  library(dplyr)
  library(simstudy)
  
  dataset = sim$res$yAB[sim$res$yAB$rep==i,c("id","time","yAB")]
  
  yAB0 = dataset[dataset$time==0,]
  colnames(yAB0)[3] <- "yAB0"
  
  dataset = merge(dataset,yAB0,by=c("id","time"),all.x = TRUE) 

  covTable = merge(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,c("id","AGE","G1","G2")],genCov(distribution,individualsNumber = 100,covariateChar = "Gen"),by="id") %>% 
    mutate(AGE = AGE - mean(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,"AGE"])) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",1000))
  
  ## 1000 save 
  dir("Files/FilesPasin/1000cov")
  dir("Files/FilesPasin/1000cov/covTable")
  dir("Files/FilesPasin/1000cov/simulation")
  
  write.csv(covTable,file=paste0("Files/FilesPasin/1000cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("Files/FilesPasin/1000cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="Files/FilesPasin/1000cov/headerTypes.RData")
  
  ## 500 save
  dir("Files/FilesPasin/500cov")
  dir("Files/FilesPasin/500cov/covTable")
  dir("Files/FilesPasin/500cov/simulation")
  
  write.csv(covTable[,1:501],file=paste0("Files/FilesPasin/500cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:504],file=paste0("Files/FilesPasin/500cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:504]
  save(headerTypes,file="Files/FilesPasin/500cov/headerTypes.RData")
  
  ## 200 save 
  dir("Files/FilesPasin/200cov")
  dir("Files/FilesPasin/200cov/covTable")
  dir("Files/FilesPasin/200cov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesPasin/200cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesPasin/200cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesPasin/200cov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesPasin/50cov")
  dir("Files/FilesPasin/50cov/covTable")
  dir("Files/FilesPasin/50cov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesPasin/50cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesPasin/50cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesPasin/50cov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesPasin/10cov")
  dir("Files/FilesPasin/10cov/covTable")
  dir("Files/FilesPasin/10cov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesPasin/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesPasin/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesPasin/10cov/headerTypes.RData")
}

## corcov
genCorMat = genCorMat(1000)
corrplot = ggcorrplot(genCorMat,ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used",
                                face="bold",size=20,color="#882255"))
ggsave("Simu/corrPasin.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

corrplot = ggcorrplot(genCorMat[1:50,1:50],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot)
ggsave("Simu/corrPasinZoom.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

def <- defData(varname="AGE",formula = "20;50", dist = "uniform")
def <- defData(def,varname="G1",formula=0,variance=1,dist="normal")
def <- defData(def,varname="G2",formula=0,variance=1,dist="normal")
for(i in 1:997){
  def <- defData(def,varname=paste0("Gen",i),formula=rnorm(1),variance=rexp(1,0.3),dist="normal")
}
setNbReplicates(1)

for(i in 1:100){
  library(dplyr)
  library(lixoftConnectors)
  library(simstudy)
  covTable = genCorFlex(100,def,corMatrix=genCorMat)
  
  write.csv(covTable [,1:4],"tmpfile.txt",quote = F,row.names = F)
  
  defineCovariateElement(name=paste0("covTable",i),
                         element = "tmpfile.txt")
  
  setGroupElement(group="simulationGroup1", elements = c(paste0("covTable",i)))
  
  runSimulation()
  sim <- getSimulationResults()
  
  dataset = sim$res$yAB[,c("id","time","yAB")]
  
  yAB0 = dataset[dataset$time==0,]
  colnames(yAB0)[3] <- "yAB0"
  
  dataset = merge(dataset,yAB0,by=c("id","time"),all.x = TRUE) 
  
  covTable = covTable %>% 
    mutate(AGE = AGE - mean(covTable$AGE)) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",1000))

  ## 1000 save
  dir("Files/FilesPasin/1000corcov")
  dir("Files/FilesPasin/1000corcov/covTable")
  dir("Files/FilesPasin/1000corcov/simulation")
  
  write.csv(covTable,file=paste0("Files/FilesPasin/1000corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("Files/FilesPasin/1000corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="Files/FilesPasin/1000corcov/headerTypes.RData")
  
  ## 500 save
  dir("Files/FilesPasin/500corcov")
  dir("Files/FilesPasin/500corcov/covTable")
  dir("Files/FilesPasin/500corcov/simulation")
  
  write.csv(covTable[,1:501],file=paste0("Files/FilesPasin/500corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:504],file=paste0("Files/FilesPasin/500corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  
  headerTypes <- headerTypes[1:504]
  save(headerTypes,file="Files/FilesPasin/500corcov/headerTypes.RData")
  
  ## 200 save 
  dir("Files/FilesPasin/200corcov")
  dir("Files/FilesPasin/200corcov/covTable")
  dir("Files/FilesPasin/200corcov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesPasin/200corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesPasin/200corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesPasin/200corcov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesPasin/50corcov")
  dir("Files/FilesPasin/50corcov/covTable")
  dir("Files/FilesPasin/50corcov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesPasin/50corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesPasin/50corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesPasin/50corcov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesPasin/10corcov")
  dir("Files/FilesPasin/10corcov/covTable")
  dir("Files/FilesPasin/10corcov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesPasin/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesPasin/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesPasin/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}