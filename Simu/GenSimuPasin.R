## On va générer l'ensemble des covariables puis créer les simulations selon celle choisi comme significative 

## Load the library required 
library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors("simulx")
library(simstudy)
library(data.table)
library(ggcorrplot)
library(ggpubr,quietly=TRUE)
source("~/Travail/00_Theme.R")
source("Simu/Fun/genCov.R")
source("Simu/Fun/randomCovariate.R")

## Generate 500 cov (correlated or not) and then create 500-200-50-10 simulation framework  
# cov
distribution = randomCovariate(n=497,NORM = TRUE)

loadProject("Simu/Pasin.smlx")

sim <- getSimulationResults()

dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

for(i in 1:100){
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
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",500))
  
  
  ## 500 save
  dir("FilesPasin/500cov")
  dir("FilesPasin/500cov/covTable")
  dir("FilesPasin/500cov/simulation")
  
  write.csv(covTable,file=paste0("FilesPasin/500cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("FilesPasin/500cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="FilesPasin/500cov/headerTypes.RData")
  
  ## 200 save 
  dir("FilesPasin/200cov")
  dir("FilesPasin/200cov/covTable")
  dir("FilesPasin/200cov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("FilesPasin/200cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("FilesPasin/200cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="FilesPasin/200cov/headerTypes.RData")
  ## 50 save 
  dir("FilesPasin/50cov")
  dir("FilesPasin/50cov/covTable")
  dir("FilesPasin/50cov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("FilesPasin/50cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("FilesPasin/50cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="FilesPasin/50cov/headerTypes.RData")
  
  ##10 save 
  dir("FilesPasin/10cov")
  dir("FilesPasin/10cov/covTable")
  dir("FilesPasin/10cov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("FilesPasin/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("FilesPasin/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="FilesPasin/10cov/headerTypes.RData")
}

## corcov
genCorMat = genCorMat(500)
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
for(i in 1:497){
  def <- defData(def,varname=paste0("Gen",i),formula=rnorm(1),variance=rexp(1,0.3),dist="normal")
}
setNbReplicates(1)

for(i in 1:100){
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
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",500))

  ## 500 save
  dir("FilesPasin/500corcov")
  dir("FilesPasin/500corcov/covTable")
  dir("FilesPasin/500corcov/simulation")
  
  write.csv(covTable,file=paste0("FilesPasin/500corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("FilesPasin/500corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="FilesPasin/500corcov/headerTypes.RData")
  
  ## 200 save 
  dir("FilesPasin/200corcov")
  dir("FilesPasin/200corcov/covTable")
  dir("FilesPasin/200corcov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("FilesPasin/200corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("FilesPasin/200corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="FilesPasin/200corcov/headerTypes.RData")
  ## 50 save 
  dir("FilesPasin/50corcov")
  dir("FilesPasin/50corcov/covTable")
  dir("FilesPasin/50corcov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("FilesPasin/50corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("FilesPasin/50corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="FilesPasin/50corcov/headerTypes.RData")
  
  ##10 save 
  dir("FilesPasin/10corcov")
  dir("FilesPasin/10corcov/covTable")
  dir("FilesPasin/10corcov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("FilesPasin/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("FilesPasin/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="FilesPasin/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}