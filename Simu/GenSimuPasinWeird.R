## On va générer l'ensemble des covariables puis créer les simulations selon celle choisi comme significative 
set.seed(81511807)


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
source("Simu/Fun/randomCovariate2.R")

doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))


## Generate 500 cov (correlated or not) and then create 500-200-50-10 simulation framework  
# cov
distribution = randomCovariate(n=197)
headerTypes = sapply(sapply(distribution,FUN=function(x){x$distribution})=="poisson",FUN=function(x){if(x){"catcov"}else{"contcov"}})


loadProject("Simu/Pasin.smlx")

runSimulation()

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
  
  ## 200 save 
  dir("Files/FilesPasinWeird/200cov")
  dir("Files/FilesPasinWeird/200cov/covTable")
  dir("Files/FilesPasinWeird/200cov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesPasinWeird/200cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesPasinWeird/200cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesPasinWeird/200cov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesPasinWeird/50cov")
  dir("Files/FilesPasinWeird/50cov/covTable")
  dir("Files/FilesPasinWeird/50cov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesPasinWeird/50cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesPasinWeird/50cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesPasinWeird/50cov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesPasinWeird/10cov")
  dir("Files/FilesPasinWeird/10cov/covTable")
  dir("Files/FilesPasinWeird/10cov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesPasinWeird/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesPasinWeird/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesPasinWeird/10cov/headerTypes.RData")
}

## corcov
# load("Simu/DataTransCoding.RData")
# 
# aux = (dataTransCoding %>% filter(visit=="M12"))[,c(7,9:18662)] 
# genesKept = names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[1:999]
# 
# aux <- aux[,c("age",genesKept)]
# colnames(aux) <- c(1:ncol(aux))
# 
# genCorMat <- cor(aux,method="spearman")
# epsilon <- 1e-10
# genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))
# 
# save(genCorMat,file="Simu/corrMatrixPasinWeird.RData")
# 
# corrplot =  ggcorrplot(genCorMat[1:200,1:200],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used, zoomed on the first 200 covariates.",
#                                 face="bold",size=30,color="#882255"))
# ggsave("Simu/corrPasinWeirdZoom200.png",
#        height = 6000, width = 6000, units = "px", bg='transparent')
# 
# corrplot =  ggcorrplot(genCorMat,ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used.",
#                                 face="bold",size=50,color="#882255"))
# ggsave("Simu/corrPasinWeird.png",
#        height = 10000, width = 10000, units = "px", bg='transparent')
load("Simu/corrMatrixPasin.RData")
genCorMat <- genCorMat[1:200,1:200]


def <- defData(varname="AGE",formula =35,variance=8, dist = "normal")
def <- defData(def,varname="G1",formula=0,variance=1,dist="normal")
def <- defData(def,varname="G2",formula=0,variance=1,dist="normal")
for(i in 1:197){
  # def <- defData(def,varname=paste0("Gen",i),formula=rnorm(1),variance=rexp(1,0.3),dist="normal")
  
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
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",200))
  
  ## 200 save 
  dir("Files/FilesPasinWeird/200corcov")
  dir("Files/FilesPasinWeird/200corcov/covTable")
  dir("Files/FilesPasinWeird/200corcov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesPasinWeird/200corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesPasinWeird/200corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesPasinWeird/200corcov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesPasinWeird/50corcov")
  dir("Files/FilesPasinWeird/50corcov/covTable")
  dir("Files/FilesPasinWeird/50corcov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesPasinWeird/50corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesPasinWeird/50corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesPasinWeird/50corcov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesPasinWeird/10corcov")
  dir("Files/FilesPasinWeird/10corcov/covTable")
  dir("Files/FilesPasinWeird/10corcov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesPasinWeird/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesPasinWeird/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesPasinWeird/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}
