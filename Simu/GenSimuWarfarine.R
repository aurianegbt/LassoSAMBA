## On va générer l'ensemble des covariables puis créer les simulations selon celle choisi comme significative 
# generateSeed <- function(){
#   t <- as.numeric(Sys.time())
#   seed <- 1e8 * (t - floor(t))
#   return(seed)
# }
set.seed(37525892)


## Load the library required 
library(dplyr)
library(lixoftConnectors)
initializeLixoftConnectors("simulx")
library(simstudy)
library(foreach)
library(data.table)
library(ggcorrplot)
library(ggpubr,quietly=TRUE)
source("~/Travail/00_Theme.R")
source("Simu/Fun/genCov.R")
source("Simu/Fun/randomCovariate.R")
doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))

## Generate 500 cov (correlated or not) and then create 500-200-50-10 simulation framework  
# cov
distribution = randomCovariate(n=998,NORM = TRUE)

loadProject("Simu/Warfarine.smlx")

runSimulation()

sim <- getSimulationResults()

dir <- function(x){if(!dir.exists(x)){dir.create(x)}}


foreach(i = 1:100) %dopar% {
  library(dplyr)
  library(simstudy)
  dataset = sim$res$y1[sim$res$y1$rep==i,c("id","time","y1")]
  
  dataset = rbind(dataset,data.frame(id=unique(dataset$id),time=0,y1=".")) %>%
    arrange(id,time) 
  
  trt = sim$doses$simulationGroup1[sim$doses$simulationGroup1$rep==i,] %>%
    dplyr::select(id,time,amt)
  
  dataset <- merge(dataset,trt,all.x = TRUE,by=c("id","time")) %>%
    arrange(id,time,y1) 
    
  covTable = merge(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,c("id","AGE","WT")],genCov(distribution,individualsNumber = 100,covariateChar = "Gen"),by="id") %>% 
    mutate(AGE = AGE - mean(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,"AGE"])) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","amount",rep("contcov",1000))
  
  
  ## 1000 save
  dir("Files/FilesWarfarine/1000cov")
  dir("Files/FilesWarfarine/1000cov/covTable")
  dir("Files/FilesWarfarine/1000cov/simulation")
  
  write.csv(covTable,file=paste0("Files/FilesWarfarine/1000cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("Files/FilesWarfarine/1000cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="Files/FilesWarfarine/1000cov/headerTypes.RData")
  
  ## 500 save
  dir("Files/FilesWarfarine/500cov")
  dir("Files/FilesWarfarine/500cov/covTable")
  dir("Files/FilesWarfarine/500cov/simulation")
  
  write.csv(covTable[,1:501],file=paste0("Files/FilesWarfarine/500cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:504],file=paste0("Files/FilesWarfarine/500cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:504]
  save(headerTypes,file="Files/FilesWarfarine/500cov/headerTypes.RData")
  
  ## 200 save 
  dir("Files/FilesWarfarine/200cov")
  dir("Files/FilesWarfarine/200cov/covTable")
  dir("Files/FilesWarfarine/200cov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesWarfarine/200cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesWarfarine/200cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesWarfarine/200cov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesWarfarine/50cov")
  dir("Files/FilesWarfarine/50cov/covTable")
  dir("Files/FilesWarfarine/50cov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesWarfarine/50cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesWarfarine/50cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesWarfarine/50cov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesWarfarine/10cov")
  dir("Files/FilesWarfarine/10cov/covTable")
  dir("Files/FilesWarfarine/10cov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesWarfarine/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesWarfarine/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesWarfarine/10cov/headerTypes.RData")
}


## corcov
# Correlation matrix
load("Simu/DataTransCoding.RData")

aux = (dataTransCoding %>% filter(visit=="D63"))[,c(7,9:18662)] 
genesKept = names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[1:999]

aux <- aux[,c("age",genesKept)]
colnames(aux) <- c(1:ncol(aux))


genCorMat <- cor(aux,method="spearman")
epsilon <- 1e-10
genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))

save(genCorMat,file="Simu/corrMatrixWarfarine.RData")

corrplot =  ggcorrplot(genCorMat[1:200,1:200],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used, zoomed on the first 200 covariates.",
                                face="bold",size=30,color="#882255"))
ggsave("Simu/corrWarfarineZoom200.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

corrplot =  ggcorrplot(genCorMat,ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used.",
                                face="bold",size=50,color="#882255"))
ggsave("Simu/corrWarfarine.png",
       height = 10000, width = 10000, units = "px", bg='transparent')


# Generation 

def <- defData(varname = "AGE", formula = 0, variance = 0.3**2, dist = "normal")
def <- defData(def,varname="WT",formula=0,variance=0.2**2,dist="normal")
for(i in 1:998){
  def <- defData(def,varname=paste0("Gen",i),formula=rnorm(1),variance=rexp(1,0.3),dist="normal")
}
setNbReplicates(1)
save(dist,file="Simu/distSave.RData")
for(i in 1:100){
  covTable =  genCorFlex(100,def,corMatrix = genCorMat)  %>% 
    mutate(AGE = 29*exp(AGE)) %>%  # to simulate log-normal distribution
    mutate(WT = 68*exp(WT))   # as it is not straightforward in simstudy
  
  write.csv(covTable [,1:3],"tmpfile.txt",quote = F,row.names = F)
  
  defineCovariateElement(name=paste0("covTable",i),
                         element = "tmpfile.txt")
  
  setGroupElement(group="simulationGroup1", elements = c(paste0("covTable",i)))
  
  runSimulation()
  sim <- getSimulationResults()
  
  dataset = sim$res$y1[,c("id","time","y1")]
  
  dataset = rbind(dataset,data.frame(id=unique(dataset$id),time=0,y1=".")) %>%
    arrange(id,time) 
  
  trt = sim$doses$simulationGroup1 %>%
    dplyr::select(id,time,amt)
  
  dataset <- merge(dataset,trt,all.x = TRUE,by=c("id","time")) %>%
    arrange(id,time,y1) 
  
  covTable = covTable %>% 
    mutate(AGE = AGE - mean(covTable$AGE)) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","amount",rep("contcov",1000))
  
  ## 1000 save
  dir("Files/FilesWarfarine/1000corcov")
  dir("Files/FilesWarfarine/1000corcov/covTable")
  dir("Files/FilesWarfarine/1000corcov/simulation")
  
  write.csv(covTable,file=paste0("Files/FilesWarfarine/1000corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("Files/FilesWarfarine/1000corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="Files/FilesWarfarine/1000corcov/headerTypes.RData")
  
  ## 500 save
  dir("Files/FilesWarfarine/500corcov")
  dir("Files/FilesWarfarine/500corcov/covTable")
  dir("Files/FilesWarfarine/500corcov/simulation")
  
  write.csv(covTable[,1:501],file=paste0("Files/FilesWarfarine/500corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:504],file=paste0("Files/FilesWarfarine/500corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:504]
  save(headerTypes,file="Files/FilesWarfarine/500corcov/headerTypes.RData")
  
  ## 200 save 
  dir("Files/FilesWarfarine/200corcov")
  dir("Files/FilesWarfarine/200corcov/covTable")
  dir("Files/FilesWarfarine/200corcov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("Files/FilesWarfarine/200corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("Files/FilesWarfarine/200corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="Files/FilesWarfarine/200corcov/headerTypes.RData")
  ## 50 save 
  dir("Files/FilesWarfarine/50corcov")
  dir("Files/FilesWarfarine/50corcov/covTable")
  dir("Files/FilesWarfarine/50corcov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("Files/FilesWarfarine/50corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("Files/FilesWarfarine/50corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="Files/FilesWarfarine/50corcov/headerTypes.RData")
  
  ##10 save 
  dir("Files/FilesWarfarine/10corcov")
  dir("Files/FilesWarfarine/10corcov/covTable")
  dir("Files/FilesWarfarine/10corcov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("Files/FilesWarfarine/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("Files/FilesWarfarine/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="Files/FilesWarfarine/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}

