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
distribution = randomCovariate(n=498,NORM = TRUE)

loadProject("Simu/Warfarine.smlx")

sim <- getSimulationResults()

dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

for(i in 1:100){
  dataset = sim$res$y1[sim$res$y1$rep==i,c("id","time","y1")]
  
  dataset = rbind(dataset,data.frame(id=unique(dataset$id),time=0,y1=".")) %>%
    arrange(id,time) 
  
  trt = sim$doses$simulationGroup1[sim$doses$simulationGroup1$rep==i,] %>%
    select(id,time,amt)
  
  dataset <- merge(dataset,trt,all.x = TRUE,by=c("id","time")) %>%
    arrange(id,time,y1) 
    
  covTable = merge(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,c("id","AGE","WT")],genCov(distribution,individualsNumber = 100,covariateChar = "Gen"),by="id") %>% 
    mutate(AGE = AGE - mean(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,"AGE"])) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","amount",rep("contcov",500))
  
  
  ## 500 save
  dir("FilesWarfarine/500cov")
  dir("FilesWarfarine/500cov/covTable")
  dir("FilesWarfarine/500cov/simulation")
  
  write.csv(covTable,file=paste0("FilesWarfarine/500cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("FilesWarfarine/500cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="FilesWarfarine/500cov/headerTypes.RData")
  
  ## 200 save 
  dir("FilesWarfarine/200cov")
  dir("FilesWarfarine/200cov/covTable")
  dir("FilesWarfarine/200cov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("FilesWarfarine/200cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("FilesWarfarine/200cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="FilesWarfarine/200cov/headerTypes.RData")
  ## 50 save 
  dir("FilesWarfarine/50cov")
  dir("FilesWarfarine/50cov/covTable")
  dir("FilesWarfarine/50cov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("FilesWarfarine/50cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("FilesWarfarine/50cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="FilesWarfarine/50cov/headerTypes.RData")
  
  ##10 save 
  dir("FilesWarfarine/10cov")
  dir("FilesWarfarine/10cov/covTable")
  dir("FilesWarfarine/10cov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("FilesWarfarine/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("FilesWarfarine/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="FilesWarfarine/10cov/headerTypes.RData")
}


## corcov
# Correlation matrix
load("Simu/DataTransCoding500SinceD63.RData")

aux = rbind(arm1,arm2,arm5,arm6,arm7)  %>%
  filter(visit=="D63") %>%
  select(AGE,WEIGHT,hla_drb4:rbpms2)

colnames(aux) <- c(1:500)

corMatrix <- cor(aux,method="spearman")
epsilon <- 1e-10
corMatrix <- corMatrix + epsilon * diag(ncol(corMatrix))

corrplot = ggcorrplot(corMatrix,ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used",
                                face="bold",size=20,color="#882255"))
ggsave("Simu/corrPK.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

corrplot = ggcorrplot(corMatrix[1:50,1:50],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain"))
annotate_figure(corrplot)
ggsave("Simu/corrPKZoom1.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

corrplot = ggcorrplot(corMatrix[1:150,1:150],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain"))
annotate_figure(corrplot)
ggsave("Simu/corrPKZoom2.png",
       height = 6000, width = 6000, units = "px", bg='transparent')


# Generation 

def <- defData(varname = "AGE", formula = 0, variance = 0.3**2, dist = "normal")
def <- defData(def,varname="WT",formula=0,variance=0.2**2,dist="normal")
for(i in 1:498){
  def <- defData(def,varname=paste0("Gen",i),formula=rnorm(1),variance=rexp(1,0.3),dist="normal")
}
setNbReplicates(1)

for(i in 1:100){
  covTable =  genCorFlex(100,def,corMatrix = corMatrix)  %>% 
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
    select(id,time,amt)
  
  dataset <- merge(dataset,trt,all.x = TRUE,by=c("id","time")) %>%
    arrange(id,time,y1) 
  
  covTable = covTable %>% 
    mutate(AGE = AGE - mean(covTable$AGE)) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","amount",rep("contcov",500))
  
  ## 500 save
  dir("FilesWarfarine/500corcov")
  dir("FilesWarfarine/500corcov/covTable")
  dir("FilesWarfarine/500corcov/simulation")
  
  write.csv(covTable,file=paste0("FilesWarfarine/500corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset,file=paste0("FilesWarfarine/500corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  save(headerTypes,file="FilesWarfarine/500corcov/headerTypes.RData")
  
  ## 200 save 
  dir("FilesWarfarine/200corcov")
  dir("FilesWarfarine/200corcov/covTable")
  dir("FilesWarfarine/200corcov/simulation")
  
  write.csv(covTable[,1:201],file=paste0("FilesWarfarine/200corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:204],file=paste0("FilesWarfarine/200corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:204]
  save(headerTypes,file="FilesWarfarine/200corcov/headerTypes.RData")
  ## 50 save 
  dir("FilesWarfarine/50corcov")
  dir("FilesWarfarine/50corcov/covTable")
  dir("FilesWarfarine/50corcov/simulation")
  
  write.csv(covTable[,1:51],file=paste0("FilesWarfarine/50corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:54],file=paste0("FilesWarfarine/50corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:54]
  save(headerTypes,file="FilesWarfarine/50corcov/headerTypes.RData")
  
  ##10 save 
  dir("FilesWarfarine/10corcov")
  dir("FilesWarfarine/10corcov/covTable")
  dir("FilesWarfarine/10corcov/simulation")
  
  write.csv(covTable[,1:11],file=paste0("FilesWarfarine/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  write.csv(dataset[,1:14],file=paste0("FilesWarfarine/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  headerTypes <- headerTypes[1:14]
  save(headerTypes,file="FilesWarfarine/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}

