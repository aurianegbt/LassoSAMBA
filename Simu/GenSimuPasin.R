## On va générer l'ensemble des covariables puis créer les simulations selon celle choisi comme significative 
set.seed(1710)
# pb avec la seed ? tous les AB0 sont identiques ? sur cov pas gênant car juste constant sur réplicats ? mais sur corcov ils sont TOUS pareil... je pense à cause de la graine de simulx ? 
# -> TEST SAY NO, juste because same group each time ? 


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

doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()-5))

# correlation matrix
load("Simu/DataTransCoding.RData")

aux = (dataTransCoding %>% filter(visit=="M12"))[,c(7,9:18662)] 
genesKept = names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[3:201]

aux <- aux[,c("age",genesKept)]
colnames(aux) <- c(1:ncol(aux))

genCorMat <- cor(aux,method="spearman")
epsilon <- 1e-10
genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))
genCorMat[1, 1] <- 16

save(genCorMat,file="Simu/corrMatrixPasin.RData")

mu = c(35,0,0,rowMeans(t(aux[,-c(1:3)])))

# Plot
corrplot =  ggcorrplot(genCorMat[1:200,1:200],ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used, zoomed on the first 200 covariates.",
                                face="bold",size=30,color="#882255"))
ggsave("Simu/corrPasinZoom200.png",
       height = 6000, width = 6000, units = "px", bg='transparent')

corrplot =  ggcorrplot(genCorMat,ggtheme=theme_riri, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))


annotate_figure(corrplot,
                top = text_grob("Theoretical Correlation Matrix used.",
                                face="bold",size=50,color="#882255"))
ggsave("Simu/corrPasin.png",
       height = 10000, width = 10000, units = "px", bg='transparent')

## Generate 500 cov (correlated or not) and then create 500-200-50-10 simulation framework  
# cov

loadProject("Simu/Pasin.smlx")

runSimulation()

sim <- getSimulationResults()

dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

foreach(i = 1:100) %dopar% {
  library(dplyr)
  
  dataset = sim$res$yAB[sim$res$yAB$rep==i,c("id","time","yAB")]
  
  yAB0 = dataset[dataset$time==0,]
  colnames(yAB0)[3] <- "yAB0"
  
  dataset = merge(dataset,yAB0,by=c("id","time"),all.x = TRUE) 

  # covTable = merge(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,c("id","AGE","G1","G2")],genCov(distribution,individualsNumber = 100,covariateChar = "Gen"),by="id") %>% 
  cov = mvtnorm::rmvnorm(100,mean=mu[-c(1:3)])
  colnames(cov) <- paste0("Gen",1:197)
  cov <- cbind(id=1:100,cov)
  covTable =  merge(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,c("id","AGE","G1","G2")],cov,by="id") %>%
    mutate(AGE = AGE - mean(sim$IndividualParameters$simulationGroup1[sim$IndividualParameters$simulationGroup1$rep==i,"AGE"])) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",200))
  
  # ## 1000 save 
  # dir("Files/FilesPasin/1000cov")
  # dir("Files/FilesPasin/1000cov/covTable")
  # dir("Files/FilesPasin/1000cov/simulation")
  # 
  # write.csv(covTable,file=paste0("Files/FilesPasin/1000cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset,file=paste0("Files/FilesPasin/1000cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # save(headerTypes,file="Files/FilesPasin/1000cov/headerTypes.RData")
  
  ## 500 save
  # dir("Files/FilesPasin/500cov")
  # dir("Files/FilesPasin/500cov/covTable")
  # dir("Files/FilesPasin/500cov/simulation")
  # 
  # write.csv(covTable[,1:501],file=paste0("Files/FilesPasin/500cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset[,1:504],file=paste0("Files/FilesPasin/500cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # headerTypes <- headerTypes[1:504]
  # save(headerTypes,file="Files/FilesPasin/500cov/headerTypes.RData")
  
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
  
  # ##10 save 
  # dir("Files/FilesPasin/10cov")
  # dir("Files/FilesPasin/10cov/covTable")
  # dir("Files/FilesPasin/10cov/simulation")
  # 
  # write.csv(covTable[,1:11],file=paste0("Files/FilesPasin/10cov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset[,1:14],file=paste0("Files/FilesPasin/10cov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # headerTypes <- headerTypes[1:14]
  # save(headerTypes,file="Files/FilesPasin/10cov/headerTypes.RData")
}

# CORCOV
setNbReplicates(1)
covTableALL = as.data.frame(mvtnorm::rmvnorm(n=100*100,mean=mu,sigma = genCorMat))
colnames(covTableALL) <- c("AGE","G1","G2",paste0("Gen",1:197))

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
  
  yAB0 = dataset[dataset$time==0,]
  colnames(yAB0)[3] <- "yAB0"
  
  dataset = merge(dataset,yAB0,by=c("id","time"),all.x = TRUE) 
  
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% 
    mutate(AGE = AGE - mean(covTableALL[(1+(i-1)*100):(i*100),]$AGE)) %>%
    mutate(id=1:100,.before = AGE) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  headerTypes = c("id","time","observation","regressor",rep("contcov",200))

  # ## 1000 save
  # dir("Files/FilesPasin/1000corcov")
  # dir("Files/FilesPasin/1000corcov/covTable")
  # dir("Files/FilesPasin/1000corcov/simulation")
  # 
  # write.csv(covTable,file=paste0("Files/FilesPasin/1000corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset,file=paste0("Files/FilesPasin/1000corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # save(headerTypes,file="Files/FilesPasin/1000corcov/headerTypes.RData")
  
  # ## 500 save
  # dir("Files/FilesPasin/500corcov")
  # dir("Files/FilesPasin/500corcov/covTable")
  # dir("Files/FilesPasin/500corcov/simulation")
  # 
  # write.csv(covTable[,1:501],file=paste0("Files/FilesPasin/500corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset[,1:504],file=paste0("Files/FilesPasin/500corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # 
  # headerTypes <- headerTypes[1:504]
  # save(headerTypes,file="Files/FilesPasin/500corcov/headerTypes.RData")
  
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
  
  # ##10 save 
  # dir("Files/FilesPasin/10corcov")
  # dir("Files/FilesPasin/10corcov/covTable")
  # dir("Files/FilesPasin/10corcov/simulation")
  # 
  # write.csv(covTable[,1:11],file=paste0("Files/FilesPasin/10corcov/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  # write.csv(dataset[,1:14],file=paste0("Files/FilesPasin/10corcov/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  # 
  # headerTypes <- headerTypes[1:14]
  # save(headerTypes,file="Files/FilesPasin/10corcov/headerTypes.RData")
  
  unlink("tmpfile.txt")
}
