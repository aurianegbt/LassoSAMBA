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


# Correlation matrix
# load("data/applicationFiles/arm1/DataTransCoding.RData")
# 
# aux = (dataTransCoding %>% filter(visit=="M12"))[,c(6,8:16894)]
# genesKept = sample(names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[3:201]) # keep only 199 - age left aside
# 
# aux <- aux[,c("age",genesKept)]
# colnames(aux) <- c(1:ncol(aux))
# 
# genCorMat <- cor(aux,method="spearman")
# mu = c(35,0,0,apply(aux[,-c(1:3)],2,mean)) # AGE, G1, G2 have already means defined in simulX
# sd = diag(c(4,1,1,apply(aux[,-c(1:3)],2,sd)))
# 
# genCovMat = sd %*% genCorMat %*% sd
# 
# save(mu,genCovMat,genCorMat,file="data/simulationSetup/distribPasin.RData")
# 
# ## Plot
# corrplot =  ggcorrplot(genCorMat, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used.",
#                                 face="bold",size=50,color="#882255"))
# ggsave("outputs/figures/explanatory/corrPasin.png",
#        height = 10000, width = 10000, units = "px", bg='transparent')
# 
# corrplot =  ggcorrplot(genCorMat[1:50,1:50], colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used, zoom on the first 50 covariates.",
#                                 face="bold",size=50,color="#882255"))
# ggsave("outputs/figures/explanatory/corrPasinzoom.png",
#        height = 10000, width = 10000, units = "px", bg='transparent')

load("data/simulationSetup/distribPasin.RData")
## Generate 200 correlated covariates and then create 100 replicates
loadProject("data/simulationSetup/Pasin.smlx")
setNbReplicates(1)

covTableALL = as.data.frame(mvtnorm::rmvnorm(n=100*100,mean=mu,sigma = genCovMat))
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
  
  covTable = covTableALL[(1+(i-1)*100):(i*100),] %>% 
    mutate(AGE = AGE - mean(covTableALL[(1+(i-1)*100):(i*100),]$AGE)) %>%
    mutate(id=1:100,.before = AGE) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  dir("data/simulationFiles/FilesGaussianPasin")
  dir("data/simulationFiles/FilesGaussianPasin/covTable")
  dir("data/simulationFiles/FilesGaussianPasin/simulation")
  
  if(i==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesGaussianPasin/covTable/covTable_",i,".txt"),quote = F,row.names = F) 
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesGaussianPasin/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}


headerTypes = c("id","time","observation",rep("contcov",200))
save(headerTypes,file="data/simulationFiles/FilesGaussianPasin/headerTypes.RData")
