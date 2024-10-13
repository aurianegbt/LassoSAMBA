dir <- function(path){if(!dir.exists(path)){dir.create(path)}}

load("data/simulationSetup/data_comp.Rdata")
load("data/simulationSetup/V_simul.Rdata")

for(i in 1:100){
  
  covTable <- cbind(id=1:200,V_simul[[i]])
  colnames(covTable)[-1] <-  paste0("C",1:500)

  dataset = merge(cbind(data[[i]]),covTable,by.x="Id",by.y="id")
  
  dir("data/simulationFiles/FilesNAVEAU2")
  dir("data/simulationFiles/FilesNAVEAU2/covTable")
  dir("data/simulationFiles/FilesNAVEAU2/simulation")
  
  if(i ==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesNAVEAU2/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesNAVEAU2/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
}
