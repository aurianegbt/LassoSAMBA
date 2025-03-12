# Generation of simulated datasets for Logistic Growth model framework with standard gaussian  covariates, completely based on Naveau et al., 2023

*This simulation setup is extracted from Marion Naveau et al. work, available on her github https://github.com/Marion-Naveau/Supp_Information_SAEMVS.*

This simulation are simply dowloaded from Marion Naveau github (see Saves folder form V_simul and data_comp Rdata object that contain respectively the covariates and the individual observations of this simulation framework).

```r
dir <- function(path){if(!dir.exists(path)){dir.create(path)}}

load("data/simulationSetup/data_comp.Rdata")
load("data/simulationSetup/V_simul.Rdata")

for(i in 1:100){
  
  covTable <- cbind(id=1:200,V_simul[[i]])
  colnames(covTable)[-1] <-  paste0("C",1:500)

  dataset = merge(cbind(data[[i]]),covTable,by.x="Id",by.y="id")
  
  dir("data/simulationFiles/FilesNaveau")
  dir("data/simulationFiles/FilesNaveau/covTable")
  dir("data/simulationFiles/FilesNaveau/simulation")
  
  if(i ==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesNaveau/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesNaveau/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
}
```
