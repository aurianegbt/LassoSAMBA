suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(webshot2)
  library(ggh4x)
  library(flextable)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(scales)
  sapply(list.files("scripts/resultsFun/graphsFun",full.names = T),FUN=function(d){source(d,echo=F)})
})

graphsGenerate <- function(project="Pasin",
                           buildMethod="all",
                           JPEG = T,
                           PNG = F){
  
  if(identical(buildMethod,"all")){
    buildMethod <- stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",project),recursive = F),paste0("outputs/buildingResults/simulation/Results",project,"/Results_"))
    
    buildMethod = c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassoSSnoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit")[which( c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassoSSnoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit") %in% buildMethod)]
    cat("For ",project," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
  }

  
  
  
  eval(parse(text=readLines("data/simulationFiles/info.txt")))
  
  source("~/Travail/00_Theme.R")
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  
  # Load data
  generalsubtitle = paste0("Among ",length(unique(errorStats$Model))," simulated datasets of ",dataType," for ",nInd," individuals, with 200 correlated covariates.\n")
  Titlelist = list(reg = "Model built with SAMBA.",
                   lasso = "Model built with a lasso approach within SAMBA, whithout stability selection.",
                   elasticnet="Model built with a lasso approach within SAMBA, whithout stability selection.",
                   lassoSS = "Model built with a lasso approach within SAMBA.",
                   elasticnetSS = "Model built with an elastic net approach within SAMBA.",
                   sharp = "Model built with a lasso approach, calibrated using sharp method.")[buildMethod[which(!stringr::str_detect(buildMethod,"regPEN") & !stringr::str_detect(buildMethod,"noCov0"))]]
  
  for(k in 1:length(buildMethod)){
    if(stringr::str_detect(buildMethod[k],"regPEN")){
      Titlelist <- append(Titlelist, paste0("Model built with SAMBA, whose criterion is ", stringr::str_remove(buildMethod[k],"regPEN")," times penalized."))
      names(Titlelist)[length(Titlelist)] <- buildMethod[k]
    }
    if(stringr::str_detect(buildMethod[k],"noCov0")){
      Titlelist <- append(Titlelist, paste0(stringr::str_sub(Titlelist[stringr::str_remove(buildMethod[k],"noCov0")],end=-2),", whithout statistical test to exclude covariates."))
      names(Titlelist)[length(Titlelist)] <- buildMethod[k]
    }
  }
  
  
  initFolder=paste0("outputs/figures/simulationResults/Plot",project)
  if(!dir.exists(initFolder)){dir.create(initFolder)}
  if(!dir.exists(paste0(initFolder,"/SOLO"))){dir.create(paste0(initFolder,"/SOLO"))}
  if(!dir.exists(paste0(initFolder,"/COMP"))){dir.create(paste0(initFolder,"/COMP"))}
  
  
  # Graphs summarizing building result
  for(meth in buildMethod){
    Folder=paste0(initFolder,"/SOLO/",meth)
    if(!dir.exists(Folder)){dir.create(Folder)}
    subtitle = paste0(generalsubtitle,Titlelist[[meth]])
    graphsParNB(Folder,subtitle,project,meth,JPEG,PNG)
    tableStats(Folder,subtitle,project,meth,JPEG,PNG)
  }
  
  # Comparison graphs
  if(length(buildMethod)>1){
    Folder=paste0(initFolder,"/COMP/",paste0(buildMethod,collapse="_"))
    if(!dir.exists(Folder)){dir.create(Folder)}
    subtitle = generalsubtitle
    graphsParCompMethod(Folder,subtitle,project,buildMethod,JPEG,PNG)
    graphsCompFR(Folder,subtitle,project,buildMethod,JPEG,PNG)
    tableStatsComp(Folder,subtitle,project,buildMethod,JPEG,PNG)
    graphsCompTime(Folder,subtitle,project,buildMethod,JPEG,PNG)
  }
}

