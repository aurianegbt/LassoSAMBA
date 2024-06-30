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
                           exclude=c("reg","lassonoCov010"),
                           JPEG = T,
                           PNG = F){
  
  if(identical(buildMethod,"all")){
    buildMethod <- stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",project),recursive = F),paste0("outputs/buildingResults/simulation/Results",proj,"/Results_"))
    buildMethod = c("reg",
                    buildMethod[stringr::str_detect(buildMethod,"lassonoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"lassonoCov0FDP"))])[which(c("reg", buildMethod[stringr::str_detect(buildMethod,"lassonoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"lassonoCov0FDP"))]) %in% buildMethod)]
    if(!is.null(exclude)){
      buildMethod = setdiff(buildMethod,exclude)
    }
    
    cat("For ",project," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
  }

  eval(parse(text=readLines(paste0("data/simulationFiles/Files",project,"/info.txt"))))
  
  source("~/Travail/00_Theme.R")
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  
  # Load data
  generalsubtitle = paste0("Among ",length(unique(errorStats$Model))," simulated datasets of ",dataType," for ",nInd," individuals, with ",nbCov," correlated covariates.\n")
  Titlelist = list(reg = "Model built with SAMBA, with statistic test to exclude covariates at each iteration.",
                   lasso = "Model built with a lasso approach within SAMBA, with statistic test to exclude covariates at each iteration.")[Reduce(union,c("reg",if(any(stringr::str_detect(buildMethod,"lasso"))){"lasso"}))]
  
  for(k in 1:length(buildMethod)){
    if(stringr::str_detect(buildMethod[k],"lassonoCov0FDP") && grepl("^[0-9]+$", stringr::str_remove(buildMethod[k],"lassonoCov0FDP"))){
      Titlelist <- append(Titlelist, paste0(stringr::str_remove(Titlelist["lasso"],", with statistic test to exclude covariates at each iteration.")," and constrained E[FDR] under",stringr::str_remove(buildMethod[k],"lassonoCov0FDP"),"%."))
      names(Titlelist)[length(Titlelist)] <- buildMethod[k]
    
  }
  Titlelist <- Titlelist[buildMethod]
  
  
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
}

