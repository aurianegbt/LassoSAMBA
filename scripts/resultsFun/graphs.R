graphsGenerate <- function(project="Pasin",
                           buildMethod="all",
                           JPEG = T,
                           PNG = F){
  
  if(identical(buildMethod,"all")){
    buildMethod <- setdiff(stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",proj),recursive = F),paste0("outputs/buildingResults/simulation/Results",proj,"/Results_")),paste0("outputs/buildingResults/simulation/Results",proj,"/gatheredResults"))
    cat("For ",project," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
  }

  eval(parse(text=readLines(paste0("data/simulationFiles/Files",project,"/info.txt"))))
  
  source("scripts/resultsFun/00_Theme.R")
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  
  # Load data
  generalsubtitle = paste0("Among ",length(unique(errorStatsPar$Model))," simulated datasets of ",dataType," for ",nInd," individuals, with ",nbCov," correlated covariates.\n")
  Titlelist = list(stepAIC = "Model built with step-SAMBA.",
                   lasso = "Model built with lasso-SAMBA method.")[Reduce(union,c("stepAIC",if(any(stringr::str_detect(buildMethod,"lasso"))){"lasso"}))] 
  
  for(k in 1:length(buildMethod)){
    if(stringr::str_detect(buildMethod[k],"lassoFDP")){
      Titlelist <- append(Titlelist, paste0(stringr::str_remove(Titlelist["lasso"],"\\.")," with constrained E[FDR] under ",stringr::str_remove(buildMethod[k],"lassoFDP"),"%."))
      names(Titlelist)[length(Titlelist)] <- buildMethod[k]
    }
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
    p1 = graphsParCompMethod(Folder,subtitle,project,buildMethod,JPEG,PNG)
    # graphsCompFR(Folder,subtitle,project,buildMethod,JPEG,PNG)
    tableStatsComp(Folder,subtitle,project,buildMethod,JPEG,PNG)
    p2 = graphsStatsComp(Folder,subtitle,project,buildMethod,JPEG,PNG)
    if(length(setdiff(buildMethod,"SAEMVS"))>1){
      p3 = graphsCompTime(Folder,subtitle,project,setdiff(buildMethod,"SAEMVS"),JPEG,PNG)
      p4 = graphsLL(Folder,subtitle,project,setdiff(buildMethod,"SAEMVS"),JPEG,PNG)
    }
  }
  
  return(list(ParComp = p1,
              StatsComp = p2,
              TimeComp = p3,
              LLComp = p4))
}

