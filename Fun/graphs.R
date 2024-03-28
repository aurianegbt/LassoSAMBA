suppressWarnings(suppressMessages(library(ggplot2,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggpubr,quietly=TRUE)))
suppressWarnings(suppressMessages(library(dplyr,quietly=TRUE)))
suppressWarnings(suppressMessages(library(webshot2,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggh4x,quietly=TRUE)))
suppressWarnings(suppressMessages(library(flextable,quietly=TRUE)))
suppressWarnings(suppressMessages(library(grid,quietly=TRUE)))
suppressWarnings(suppressMessages(library(gtable,quietly=TRUE)))
suppressWarnings(suppressMessages(library(gridExtra,quietly=TRUE)))
suppressWarnings(suppressMessages(library(scales,quietly=TRUE)))
sapply(paste0("Fun/graphs/",list.files("Fun/graphs")),source)

graphsGenerate <- function(project="Pasin",
                           covariateSize=200,
                           buildMethod=c("reg","lassoSS"),
                           JPEG = T,
                           PNG = F){
  eval(parse(text=readLines(paste0("Files/Files",project,"/info.txt"))))
  
  source("~/Travail/00_Theme.R")
  load(paste0("Save/BuildResults_",project,".RData"))
  
  # Load data
  
  for(sim in covariateSize){

    generalsubtitle = paste0("Among ",length(unique(errorStats$Model))," simulated datasets of ",dataType," for ",nInd," individuals, with ",sim," covariates.\n")
    Titlelist = list(reg = "Model built with SAMBA.",
                     lasso = "Model built with a lasso approach within SAMBA, whithout stability selection.",
                     elasticnet="Model built with a lasso approach within SAMBA, whithout stability selection.",
                     lassoSS = "Model built with a lasso approach within SAMBA.",
                     elasticnetSS = "Model built with an elastic net approach within SAMBA.",
                     rlasso= "Model built with a lasso approach within SAMBA, and s.s. on replicates.",
                     relasticnet="Model built with an elastic net approach within SAMBA, and on replicates.",
                     lassoCrit = "Model built with a lasso approach within SAMBA, and multiple thresholds, whithout stability selection.",
                     elasticnetCrit = "Model built with an elastic net approach within SAMBA, and multiple thresholds, whithout stability selection.",
                     lassoSSCrit = "Model built with a lasso approach within SAMBA, and multiple thresholds.",
                     elasticnetSSCrit ="Model built with an elastic net approach within SAMBA, and multiple thresholds",
                     rlassoCrit= "Model built with a lasso approach within SAMBA, and multiple thresholds and s.s. on replicates.",
                     relasticnetCrit="Model built with an elastic net approach within SAMBA, and multiple thresholds and s.s. on replicates.",
                     rsharp = "Model built with sharp method on replicates.",
                     sharp="Model built with sharp method.")[buildMethod[which(!stringr::str_detect(buildMethod,"regPEN") & !stringr::str_detect(buildMethod,"noCov0"))]]
    
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
    
    
    initFolder=paste0("~/Travail/Presentation/Plot/LassoSAMBA/Plot",project)
    if(!dir.exists(initFolder)){dir.create(initFolder)}
    initFolder=paste0(initFolder,"/Results",sim)
    if(!dir.exists(initFolder)){dir.create(initFolder)}
    
    
    
    # Graphs summarizing building result
    for(meth in buildMethod){
      Folder=paste0(initFolder,"/",meth)
      if(!dir.exists(Folder)){dir.create(Folder)}
      subtitle = paste0(generalsubtitle,Titlelist[[meth]])
      # graphsTotalNB(Folder,subtitle,project,sim,meth,JPEG,PNG)
      # graphsParNB(Folder,subtitle,project,sim,meth,JPEG,PNG)
      tableStats(Folder,subtitle,project,sim,meth,JPEG,PNG)
    }
    
    # Comparison graphs
    if(length(buildMethod)>1){
      Folder=initFolder
      subtitle = generalsubtitle
      # graphsCompMethod(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      # graphsParCompMethod(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      # graphsCompFR(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      tableStatsComp(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      graphsCompTime(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
    }
  }
}

