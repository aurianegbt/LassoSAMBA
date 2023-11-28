graphsGenerate <- function(project="Pasin",
                           covariateSize=c(10,50,200,500),
                           buildMethod=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"),
                           JPEG = FALSE,
                           PNG = TRUE){
  eval(parse(text=readLines(paste0("Files",project,"/info.txt"))))
  
  library(ggplot2,quietly=TRUE)
  library(ggpubr,quietly=TRUE)
  library(webshot2,quietly=TRUE)
  library(ggh4x,quietly=TRUE)
  suppressMessages(library(flextable,quietly=TRUE))
  library(grid,quietly=TRUE)
  library(gtable,quietly=TRUE)
  library(gridExtra,quietly=TRUE)
  library(scales,quietly=TRUE)
  source("~/Travail/00_Theme.R")
  
  # Load data
  sapply(paste0("Fun/graphs/",list.files("Fun/graphs")),source)
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files",project,"/H1.all.R"))
  
  for(sim in covariateSize){

    generalsubtitle = paste0("Among 100 simulated datasets of ",dataType," for ",nInd," individuals, with ",sim," covariates.\n")
    Titlelist = list(reg = "Model built with SAMBA.",
                     lasso = "Model built with a lasso approach within SAMBA, whithout stability selection and a 100 search scoope.",
                     lassoSS = "Model built with a lasso approach within SAMBA.",
                     elasticnetSS = "Model built with an elastic net approach within SAMBA.",
                     lassoSSCrit = "Model built with a lasso approach within SAMBA, and multiple thresholds.",
                     elasticnetSSCrit ="Model built with an elastic net approach within SAMBA, and multiple thresholds")
    
    initFolder=paste0("~/Travail/Presentation/Plot",project)
    if(!dir.exists(initFolder)){dir.create(initFolder)}
    initFolder=paste0(initFolder,"/Results",sim)
    if(!dir.exists(initFolder)){dir.create(initFolder)}
    
    
    
    # Graphs summarizing building result
    for(meth in buildMethod){
      Folder=paste0(initFolder,"/",meth)
      if(!dir.exists(Folder)){dir.create(Folder)}
      subtitle = paste0(generalsubtitle,Titlelist[[meth]])
      graphsTotalNB(Folder,subtitle,project,sim,meth,JPEG,PNG)
      graphsParNB(Folder,subtitle,project,sim,meth,JPEG,PNG)
      tableStats(Folder,subtitle,project,sim,meth,JPEG,PNG)
    }
    
    # Comparison graphs
    if(length(buildMethod)>1){
      Folder=initFolder
      subtitle = generalsubtitle
      graphsCompMethod(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      graphsParCompMethod(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      graphsCompFR(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      tableStatsComp(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
      graphsCompTime(Folder,subtitle,project,sim,buildMethod,JPEG,PNG)
    }
  }
}

