suppressWarnings(suppressMessages(library(ggplot2,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggpubr,quietly=TRUE)))
suppressWarnings(suppressMessages(library(webshot2,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggh4x,quietly=TRUE)))
suppressWarnings(suppressMessages(library(flextable,quietly=TRUE)))
suppressWarnings(suppressMessages(library(grid,quietly=TRUE)))
suppressWarnings(suppressMessages(library(gtable,quietly=TRUE)))
suppressWarnings(suppressMessages(library(gridExtra,quietly=TRUE)))
suppressWarnings(suppressMessages(library(scales,quietly=TRUE)))

graphsGenerate <- function(project="Pasin",
                           covariateSize=c(10,50,200,500),
                           buildMethod=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"),
                           JPEG = FALSE,
                           PNG = TRUE){
  eval(parse(text=readLines(paste0("Files/Files",project,"/info.txt"))))
  
  source("~/Travail/00_Theme.R")
  
  # Load data
  sapply(paste0("Fun/graphs/",list.files("Fun/graphs")),source)
  
  for(sim in covariateSize){

    generalsubtitle = paste0("Among 100 simulated datasets of ",dataType," for ",nInd," individuals, with ",sim," covariates.\n")
    Titlelist = list(reg = "Model built with SAMBA.",
                     lasso = "Model built with a lasso approach within SAMBA, whithout stability selection and a 100 search scoope.",
                     lassoSS = "Model built with a lasso approach within SAMBA.",
                     elasticnetSS = "Model built with an elastic net approach within SAMBA.",
                     lassoSSCrit = "Model built with a lasso approach within SAMBA, and multiple thresholds.",
                     elasticnetSSCrit ="Model built with an elastic net approach within SAMBA, and multiple thresholds")
    
    initFolder=paste0("~/Travail/Presentation/Plot/LassoSAMBA/Plot",project)
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

